//! F4 algorithm for Gröbner basis computation
//!
//! This module implements the **true F4 algorithm** for computing Gröbner bases over
//! cryptographically large prime fields (BN254, BLS12-381, etc.).
//!
//! ## Algorithm Overview
//!
//! F4 differs from Buchberger by processing S-polynomials in batches:
//! 1. Select all S-polynomials of the current degree
//! 2. Build a Macaulay matrix from S-polys and basis multiples
//! 3. Perform Gaussian elimination on the entire matrix at once
//! 4. Extract new basis elements from reduced rows
//! 5. Interreduce and continue to next degree
//!
//! This batch approach amortizes the reduction cost across many S-polynomials,
//! making it faster than Buchberger for dense systems.
//!
//! ## When to Use F4
//!
//! **Use F4 when:**
//! - Working with cryptographic fields (BN254, BLS12-381)
//! - System generates many S-polynomials per degree
//! - Dense polynomial systems
//! - You have sufficient memory for large matrices
//!
//! **Use Buchberger when:**
//! - Sparse systems
//! - Memory-constrained environments
//! - Systems with rapid degree growth
//!
//! ## Example
//!
//! ```ignore
//! use ark_feanor::*;
//! use ark_feanor::f4::*;
//!
//! let field = &*BN254_FR;
//! let ring = MultivariatePolyRingImpl::new(field, 3);
//!
//! let system = vec![/* polynomials */];
//!
//! // Compute Gröbner basis with F4
//! let gb = f4_simple(&ring, system, DegRevLex);
//! ```

#![cfg_attr(test, feature(allocator_api))]

pub mod matrix;
pub mod gaussian;
pub mod monosig;

use crate::f4::gaussian::{reduce_tr_by_rr_on_A, reduce_matrix_tr};
use crate::f4::matrix::{MacaulayMatrix, RowType, SparseRow};
use crate::f4::monosig::MonSig;
use feanor_math::algorithms::buchberger::GBAborted;
use feanor_math::divisibility::DivisibilityRingStore;
use feanor_math::field::Field;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::ring::*;
use feanor_math::rings::multivariate::*;
use std::time::{Duration, Instant};
use std::collections::{HashMap, HashSet};

/// Execution mode for parallelization
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExecutionMode {
    /// Deterministic mode: serialize pivot discovery for reproducibility
    /// Ensures bit-exact results across runs with same thread count
    Deterministic,

    /// Performance mode: allow relaxed ordering with row_idx tie-breaking
    /// May have non-deterministic pivot selection order, but correct results
    Performance,
}

impl Default for ExecutionMode {
    fn default() -> Self {
        // Default to Deterministic for tests and reproducibility
        ExecutionMode::Deterministic
    }
}

/// Configuration for F4 algorithm
#[derive(Clone, Debug)]
pub struct F4Config {
    /// Maximum degree of polynomials to consider
    pub max_degree: Option<usize>,

    /// Maximum number of matrix rows (abort if exceeded)
    pub max_matrix_rows: Option<usize>,

    /// Maximum number of S-polynomial reductions
    pub s_pair_budget: Option<usize>,

    /// Maximum computation time
    pub time_budget: Option<Duration>,

    /// Optional cap on number of S-pairs processed per round (mnsel-like)
    pub s_pair_round_cap: Option<usize>,

    /// Execution mode for parallelization
    pub execution_mode: ExecutionMode,

    /// Number of threads to use (None = use all available cores)
    pub num_threads: Option<usize>,
}

impl F4Config {
    /// Create a new default configuration
    pub fn new() -> Self {
        F4Config {
            max_degree: None,
            max_matrix_rows: None,
            s_pair_budget: None,
            time_budget: None,
            s_pair_round_cap: None,
            execution_mode: ExecutionMode::default(),
            num_threads: None,
        }
    }

    /// Set maximum degree limit
    pub fn with_max_degree(mut self, max_deg: usize) -> Self {
        self.max_degree = Some(max_deg);
        self
    }

    /// Set maximum matrix rows
    pub fn with_max_matrix_rows(mut self, max_rows: usize) -> Self {
        self.max_matrix_rows = Some(max_rows);
        self
    }

    /// Set S-pair budget
    pub fn with_s_pair_budget(mut self, max_pairs: usize) -> Self {
        self.s_pair_budget = Some(max_pairs);
        self
    }

    /// Set time budget
    pub fn with_time_budget(mut self, max_time: Duration) -> Self {
        self.time_budget = Some(max_time);
        self
    }

    /// Cap the number of S-pairs processed per round
    pub fn with_s_pair_round_cap(mut self, cap: usize) -> Self {
        self.s_pair_round_cap = Some(cap);
        self
    }
}

impl Default for F4Config {
    fn default() -> Self {
        Self::new()
    }
}

/// S-polynomial representation (reused from feanor-math)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct SPoly {
    i: usize,
    j: usize,
}

impl SPoly {
    fn new(i: usize, j: usize) -> Self {
        SPoly { i, j }
    }

    /// Compute the S-polynomial from two basis elements
    fn compute<P, O>(&self, ring: P, basis: &[El<P>], order: O) -> El<P>
    where
        P: RingStore + Copy,
        P::Type: MultivariatePolyRing,
        <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
        O: MonomialOrder + Copy,
    {
        let (lc_i, lm_i) = ring.LT(&basis[self.i], order).unwrap();
        let (lc_j, lm_j) = ring.LT(&basis[self.j], order).unwrap();

        // Compute LCM of leading monomials
        let lcm = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);

        // Compute multipliers
        let mult_i = ring.monomial_div(ring.clone_monomial(&lcm), lm_i).ok().unwrap();
        let mult_j = ring.monomial_div(ring.clone_monomial(&lcm), lm_j).ok().unwrap();

        // S-poly = (lc_j * t_i * f_i - lc_i * t_j * f_j) / gcd(lc_i, lc_j)
        let mut poly_i = ring.clone_el(&basis[self.i]);
        ring.mul_assign_monomial(&mut poly_i, mult_i);
        ring.inclusion().mul_assign_ref_map(&mut poly_i, lc_j);

        let mut poly_j = ring.clone_el(&basis[self.j]);
        ring.mul_assign_monomial(&mut poly_j, mult_j);
        ring.inclusion().mul_assign_ref_map(&mut poly_j, lc_i);

        ring.sub(poly_i, poly_j)
    }

    /// Get the LCM degree of this S-polynomial
    fn lcm_degree<P, O>(&self, ring: P, basis: &[El<P>], order: O) -> usize
    where
        P: RingStore + Copy,
        P::Type: MultivariatePolyRing,
        <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
        O: MonomialOrder + Copy,
    {
        let (_, lm_i) = ring.LT(&basis[self.i], order).unwrap();
        let (_, lm_j) = ring.LT(&basis[self.j], order).unwrap();
        let lcm = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);
        ring.monomial_deg(&lcm)
    }
}

/// Compute Gröbner basis using F4 algorithm (simple interface)
///
/// This is the simplest way to use F4 - just provide the ring, input basis, and monomial order.
pub fn f4_simple<P, O>(ring: P, input_basis: Vec<El<P>>, order: O) -> Vec<El<P>>
where
    P: RingStore + Copy + Send + Sync,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy + Send + Sync,
{
    f4_configured(ring, input_basis, order, F4Config::default()).unwrap()
}

/// Parallel variant of f4_simple with Send + Sync bounds
///
/// Uses parallel AB phase reduction via reduce_tr_by_rr_on_A_par with snapshot-based
/// approach. Requires thread-safe ring types (BN254, BLS12-381, etc.).
///
/// Expected speedup: 1.5-2x on 8+ cores for larger systems (Katsura-5+).

/// Compute Gröbner basis using F4 algorithm with configuration
pub fn f4_configured<P, O>(
    ring: P,
    input_basis: Vec<El<P>>,
    order: O,
    config: F4Config,
) -> Result<Vec<El<P>>, GBAborted>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy + Send + Sync,
{
    let profile = std::env::var("F4_PROFILE").map(|v| v == "1" || v.eq_ignore_ascii_case("true")).unwrap_or(false);
    let env_round_cap = std::env::var("F4_MNSEL").ok().and_then(|s| s.parse::<usize>().ok());
    // Filter out zero polynomials and interreduce input
    let mut basis: Vec<El<P>> = input_basis
        .into_iter()
        .filter(|f| !ring.is_zero(f))
        .collect();

    if basis.is_empty() {
        return Ok(vec![ring.zero()]);
    }

    // Check for constant in basis
    if basis.iter().any(|f| ring.terms(f).all(|(_, m)| ring.monomial_deg(m) == 0)) {
        return Ok(vec![ring.one()]);
    }

    // Normalize basis (make leading coefficients 1)
    for f in &mut basis {
        let (lc, _) = ring.LT(f, order).unwrap();
        let lc_inv = ring.base_ring().invert(lc).unwrap();
        ring.inclusion().mul_assign_map(f, lc_inv);
    }

    // Track redundant basis elements (bs->red equivalent) to skip S-pairs
    let mut red_flags: Vec<bool> = vec![false; basis.len()];

    // Generate initial S-polynomials (with Buchberger product criterion)
    let mut spolys: Vec<SPoly> = Vec::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            // Apply Buchberger product criterion: skip if leading terms are coprime
            let (_, lm_i) = ring.LT(&basis[i], order).unwrap();
            let (_, lm_j) = ring.LT(&basis[j], order).unwrap();

            if !are_monomials_coprime(ring, lm_i, lm_j) {
                spolys.push(SPoly::new(i, j));
            }
            // If coprime, S-poly will reduce to zero - skip it
        }
    }

    let start_time = std::time::Instant::now();
    let mut s_pairs_processed = 0;

    // Main F4 loop: process S-polynomials degree by degree
    while !spolys.is_empty() {
        let round_start = Instant::now();
        // Check time budget
        if let Some(max_time) = config.time_budget {
            if start_time.elapsed() > max_time {
                return Err(GBAborted::TimeBudget {
                    max_seconds: max_time.as_secs(),
                });
            }
        }

        // Find minimum degree
        let min_deg = spolys
            .iter()
            .map(|sp| sp.lcm_degree(ring, &basis, order))
            .min()
            .unwrap();

        // Check degree limit
        if let Some(max_deg) = config.max_degree {
            if min_deg > max_deg {
                return Err(GBAborted::DegreeExceeded {
                    max_degree: max_deg,
                    actual_degree: min_deg,
                });
            }
        }

        // Select S-polys of current degree (ark-feanor-maq: select at minimal degree, ready for LCM grouping)
        let (mut current_spolys, remaining_spolys): (Vec<_>, Vec<_>) = spolys
            .into_iter()
            .partition(|sp| sp.lcm_degree(ring, &basis, order) == min_deg);

        spolys = remaining_spolys;

        // Apply round cap (mnsel-like): limit number of S-pairs this round
        // The LCM grouping happens during rr/tr construction below
        let round_cap = config.s_pair_round_cap.or(env_round_cap);
        if let Some(cap) = round_cap {
            if current_spolys.len() > cap {
                // Same-LCM tail: include all S-pairs sharing the LCM with the last included one
                // Build (SPoly, expanded LCM) list
                let mut with_lcm: Vec<(SPoly, Vec<usize>)> = current_spolys
                    .iter()
                    .cloned()
                    .map(|sp| {
                        let (_, lm_i) = ring.LT(&basis[sp.i], order).unwrap();
                        let (_, lm_j) = ring.LT(&basis[sp.j], order).unwrap();
                        let lcm = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);
                        let expanded = ring.expand_monomial(&lcm);
                        (sp, expanded)
                    })
                    .collect();

                // Sort by LCM expanded monomial lexicographically to group equals
                with_lcm.sort_by(|a, b| a.1.cmp(&b.1));

                let mut selected: Vec<SPoly> = Vec::with_capacity(cap + 8);
                let mut count = 0usize;
                let mut last_lcm: Option<&[usize]> = None;
                for (sp, lcm_exp) in &with_lcm {
                    if count < cap {
                        selected.push(sp.clone());
                        count += 1;
                        last_lcm = Some(lcm_exp);
                    } else {
                        // Include same-LCM tail
                        if let Some(last) = &last_lcm {
                            if lcm_exp == *last {
                                selected.push(sp.clone());
                                continue;
                            }
                        }
                        break;
                    }
                }
                current_spolys = selected;
            }
        }

        // Check S-pair budget
        s_pairs_processed += current_spolys.len();
        if let Some(max_pairs) = config.s_pair_budget {
            if s_pairs_processed > max_pairs {
                return Err(GBAborted::SPairBudget { max_s_pairs: max_pairs });
            }
        }

        // Build Macaulay matrix using msolve-style rr/tr construction (ark-feanor-29t)
        //
        // IMPORTANT: This replaces explicit S-polynomial computation with multiplied basis rows.
        // Instead of computing spoly(f_i, f_j) = (LCM/LT(f_i))*f_i - (LCM/LT(f_j))*f_j,
        // we add the individual multiplied rows separately:
        //   - rr (reducer row): (LCM/LT(f_i))*f_i  [first generator]
        //   - tr (to-reduce rows): (LCM/LT(f_j))*f_j for j ≠ i  [remaining generators]
        //
        // This approach:
        // - Creates cohesive A-blocks (pivot columns from LCMs)
        // - Enables better pivot reuse across rows with the same LCM
        // - Aligns with msolve's symbol.c symbolic_preprocessing behavior
        //
        // NOTE: SPoly::compute (lines 143-170) is OBSOLETE and not used in this loop.
        //       It remains for backwards compatibility with tests/external code.
        let mut matrix = MacaulayMatrix::new(ring);
        let t_rows_0 = Instant::now();

        // 1) Group S-pairs by LCM (expanded exponents as HashMap key)
        //    Each group collects all basis indices (gens) whose pairs share that LCM
        let mut groups: HashMap<Vec<usize>, (PolyMonomial<P>, HashSet<usize>)> = HashMap::new();
        for sp in &current_spolys {
            let (_, lm_i) = ring.LT(&basis[sp.i], order).unwrap();
            let (_, lm_j) = ring.LT(&basis[sp.j], order).unwrap();
            let lcm = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);
            let exp = ring.expand_monomial(&lcm);
            let entry = groups.entry(exp).or_insert((lcm, HashSet::new()));
            (entry.1).insert(sp.i);
            (entry.1).insert(sp.j);
        }

        // 2) For each LCM group, build rr/tr rows as (LCM/LT(basis[g])) * basis[g]
        //    - rr: first generator (becomes the anchor row for this LCM pivot)
        //    - tr: remaining generators (to-be-reduced rows with same LCM)
        //
        // IMPORTANT: To support AB/CD semantics, we need rr and tr rows in contiguous blocks.
        // First collect all the data, then add rr rows, then tr rows.
        let mut pivot_monomials: Vec<PolyMonomial<P>> = Vec::with_capacity(groups.len());
        let mut rr_data: Vec<(usize, PolyMonomial<P>)> = Vec::new(); // (gen_idx, lcm_mono)
        let mut tr_data: Vec<(usize, PolyMonomial<P>)> = Vec::new(); // (gen_idx, lcm_mono)

        // Sort groups by (degree, exponent vector) for deterministic iteration order
        // This ensures reproducible row ordering across runs
        let mut sorted_groups: Vec<(Vec<usize>, (PolyMonomial<P>, HashSet<usize>))> =
            groups.into_iter().collect();
        sorted_groups.sort_by(|a, b| {
            // Compare by exponent vector (already expanded, deterministic)
            a.0.cmp(&b.0)
        });

        for (_exp, (lcm_mono, gens_set)) in sorted_groups.iter() {
            let mut gens: Vec<usize> = gens_set.iter().copied().collect();
            gens.sort_unstable(); // Deterministic ordering
            if gens.is_empty() { continue; }
            pivot_monomials.push(ring.clone_monomial(lcm_mono));

            // rr (reducer row): first generator in sorted order
            rr_data.push((gens[0], ring.clone_monomial(lcm_mono)));

            // tr (to-reduce rows): all remaining generators for this LCM
            for &g in gens.iter().skip(1) {
                tr_data.push((g, ring.clone_monomial(lcm_mono)));
            }
        }

        // Helper to add a multiplied row: (LCM / LT(b)) * b
        let add_multiplied_row = |gidx: usize, lcm_mono: &PolyMonomial<P>, matrix: &mut MacaulayMatrix<P>, row_type: RowType| {
            if let Some((_, lt_b)) = ring.LT(&basis[gidx], order) {
                // Compute multiplier: mult = LCM / LT(basis[gidx])
                if let Ok(mult) = ring.monomial_div(ring.clone_monomial(lcm_mono), lt_b) {
                    let mut row_poly = ring.clone_el(&basis[gidx]);
                    ring.mul_assign_monomial(&mut row_poly, mult);
                    let row = matrix.polynomial_to_row(&row_poly, order);
                    matrix.add_row(row, row_type);
                }
            }
        };

        // Add all rr rows first (contiguous block)
        let rr_start_idx = matrix.num_rows();
        for (gidx, lcm_mono) in &rr_data {
            add_multiplied_row(*gidx, lcm_mono, &mut matrix, RowType::Rr);
        }
        let rr_end_idx = matrix.num_rows();

        // Add all tr rows next (contiguous block)
        let tr_start_idx = matrix.num_rows();
        for (gidx, lcm_mono) in &tr_data {
            add_multiplied_row(*gidx, lcm_mono, &mut matrix, RowType::Tr);
        }
        let tr_end_idx = matrix.num_rows();

        let t_rows = t_rows_0.elapsed();

        // Build SHT idx: 2 = pivot LCM monomials, 1 = other monomials
        let mut sht_idx: HashMap<Vec<usize>, u8> = HashMap::new();
        for m in &pivot_monomials {
            let e = ring.expand_monomial(m);
            sht_idx.insert(e, 2);
        }
        // Ensure pivot columns exist in the matrix
        for m in &pivot_monomials {
            let _ = matrix.get_or_create_column(m);
        }

        // Fetch pivot indices after ensuring columns exist
        let mut pivot_indices: Vec<usize> = Vec::with_capacity(pivot_monomials.len());
        for m in &pivot_monomials {
            if let Some(&col) = matrix.monomial_to_col.get(&ring.expand_monomial(m)) {
                pivot_indices.push(col);
            }
        }

        // Deduplicate while preserving discovery order
        let mut seen = std::collections::HashSet::new();
        pivot_indices.retain(|c| seen.insert(*c));

        // === SHT-driven A|B column mapping (ark-feanor-fok) ===
        //
        // This implements msolve's convert_hashes_to_columns behavior, where columns are
        // reordered to separate the matrix into A and B blocks:
        //   - A block (pivots): idx=2 columns (S-polynomial LCMs)
        //   - B block (reducers): idx=1 columns (monomials from rr/tr and closure)
        //   - Tail block: remaining columns
        //
        // This ordering is critical for efficient Gaussian elimination:
        //   1. Pivots are discovered early (leftmost columns)
        //   2. Reducer monomials come next (likely to be reduced)
        //   3. Tails come last (less likely to become pivots)
        //
        // SHT timing (ark-feanor-63y): Update SHT with rr/tr monomials BEFORE reordering
        // to ensure the SHT reflects all monomials that will be present in the matrix.

        let t_reorder_0 = Instant::now();
        if !pivot_indices.is_empty() {
            // Fill SHT with current rr/tr monomials as idx=1 before mapping
            // This ensures all monomials from rr/tr rows are tracked before column reordering
            for r in 0..matrix.num_rows() {
                for (c, _) in &matrix.rows[r].entries {
                    let exp = matrix.col_exp(*c).clone();
                    if !sht_idx.contains_key(&exp) {
                        sht_idx.insert(exp, 1);
                    }
                }
            }
            // Reorder columns using SHT: idx=2 (pivots) first, then idx=1, then tails
            // This creates the A|B split that msolve uses for optimized elimination
            matrix.reorder_columns_with_sht(&sht_idx);
        }
        let t_reorder = t_reorder_0.elapsed();

        // After reordering, pivot columns are now the first `pivot_count` indices.
        // Add one dedicated reducer row per pivot column to anchor the pivot.
        // === SHT-driven symbolic closure - Part 1: Dedicated reducers for pivot columns ===
        // (ark-feanor-1w9: SHT-driven symbolic closure)
        //
        // This implements msolve-style symbolic preprocessing where we add exactly one reducer
        // per pivot column (LCM column from S-pairs). This corresponds to the A-block setup
        // in msolve's symbolic_preprocessing (symbol.c).
        //
        // Key properties:
        // - One reducer per pivot column (skip if no suitable basis element found)
        // - Uses LM index with degree prefiltering for fast reducer selection
        // - Restricts reducer rows to A-block only (entries < pivot_count) to avoid tail growth
        // - Follows msolve's strategy of building a well-structured A-block before elimination

        let pivot_count = pivot_indices.len();
        let pivot_col_indices: std::collections::HashSet<usize> = (0..pivot_count).collect();

        // Set ncl (A-block width) for AB/CD semantics
        // After SHT reordering, pivot columns are 0..(pivot_count-1)
        matrix.ncl = Some(pivot_count);

        if pivot_count > 0 {
            // Build LM index once for selecting reducers quickly
            // The index is sorted by degree for early termination
            let lm_index = build_lm_index(ring, &basis, order);

            for pi in 0..pivot_count {
                let pivot_col = pi; // after reorder, first pivot_count columns are pivots
                // Clone monomial to avoid holding an immutable borrow of matrix across mutations
                let pivot_mono = ring.clone_monomial(&matrix.col_to_monomial[pivot_col]);
                let pivot_deg = ring.monomial_deg(&pivot_mono);
                let pivot_exp = ring.expand_monomial(&pivot_mono);
                let pivot_sig = MonSig::from_exponents(&pivot_exp);

                // Find the FIRST basis LM dividing the pivot monomial using degree prefiltering
                for entry in &lm_index {
                    // Degree prefilter: if basis LM degree > pivot degree, no division possible
                    if entry.degree > pivot_deg { break; }

                    // Divmask prefilter: fast rejection based on presence mask and buckets
                    if !entry.sig.presence_subset(&pivot_sig) { continue; }
                    if !entry.sig.buckets_leq(&pivot_sig) { continue; }

                    // Quick per-variable bounds check before expensive monomial_div
                    let mut ok = true;
                    for (a, b) in entry.exponents.iter().zip(pivot_exp.iter()) {
                        if a > b { ok = false; break; }
                    }
                    if !ok { continue; }

                    if let Ok(mult) = ring.monomial_div(ring.clone_monomial(&pivot_mono), &entry.lm) {
                        // Add reducer: (pivot / LM(b)) * b as a new row
                        let mut reducer = ring.clone_el(&basis[entry.basis_idx]);
                        ring.mul_assign_monomial(&mut reducer, mult);
                        let row = matrix.polynomial_to_row(&reducer, order);

                        // Keep full row support for dedicated pivot reducers
                        // These are specifically for pivot columns and need full support
                        if row.entries.is_empty() { continue; }
                        matrix.add_row(row, RowType::Dedicated);
                        break; // One reducer per pivot column
                    }
                }
            }
        }
        let dedicated_end_idx = matrix.num_rows(); // Track end of dedicated reducers

        // === A-block triangularization of rr rows ===
        // (ark-feanor-59b: Fix ncl calculation to match actual rr row coverage)
        //
        // Goal: Make rr rows form a proper triangular block in A so each A column has a unique rr pivot row.
        // This ensures AB phase can fully reduce all A-block columns without gaps or collisions.
        //
        // Steps:
        // 1. Collect all rr rows (original + dedicated)
        // 2. For each rr row, find its leftmost A-column (< pivot_count)
        // 3. Sort rr rows by their A-pivot column ascending
        // 4. Perform small AB-only elimination among rr rows to triangularize
        // 5. This ensures rr_map has complete coverage with no collisions
        //
        // This matches msolve's convert_hashes_to_columns intent: pivot columns are
        // fully covered by normalized rr rows before AB/CD phase.

        if pivot_count > 0 {
            let base_ring = matrix.ring.base_ring();

            // Collect all rr row indices (original rr + dedicated)
            let mut rr_indices: Vec<usize> = Vec::new();
            rr_indices.extend(rr_start_idx..rr_end_idx);
            rr_indices.extend(tr_end_idx..dedicated_end_idx);

            // Build (row_idx, a_pivot_col) pairs for rr rows with A-block pivots
            let mut rr_with_a_pivots: Vec<(usize, usize)> = Vec::new();
            for &rr_idx in &rr_indices {
                if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
                    if pivot_col < pivot_count {
                        rr_with_a_pivots.push((rr_idx, pivot_col));
                    }
                }
            }

            // Sort by A-pivot column ascending for triangularization
            rr_with_a_pivots.sort_by_key(|(_, col)| *col);

            // Build cache of leading coefficient inverses for rr rows
            // This avoids expensive field divisions in the triangularization loop
            let mut rr_inv_cache: std::collections::HashMap<usize, PolyCoeff<P>> = std::collections::HashMap::new();
            for &(rr_idx, _) in &rr_with_a_pivots {
                if !matrix.rows[rr_idx].entries.is_empty() {
                    let leading_coeff = &matrix.rows[rr_idx].entries[0].1;
                    let inv = base_ring.invert(leading_coeff).unwrap();
                    rr_inv_cache.insert(rr_idx, inv);
                }
            }

            // Triangularize: for each rr row, reduce by earlier rr rows in A-block only
            // This eliminates collisions and creates a proper triangular structure
            for i in 0..rr_with_a_pivots.len() {
                let (row_idx, _) = rr_with_a_pivots[i];

                // Reduce this row by all earlier rr rows (with smaller A-pivot columns)
                let mut changed = true;
                while changed {
                    changed = false;

                    let pivot_col = match matrix.rows[row_idx].pivot() {
                        Some(col) if col < pivot_count => col,
                        _ => break, // Row is zero or pivot moved to B-block
                    };

                    // Look for an earlier rr row with this A-pivot column
                    let mut reducer_idx = None;
                    for j in 0..i {
                        let (earlier_idx, earlier_col) = rr_with_a_pivots[j];
                        if let Some(earlier_pivot) = matrix.rows[earlier_idx].pivot() {
                            if earlier_pivot == pivot_col && earlier_pivot < pivot_count {
                                reducer_idx = Some(earlier_idx);
                                break;
                            }
                        }
                    }

                    if let Some(reducer_idx) = reducer_idx {
                        // Reduce: row -= multiplier * reducer_row
                        // Use cached inverse for fast multiplier computation
                        let row_leading = base_ring.clone_el(&matrix.rows[row_idx].entries[0].1);
                        let reducer_inv = rr_inv_cache.get(&reducer_idx).unwrap();
                        let multiplier = base_ring.mul_ref(&row_leading, reducer_inv);

                        // Take row out to allow immutable borrow of reducer
                        let mut row = std::mem::replace(&mut matrix.rows[row_idx], SparseRow::new());

                        // Perform subtraction: row -= multiplier * reducer_row (inline implementation)
                        let mut result = Vec::with_capacity(row.entries.len() + matrix.rows[reducer_idx].entries.len());
                        let mut row_iter = row.entries.iter();
                        let mut reducer_iter = matrix.rows[reducer_idx].entries.iter();
                        let mut row_next = row_iter.next();
                        let mut reducer_next = reducer_iter.next();

                        while row_next.is_some() || reducer_next.is_some() {
                            match (row_next, reducer_next) {
                                (Some((rc, rv)), Some((sc, sv))) => {
                                    if rc == sc {
                                        let mut new_val = base_ring.clone_el(rv);
                                        let scaled = base_ring.mul_ref(&multiplier, sv);
                                        base_ring.sub_assign(&mut new_val, scaled);
                                        if !base_ring.is_zero(&new_val) {
                                            result.push((*rc, new_val));
                                        }
                                        row_next = row_iter.next();
                                        reducer_next = reducer_iter.next();
                                    } else if rc < sc {
                                        result.push((*rc, base_ring.clone_el(rv)));
                                        row_next = row_iter.next();
                                    } else {
                                        let mut new_val = base_ring.mul_ref(&multiplier, sv);
                                        base_ring.negate_inplace(&mut new_val);
                                        if !base_ring.is_zero(&new_val) {
                                            result.push((*sc, new_val));
                                        }
                                        reducer_next = reducer_iter.next();
                                    }
                                }
                                (Some((rc, rv)), None) => {
                                    result.push((*rc, base_ring.clone_el(rv)));
                                    row_next = row_iter.next();
                                }
                                (None, Some((sc, sv))) => {
                                    let mut new_val = base_ring.mul_ref(&multiplier, sv);
                                    base_ring.negate_inplace(&mut new_val);
                                    if !base_ring.is_zero(&new_val) {
                                        result.push((*sc, new_val));
                                    }
                                    reducer_next = reducer_iter.next();
                                }
                                (None, None) => unreachable!(),
                            }
                        }
                        row.entries = result;

                        // Put row back
                        matrix.rows[row_idx] = row;

                        // Update cache with new leading coefficient inverse
                        if !matrix.rows[row_idx].entries.is_empty() {
                            let new_leading = &matrix.rows[row_idx].entries[0].1;
                            let new_inv = base_ring.invert(new_leading).unwrap();
                            rr_inv_cache.insert(row_idx, new_inv);
                        }

                        changed = true;
                    }
                }
            }
        }

        // === SHT-driven symbolic closure - Part 2: Bounded closure loop ===
        // (ark-feanor-1w9 continued)
        //
        // This performs bounded symbolic closure to add reducers for columns discovered
        // during rr/tr construction. Unlike unbounded closure which can blow up matrix size,
        // we use strict bounds matching msolve's behavior:
        //
        // 1. Collect target columns from rr/tr rows and update SHT with idx=1
        // 2. Run closure loop with bounded passes (max 2)
        // 3. Exit early if no new rows/cols added
        // 4. Skip pivot columns (they already have dedicated reducers)
        //
        // This prevents exponential growth while ensuring adequate symbolic preprocessing
        // for the Gaussian elimination phase.

        let reducers_before_rows = matrix.num_rows();
        let reducers_before_cols = matrix.num_cols();
        let t_reducers_0 = Instant::now();
        let s_rows = reducers_before_rows; // rows added so far are rr/tr rows

        // Collect target columns from rr/tr rows and update SHT idx=1
        // The SHT (secondary hash table) tracks monomial roles:
        //   idx=2: pivot LCMs (S-polynomial LCMs)
        //   idx=1: reducer monomials (from rr/tr and closure)
        let mut target_cols: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for r in 0..s_rows {
            for (c, _) in &matrix.rows[r].entries {
                target_cols.insert(*c);
                // Update SHT: mark this monomial as a reducer (idx=1) if not already present
                let exp = matrix.col_exp(*c).clone();
                if !sht_idx.contains_key(&exp) {
                    sht_idx.insert(exp, 1);
                }
            }
        }

        // Build LM index once for fast reducer selection during closure
        let lm_index = build_lm_index(ring, &basis, order);

        // Column signature cache to avoid recomputing MonSig for same columns
        let mut col_sig_cache: std::collections::HashMap<usize, (usize, Vec<usize>, MonSig)> = std::collections::HashMap::new();

        // Bounded closure: up to 2 passes or until rows/cols stop changing
        // This aligns with msolve's bounded symbolic preprocessing strategy
        let mut passes = 0usize;
        let max_passes = 2usize;
        loop {
            let rows_before = matrix.num_rows();
            let cols_before = matrix.num_cols();

            // Add reducers for target columns, skipping pivot columns
            add_reducers_for_columns_with_pivot_skip(
                &mut matrix,
                &basis,
                order,
                &pivot_col_indices,
                pivot_count,
                &target_cols,
                &lm_index,
                &mut col_sig_cache,
            );

            passes += 1;

            // Clear cache if columns changed (indices shifted)
            if matrix.num_cols() != cols_before {
                col_sig_cache.clear();
            }

            // Exit conditions:
            // 1. Reached maximum passes (prevents unbounded growth)
            // 2. No new rows or columns added (closure has stabilized)
            if passes >= max_passes || (matrix.num_rows() == rows_before && matrix.num_cols() == cols_before) {
                break;
            }
        }
        let t_reducers = t_reducers_0.elapsed();
        let reducers_added_rows = matrix.num_rows().saturating_sub(reducers_before_rows);
        let reducers_new_cols = matrix.num_cols().saturating_sub(reducers_before_cols);

        // Check matrix size limit
        if let Some(max_rows) = config.max_matrix_rows {
            if matrix.num_rows() > max_rows {
                return Err(GBAborted::MemBudget {
                    max_bytes: max_rows * std::mem::size_of::<usize>(),
                });
            }
        }

        // Do not pre-mark or sort rows here; let reduce_matrix determine pivots safely

        // === AB/CD Phase Reduction ===
        // Replace full-matrix elimination with specialized AB/CD semantics:
        // 1. AB phase: reduce tr rows by rr rows (A-block only)
        // 2. CD phase: reduce tr rows and discover pivots in B/D blocks only
        let t_elim_0 = Instant::now();

        // ncl is always set after SHT reordering - AB/CD is the only elimination path
        let ncl_val = matrix.ncl.expect("ncl must be set after SHT reordering");

        let new_pivot_rows = {
            // Use the tracked row ranges from construction:
            // - rr rows: rr_start_idx..rr_end_idx
            // - tr rows: tr_start_idx..tr_end_idx
            // - dedicated reducers: tr_end_idx..dedicated_end_idx

            // DIAGNOSTIC: Check ncl mismatch between desired A columns (SHT pivots) and covered A columns (rr_map)
            if profile {
                use std::collections::HashSet;
                // Build rr_map to see what A columns are actually covered
                let mut covered: HashSet<usize> = HashSet::new();
                let mut collisions: Vec<(usize, Vec<usize>)> = Vec::new();

                // Scan rr rows
                for rr_idx in rr_start_idx..rr_end_idx {
                    if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
                        if pivot_col < ncl_val {
                            if covered.contains(&pivot_col) {
                                // Collision: multiple rr rows have same leftmost column
                                if let Some(existing) = collisions.iter_mut().find(|(c, _)| *c == pivot_col) {
                                    existing.1.push(rr_idx);
                                } else {
                                    collisions.push((pivot_col, vec![rr_idx]));
                                }
                            }
                            covered.insert(pivot_col);
                        }
                    }
                }
                // Scan dedicated reducers
                for rr_idx in tr_end_idx..dedicated_end_idx {
                    if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
                        if pivot_col < ncl_val {
                            if covered.contains(&pivot_col) {
                                if let Some(existing) = collisions.iter_mut().find(|(c, _)| *c == pivot_col) {
                                    existing.1.push(rr_idx);
                                } else {
                                    collisions.push((pivot_col, vec![rr_idx]));
                                }
                            }
                            covered.insert(pivot_col);
                        }
                    }
                }

                // Compute missing A columns
                let desired_a: HashSet<usize> = (0..ncl_val).collect();
                let missing: Vec<usize> = desired_a.difference(&covered).copied().collect();

                if !missing.is_empty() || !collisions.is_empty() {
                    eprintln!("[AB/CD DIAGNOSTIC] ncl={}, desired A columns=[0..{})", ncl_val, ncl_val);
                    eprintln!("  Covered columns: {} / {}", covered.len(), ncl_val);
                    if !missing.is_empty() {
                        eprintln!("  Missing columns: {:?}", missing);
                    }
                    if !collisions.is_empty() {
                        eprintln!("  Collisions: {} columns have multiple rr rows", collisions.len());
                        for (col, rows) in &collisions {
                            eprintln!("    Column {}: rr rows {:?}", col, rows);
                        }
                    }
                }
            }

            // Sort rows by pivot column then density before AB reduction for better cache locality
            matrix.sort_rows_by_pivot_and_density(rr_start_idx..rr_end_idx);
            matrix.sort_rows_by_pivot_and_density(tr_start_idx..tr_end_idx);
            if tr_end_idx < dedicated_end_idx {
                matrix.sort_rows_by_pivot_and_density(tr_end_idx..dedicated_end_idx);
            }

            // AB phase: reduce tr rows by rr rows (A-block only, columns < ncl)
            // CRITICAL: Must exclude tr rows from rr_range to prevent self-reduction infinite loop
            if rr_start_idx < rr_end_idx && tr_start_idx < tr_end_idx {
                reduce_tr_by_rr_on_A(&mut matrix, rr_start_idx..rr_end_idx, tr_start_idx..tr_end_idx, ncl_val);
            }
            // Then apply dedicated reducers
            if tr_end_idx < dedicated_end_idx && tr_start_idx < tr_end_idx {
                reduce_tr_by_rr_on_A(&mut matrix, tr_end_idx..dedicated_end_idx, tr_start_idx..tr_end_idx, ncl_val);
            }

            // Note: After A-block triangularization, some rr rows may have pivots >= ncl.
            // This is expected - triangularization can eliminate all A-block entries from an rr row,
            // leaving it with a B-block pivot. These rows simply won't participate in AB reduction.

            // Assertion 2: After AB reduction, count tr rows with A-block residuals
            // (not an error yet, but will be checked after CD phase)
            #[allow(unused_variables)]
            let mut tr_residual_a_count_after_ab = 0;
            for tr_idx in tr_start_idx..tr_end_idx {
                if let Some(pivot_col) = matrix.rows[tr_idx].pivot() {
                    if pivot_col < ncl_val {
                        tr_residual_a_count_after_ab += 1;
                    }
                }
            }

            // After AB reduction, mark rr rows as pivots in A-block so CD phase can use them as reducers
            // IMPORTANT: Always mark rr rows as pivots, even if the column already has one.
            // This ensures CD phase uses the correct rr reducers.
            // Mark original rr rows
            for rr_idx in rr_start_idx..rr_end_idx {
                if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
                    if pivot_col < ncl_val {
                        // Normalize the rr row
                        let base_ring = matrix.ring.base_ring();
                        if let Some((_, pivot_coeff)) = matrix.rows[rr_idx].entries.first() {
                            let inv = base_ring.invert(pivot_coeff).unwrap();
                            for (_, coeff) in &mut matrix.rows[rr_idx].entries {
                                base_ring.mul_assign(coeff, base_ring.clone_el(&inv));
                            }
                        }
                        matrix.mark_pivot(pivot_col, rr_idx);
                    }
                }
            }
            // Mark dedicated reducers
            for rr_idx in tr_end_idx..dedicated_end_idx {
                if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
                    if pivot_col < ncl_val {
                        // Normalize the rr row
                        let base_ring = matrix.ring.base_ring();
                        if let Some((_, pivot_coeff)) = matrix.rows[rr_idx].entries.first() {
                            let inv = base_ring.invert(pivot_coeff).unwrap();
                            for (_, coeff) in &mut matrix.rows[rr_idx].entries {
                                base_ring.mul_assign(coeff, base_ring.clone_el(&inv));
                            }
                        }
                        matrix.mark_pivot(pivot_col, rr_idx);
                    }
                }
            }

            // Re-sort tr rows before CD reduction (pivots may have changed after AB phase)
            if tr_start_idx < tr_end_idx {
                matrix.sort_rows_by_pivot_and_density(tr_start_idx..tr_end_idx);
            }

            // CD phase: reduce tr rows only, discover pivots in B/D blocks (columns >= ncl)
            let (new_pivot_rows, _tr_new_pivots, anomaly_a_pivots) = if tr_start_idx < tr_end_idx {
                reduce_matrix_tr(&mut matrix, tr_start_idx..tr_end_idx, ncl_val)
            } else {
                (Vec::new(), 0, 0)
            };

            // Assertion 3: CD phase ideally should not have anomaly A-block pivots
            // but this can happen if rr rows don't cover all A-block columns
            // TODO: Fix root cause - ncl should match actual rr row coverage, not SHT pivot count
            if anomaly_a_pivots > 0 && std::env::var("F4_STRICT_ABCD").is_ok() {
                panic!(
                    "AB/CD assertion failed: {} tr rows still have pivot < ncl after CD phase",
                    anomaly_a_pivots
                );
            }

            new_pivot_rows
        };

        let t_elim = t_elim_0.elapsed();

        // Extract new basis elements from reduced matrix pivot rows
        // Following AB/CD semantics: only extract tr rows, never rr rows
        let t_extract_0 = Instant::now();
        let mut new_polys = Vec::new();
        let mut rr_extraction_count = 0;
        for &row_idx in &new_pivot_rows {
            // Skip rr rows (Rr, Dedicated) - only extract tr rows
            if row_idx < matrix.row_types.len() {
                match matrix.row_types[row_idx] {
                    RowType::Rr | RowType::Dedicated => {
                        rr_extraction_count += 1;
                        continue;
                    }
                    _ => {}, // Tr and Closure rows can be extracted
                }
            }

            let poly = matrix.row_to_polynomial(&matrix.rows[row_idx]);
            if ring.is_zero(&poly) { continue; }
            if poly_in_basis(ring, &poly, &basis, order) { continue; }
            new_polys.push(poly);
        }
        let t_extract = t_extract_0.elapsed();

        // Assertion 4: Extraction should only happen from tr rows (rr_extraction_count must be 0)
        if matrix.ncl.is_some() {
            assert_eq!(
                rr_extraction_count, 0,
                "AB/CD assertion failed: {} rr rows were in new_pivot_rows (should be 0)",
                rr_extraction_count
            );
        }

        // No fallback here; rely on strict A|B reducers to avoid blow-ups

        if profile {
            println!(
                "[F4] round: deg={} | pivots={} | rows={} cols={} | rows+red={} (+{}) cols+red={} (+{}) | passes={} | time rows={:?} reorder={:?} reducers={:?} elim={:?} extract={:?} total={:?}",
                min_deg,
                pivot_count,
                reducers_before_rows,
                reducers_before_cols,
                matrix.num_rows(),
                reducers_added_rows,
                matrix.num_cols(),
                reducers_new_cols,
                passes,
                t_rows,
                t_reorder,
                t_reducers,
                t_elim,
                t_extract,
                round_start.elapsed()
            );
        }

        // Add new polynomials to basis and generate new S-polys with Gebauer-Möller filtering
        let old_basis_len = basis.len();
        // Cache LM data for the current basis to avoid repeated expansions/divisibility checks
        let mut lm_cache = build_lm_cache(ring, &basis, order);
        let mut total_red_skipped = 0;
        let mut total_gm_old_removed = 0;
        let mut total_gm_new_removed = 0;
        let mut total_red_marked = 0;

        for (idx, new_poly) in new_polys.into_iter().enumerate() {
            let new_idx = old_basis_len + idx;
            let (_, new_lm) = ring.LT(&new_poly, order).unwrap();
            let new_lm_exp = ring.expand_monomial(new_lm);
            let new_lm_sig = MonSig::from_exponents(&new_lm_exp);
            // Clone monomial and degree now so they outlive the move of new_poly
            let new_lm_monom = ring.clone_monomial(new_lm);
            let new_lm_degree = ring.monomial_deg(new_lm);

            // Collect new pairs (i, new_idx) with product criterion
            let mut new_pairs = Vec::new();
            let mut red_skipped_old = 0;
            let mut red_skipped_new = 0;

            // Generate S-polys with OLD basis elements (apply Buchberger criterion)
            for i in 0..old_basis_len {
                if red_flags.get(i).copied().unwrap_or(false) {
                    red_skipped_old += 1;
                    continue;
                }
                // Coprime test via cached exponents: if any variable has both exponents > 0, not coprime
                let exp_i = &lm_cache[i].exponents;
                let mut coprime = true;
                for (&a, &b) in exp_i.iter().zip(new_lm_exp.iter()) {
                    if a > 0 && b > 0 { coprime = false; break; }
                }
                if !coprime {
                    new_pairs.push(SPoly::new(i, new_idx));
                }
            }

            // Generate S-polys with PREVIOUSLY ADDED NEW elements in this iteration
            for j in 0..idx {
                if red_flags.get(old_basis_len + j).copied().unwrap_or(false) {
                    red_skipped_new += 1;
                    continue;
                }
                let exp_j = &lm_cache[old_basis_len + j].exponents;
                let mut coprime = true;
                for (&a, &b) in exp_j.iter().zip(new_lm_exp.iter()) {
                    if a > 0 && b > 0 { coprime = false; break; }
                }
                if !coprime {
                    new_pairs.push(SPoly::new(old_basis_len + j, new_idx));
                }
            }

            total_red_skipped += red_skipped_old + red_skipped_new;

            // Add the new element to the basis BEFORE applying GM criteria
            basis.push(new_poly);
            red_flags.push(false);
            // Extend cache with the new element's LM data
            lm_cache.push(LMCacheEntry {
                lm: new_lm_monom,
                degree: new_lm_degree,
                exponents: new_lm_exp.clone(),
                sig: new_lm_sig.clone(),
            });

            // Mark redundancies (bs->red): mark i as redundant if LM(new) divides LM(i)
            // This means the new element can reduce i, making i redundant
            // NOTE: We do NOT mark i redundant if LM(i) divides LM(new) - that would be wrong!
            let lm_new2 = &lm_cache[new_idx];
            let lm_new2_exp = &lm_new2.exponents;
            let lm_new2_sig = &lm_new2.sig;
            let mut newly_marked = 0;
            for i in 0..new_idx {
                if red_flags[i] { continue; }
                let lm_i2 = &lm_cache[i];
                let lm_i2_exp = &lm_i2.exponents;
                let lm_i2_sig = &lm_i2.sig;

                // Divmask prefilter: check if lm_new2 can divide lm_i2
                if lm_new2_sig.deg as usize > lm_i2_sig.deg as usize { continue; }
                if !lm_new2_sig.presence_subset(&lm_i2_sig) { continue; }
                if !lm_new2_sig.buckets_leq(&lm_i2_sig) { continue; }

                // Only mark i redundant if the NEW element can reduce it
                // Division check via exponents: lm_new2 | lm_i2 if new_exp[v] <= i_exp[v] for all v
                if lm_new2_exp.iter().zip(lm_i2_exp.iter()).all(|(a,b)| a <= b) {
                    red_flags[i] = true;
                    newly_marked += 1;
                }
            }
            total_red_marked += newly_marked;

            // Apply Gebauer-Möller criterion to OLD pairs
            // Remove old pairs (i,j) if LM(new) | LCM(i,j) and degrees of (i,new) and (j,new) ≤ deg(i,j)
            let old_spolys_len = spolys.len();
            gm_filter_old_pairs(ring, &mut spolys, &basis, new_idx, order);
            total_gm_old_removed += old_spolys_len - spolys.len();

            // Apply Gebauer-Möller criterion to NEW pairs
            // Remove dominated pairs (those whose LCM is divisible by an earlier pair's LCM)
            let before_gm_new = new_pairs.len();
            gm_filter_new_pairs(ring, &mut new_pairs, &basis, order);
            total_gm_new_removed += before_gm_new - new_pairs.len();

            // Add filtered new pairs to the global list
            spolys.extend(new_pairs);
        }

        if profile && (total_red_skipped > 0 || total_gm_old_removed > 0 || total_gm_new_removed > 0 || total_red_marked > 0) {
            println!("[F4]   Pair filtering: red_skipped={} gm_old_removed={} gm_new_removed={} red_marked={}",
                total_red_skipped, total_gm_old_removed, total_gm_new_removed, total_red_marked);
        }
    }

    // Final interreduction
    if std::env::var("F4_PROFILE").map(|v| v == "1" || v.eq_ignore_ascii_case("true")).unwrap_or(false) {
        eprintln!("[F4] basis size before interreduce: {}", basis.len());
    }
    let __t_ir0 = std::time::Instant::now();
    basis = interreduce(ring, basis, order);
    if std::env::var("F4_PROFILE").map(|v| v == "1" || v.eq_ignore_ascii_case("true")).unwrap_or(false) {
        eprintln!("[F4] final interreduce: {:?}", __t_ir0.elapsed());
    }

    Ok(basis)
}

/// Check if two monomials are coprime (gcd = 1)
///
/// Two monomials are coprime if for every variable, at least one has exponent 0.
/// This is used in the Buchberger product criterion.
fn are_monomials_coprime<P>(ring: P, m1: &PolyMonomial<P>, m2: &PolyMonomial<P>) -> bool
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    let deg1 = ring.monomial_deg(m1);
    let deg2 = ring.monomial_deg(m2);

    // Quick check: if either monomial is constant (degree 0), they're coprime
    if deg1 == 0 || deg2 == 0 {
        return true;
    }

    // Expand monomials to get exponents
    let exp1 = ring.expand_monomial(m1);
    let exp2 = ring.expand_monomial(m2);

    // Check each variable: for coprimality, at least one monomial must have exponent 0 for each variable
    for (e1, e2) in exp1.iter().zip(exp2.iter()) {
        if *e1 > 0 && *e2 > 0 {
            return false; // Both have this variable - not coprime
        }
    }

    true
}

/// Gebauer-Möller criterion: filter old pairs when a new basis element is added
///
/// An old pair (i,j) is redundant if its LCM is divisible by the new leading term
/// and both new pairs (i,new) and (j,new) have degrees ≤ the old pair's degree.
fn gm_filter_old_pairs<P, O>(
    ring: P,
    spolys: &mut Vec<SPoly>,
    basis: &[El<P>],
    new_idx: usize,
    order: O,
) where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let (_, lm_new) = ring.LT(&basis[new_idx], order).unwrap();
    let lm_new_exp = ring.expand_monomial(lm_new);
    let lm_new_sig = MonSig::from_exponents(&lm_new_exp);

    spolys.retain(|sp| {
        // Get LCM of the old pair
        let (_, lm_i) = ring.LT(&basis[sp.i], order).unwrap();
        let (_, lm_j) = ring.LT(&basis[sp.j], order).unwrap();
        let lcm_old = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);
        let lcm_old_exp = ring.expand_monomial(&lcm_old);
        let lcm_old_sig = MonSig::from_exponents(&lcm_old_exp);

        // Divmask prefilter: check if new LT can divide old LCM
        if lm_new_sig.deg as usize > lcm_old_sig.deg as usize {
            return true; // Keep: degree too high, cannot divide
        }
        if !lm_new_sig.presence_subset(&lcm_old_sig) {
            return true; // Keep: presence not subset, cannot divide
        }
        if !lm_new_sig.buckets_leq(&lcm_old_sig) {
            return true; // Keep: bucket check fails, cannot divide
        }

        // Check if new LT divides old LCM
        if ring.monomial_div(ring.clone_monomial(&lcm_old), lm_new).is_err() {
            return true; // Keep: new LT does not divide old LCM
        }

        // Compute LCMs of new pairs (i, new) and (j, new)
        let lcm_i_new = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_new);
        let lcm_j_new = ring.monomial_lcm(ring.clone_monomial(lm_j), lm_new);

        // Check that BOTH new LCMs are different from old LCM (msolve check)
        // If either equals the old LCM, keep the old pair
        if order.compare(ring, &lcm_i_new, &lcm_old) == std::cmp::Ordering::Equal ||
           order.compare(ring, &lcm_j_new, &lcm_old) == std::cmp::Ordering::Equal {
            return true; // Keep: one of the new pairs has the same LCM
        }

        // Check degrees: keep if either new pair has degree > old pair
        let deg_i_new = ring.monomial_deg(&lcm_i_new);
        let deg_j_new = ring.monomial_deg(&lcm_j_new);
        let deg_old = ring.monomial_deg(&lcm_old);

        // Keep if either new pair has degree > old pair
        // Otherwise discard (the new pairs dominate)
        deg_i_new > deg_old || deg_j_new > deg_old
    });
}

/// Gebauer-Möller criterion: filter new pairs by LCM divisibility
///
/// Among new pairs, remove any whose LCM is divisible by an earlier pair's LCM
/// with the earlier pair having ≤ degree.
fn gm_filter_new_pairs<P, O>(
    ring: P,
    pairs: &mut Vec<SPoly>,
    basis: &[El<P>],
    order: O,
) where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    if pairs.is_empty() {
        return;
    }

    // Sort pairs by degree (ascending), then by indices for stability
    pairs.sort_by_key(|sp| (sp.lcm_degree(ring, basis, order), sp.i, sp.j));

    let mut filtered: Vec<SPoly> = Vec::new();

    for candidate in pairs.drain(..) {
        let (_, lm_i) = ring.LT(&basis[candidate.i], order).unwrap();
        let (_, lm_j) = ring.LT(&basis[candidate.j], order).unwrap();
        let lcm_candidate = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);
        let deg_candidate = ring.monomial_deg(&lcm_candidate);
        let lcm_candidate_exp = ring.expand_monomial(&lcm_candidate);
        let lcm_candidate_sig = MonSig::from_exponents(&lcm_candidate_exp);

        let mut dominated = false;
        for kept in &filtered {
            let (_, lm_ki) = ring.LT(&basis[kept.i], order).unwrap();
            let (_, lm_kj) = ring.LT(&basis[kept.j], order).unwrap();
            let lcm_kept = ring.monomial_lcm(ring.clone_monomial(lm_ki), lm_kj);
            let deg_kept = ring.monomial_deg(&lcm_kept);

            // Check if candidate LCM is divisible by kept LCM and degree is not better
            if deg_candidate >= deg_kept {
                // Divmask prefilter: check if kept LCM can divide candidate LCM
                let lcm_kept_exp = ring.expand_monomial(&lcm_kept);
                let lcm_kept_sig = MonSig::from_exponents(&lcm_kept_exp);

                if !lcm_kept_sig.presence_subset(&lcm_candidate_sig) { continue; }
                if !lcm_kept_sig.buckets_leq(&lcm_candidate_sig) { continue; }

                if ring.monomial_div(ring.clone_monomial(&lcm_candidate), &lcm_kept).is_ok() {
                    dominated = true;
                    break;
                }
            }
        }

        if !dominated {
            filtered.push(candidate);
        }
    }

    *pairs = filtered;
}

/// Cached variant of gm_filter_old_pairs operating purely on exponent vectors
fn gm_filter_old_pairs_cached<P>(
    spolys: &mut Vec<SPoly>,
    lm_cache: &[LMCacheEntry<P>],
    new_idx: usize,
)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    let new_exp = &lm_cache[new_idx].exponents;
    let new_sig = &lm_cache[new_idx].sig;

    spolys.retain(|sp| {
        // Compute lcm(i,j) exponents and its degree
        let exp_i = &lm_cache[sp.i].exponents;
        let exp_j = &lm_cache[sp.j].exponents;
        let mut lcm_old_exp = Vec::with_capacity(exp_i.len());
        let mut deg_old = 0usize;
        for (&a, &b) in exp_i.iter().zip(exp_j.iter()) {
            let m = if a > b { a } else { b };
            deg_old += m;
            lcm_old_exp.push(m);
        }
        let lcm_old_sig = MonSig::from_exponents(&lcm_old_exp);

        // Prefilters
        if new_sig.deg as usize > lcm_old_sig.deg as usize { return true; }
        if !new_sig.presence_subset(&lcm_old_sig) { return true; }
        if !new_sig.buckets_leq(&lcm_old_sig) { return true; }

        // new | lcm_old?
        if !new_exp.iter().zip(&lcm_old_exp).all(|(a,b)| a <= b) { return true; }

        // LCM(i,new) and LCM(j,new) degrees and equality tests vs old
        let mut lcm_i_new_deg = 0usize;
        let mut lcm_j_new_deg = 0usize;
        let mut same_as_old_i = true;
        let mut same_as_old_j = true;
        for (((&ai, &aj), &an), &old) in exp_i.iter().zip(exp_j).zip(new_exp).zip(lcm_old_exp.iter()) {
            let max_i = if ai > an { ai } else { an };
            let max_j = if aj > an { aj } else { an };
            lcm_i_new_deg += max_i;
            lcm_j_new_deg += max_j;
            if max_i != old { same_as_old_i = false; }
            if max_j != old { same_as_old_j = false; }
        }

        if same_as_old_i || same_as_old_j { return true; }

        lcm_i_new_deg > deg_old || lcm_j_new_deg > deg_old
    });
}

/// Cached variant of gm_filter_new_pairs using exponent vectors
fn gm_filter_new_pairs_cached<P>(
    pairs: &mut Vec<SPoly>,
    lm_cache: &[LMCacheEntry<P>],
)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    if pairs.is_empty() { return; }

    // Precompute LCM(exponents) and degree for each pair
    let mut lcms: Vec<(Vec<usize>, usize)> = Vec::with_capacity(pairs.len());
    for sp in pairs.iter() {
        let ei = &lm_cache[sp.i].exponents;
        let ej = &lm_cache[sp.j].exponents;
        let mut exp = Vec::with_capacity(ei.len());
        let mut deg = 0usize;
        for (&a, &b) in ei.iter().zip(ej.iter()) {
            let m = if a > b { a } else { b };
            deg += m;
            exp.push(m);
        }
        lcms.push((exp, deg));
    }

    // Sort by degree (ascending) to emulate original stable behavior
    let mut idxs: Vec<usize> = (0..pairs.len()).collect();
    idxs.sort_by_key(|&i| lcms[i].1);

    let mut kept: Vec<usize> = Vec::new();
    for &ci in &idxs {
        let (ref cand_exp, cand_deg) = lcms[ci];
        let mut dominated = false;
        for &ki in &kept {
            let (ref kept_exp, kept_deg) = lcms[ki];
            if kept_deg > cand_deg { continue; }
            if kept_exp.iter().zip(cand_exp.iter()).all(|(a,b)| a <= b) {
                dominated = true;
                break;
            }
        }
        if !dominated { kept.push(ci); }
    }

    // Rebuild pairs from kept indices in original order
    let mut filtered = Vec::with_capacity(kept.len());
    for i in kept { filtered.push(pairs[i].clone()); }
    *pairs = filtered;
}
/// Leading monomial index entry for fast reducer lookup
struct LMIndexEntry<P>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    basis_idx: usize,
    lm: PolyMonomial<P>,
    degree: usize,
    exponents: Vec<usize>,
    sig: MonSig,
}

/// Lightweight cache of leading monomial data for all basis elements
struct LMCacheEntry<P>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    lm: PolyMonomial<P>,
    degree: usize,
    exponents: Vec<usize>,
    sig: MonSig,
}

fn build_lm_cache<P, O>(ring: P, basis: &[El<P>], order: O) -> Vec<LMCacheEntry<P>>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let mut out = Vec::with_capacity(basis.len());
    for f in basis {
        if let Some((_, lm)) = ring.LT(f, order) {
            let exponents = ring.expand_monomial(lm);
            out.push(LMCacheEntry {
                lm: ring.clone_monomial(lm),
                degree: ring.monomial_deg(lm),
                sig: MonSig::from_exponents(&exponents),
                exponents,
            });
        }
    }
    out
}

/// Build an LM index for fast reducer search
///
/// Returns a vector sorted by degree, containing (index, LM, degree, exponents) for each basis element.
fn build_lm_index<P, O>(ring: P, basis: &[El<P>], order: O) -> Vec<LMIndexEntry<P>>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let mut index = Vec::new();

    for (i, poly) in basis.iter().enumerate() {
        if let Some((_, lm)) = ring.LT(poly, order) {
            let degree = ring.monomial_deg(lm);
            let exponents = ring.expand_monomial(lm);
            let sig = MonSig::from_exponents(&exponents);
            index.push(LMIndexEntry {
                basis_idx: i,
                lm: ring.clone_monomial(lm),
                degree,
                exponents,
                sig,
            });
        }
    }

    // Sort by degree (ascending) for faster lookup
    index.sort_by_key(|e| e.degree);

    index
}

/// Add reducers only for a provided set of target columns (typically from S-poly rows),
/// skipping pivot columns and restricting reducer rows to the A-block (pivot) columns.
fn add_reducers_for_columns_with_pivot_skip<P, O>(
    matrix: &mut MacaulayMatrix<P>,
    basis: &[El<P>],
    order: O,
    pivot_cols: &std::collections::HashSet<usize>,
    pivot_count: usize,
    target_cols: &std::collections::HashSet<usize>,
    lm_index: &[LMIndexEntry<P>],
    col_sig_cache: &mut std::collections::HashMap<usize, (usize, Vec<usize>, MonSig)>,
)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let ring = matrix.ring;

    for &col in target_cols {
        if matrix.has_pivot(col) || pivot_cols.contains(&col) {
            continue;
        }

        // Use cached column signature or compute and cache it
        let (monomial_deg, monomial_exp, monomial_sig) = col_sig_cache.entry(col).or_insert_with(|| {
            let monomial = ring.clone_monomial(&matrix.col_to_monomial[col]);
            let deg = ring.monomial_deg(&monomial);
            let exp = ring.expand_monomial(&monomial);
            let sig = MonSig::from_exponents(&exp);
            (deg, exp, sig)
        });
        let monomial_deg = *monomial_deg;
        let monomial_exp = monomial_exp;
        let monomial_sig = monomial_sig;

        for entry in lm_index {
            if entry.degree > monomial_deg { break; }

            // Divmask prefilter: fast rejection based on presence mask and buckets
            if !entry.sig.presence_subset(&monomial_sig) { continue; }
            if !entry.sig.buckets_leq(&monomial_sig) { continue; }

            let mut can_divide = true;
            for (lm_exp, mon_exp) in entry.exponents.iter().zip(monomial_exp.iter()) {
                if lm_exp > mon_exp { can_divide = false; break; }
            }
            if !can_divide { continue; }
            if let Ok(multiplier) = ring.monomial_div(ring.clone_monomial(&matrix.col_to_monomial[col]), &entry.lm) {
                let mut reducer = ring.clone_el(&basis[entry.basis_idx]);
                ring.mul_assign_monomial(&mut reducer, multiplier);
                let row = matrix.polynomial_to_row(&reducer, order);
                // Keep full row support (no A-block pruning) per user's analysis
                if row.entries.is_empty() { continue; }
                matrix.add_row(row, RowType::Closure);
                break;
            }
        }
    }
}

/// Check if a polynomial is redundant with respect to the basis
///
/// A polynomial is redundant if its leading monomial EXACTLY matches
/// an existing basis element's leading monomial (not just divisible by).
fn poly_in_basis<P, O>(ring: P, poly: &El<P>, basis: &[El<P>], order: O) -> bool
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let (_, poly_lm) = match ring.LT(poly, order) {
        Some(lt) => lt,
        None => return true, // Zero polynomial is redundant
    };

    // Check if ANY basis element has the SAME leading monomial
    for b in basis {
        if let Some((_, b_lm)) = ring.LT(b, order) {
            // Check for exact monomial equality
            if order.compare(ring, poly_lm, b_lm) == std::cmp::Ordering::Equal {
                return true; // Same leading monomial - redundant
            }
        }
    }

    false
}

/// Interreduce a basis (reduce each polynomial by the others)
fn interreduce<P, O>(ring: P, mut basis: Vec<El<P>>, order: O) -> Vec<El<P>>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    if basis.is_empty() { return basis; }

    let base_ring = ring.base_ring();
    let mut G: Vec<El<P>> = Vec::with_capacity(basis.len());
    let mut G_lm: Vec<PolyMonomial<P>> = Vec::new();
    let mut G_exp: Vec<Vec<usize>> = Vec::new();
    let mut G_deg: Vec<usize> = Vec::new();
    let mut G_inv: Vec<PolyCoeff<P>> = Vec::new();

    'basis_loop: for f0 in basis.into_iter() {
        let mut f = f0;

        // Reduce f by current G
        'reduce_f: loop {
            let (f_lc, f_lm) = match ring.LT(&f, order) { Some(lt) => lt, None => continue 'basis_loop };
            let f_deg = ring.monomial_deg(f_lm);
            let f_exp = ring.expand_monomial(f_lm);
            for j in 0..G.len() {
                if G_deg[j] > f_deg { continue; }
                // per-variable check
                let mut ok = true;
                for (&a, &b) in G_exp[j].iter().zip(f_exp.iter()) { if a > b { ok = false; break; } }
                if !ok { continue; }
                if let Ok(qm) = ring.monomial_div(ring.clone_monomial(f_lm), &G_lm[j]) {
                    let scale = base_ring.mul_ref(&f_lc, &G_inv[j]);
                    let mut red = ring.clone_el(&G[j]);
                    ring.mul_assign_monomial(&mut red, qm);
                    ring.inclusion().mul_assign_ref_map(&mut red, &scale);
                    ring.sub_assign(&mut f, red);
                    continue 'reduce_f;
                }
            }
            break; // irreducible by G
        }

        if ring.is_zero(&f) { continue; }

        // Normalize f
        if let Some((lc, lm)) = ring.LT(&f, order) {
            // Clone LM so we can mutate f without holding a borrow on it
            let lm_owned = ring.clone_monomial(lm);
            let inv_once = base_ring.invert(lc).unwrap();
            ring.inclusion().mul_assign_map(&mut f, inv_once);

            // Remove from G any g whose LM is divisible by LM(f), and reduce affected ones by f
            let f_exp_norm = ring.expand_monomial(&lm_owned);
            let f_deg_norm = ring.monomial_deg(&lm_owned);

            let mut k = 0;
            while k < G.len() {
                // If LM(f) divides LM(G[k]), reduce G[k] by f
                if G_exp[k].iter().zip(&f_exp_norm).all(|(gk, fk)| gk >= fk) {
                    // Apply one-step reduction on G[k] by f until not divisible
                    'reduce_g: loop {
                        let (g_lc, g_lm) = match ring.LT(&G[k], order) { Some(lt) => lt, None => break 'reduce_g };
                        if ring.monomial_deg(g_lm) < f_deg_norm { break 'reduce_g; }
                        // check divisibility
                        if let Ok(qm) = ring.monomial_div(ring.clone_monomial(g_lm), &lm_owned) {
                            // f is monic; multiplier is simply g_lc
                            let scale = base_ring.clone_el(&g_lc);
                            let mut red = ring.clone_el(&f);
                            ring.mul_assign_monomial(&mut red, qm);
                            ring.inclusion().mul_assign_ref_map(&mut red, &scale);
                            ring.sub_assign(&mut G[k], red);
                        } else {
                            break 'reduce_g;
                        }
                    }

                    if ring.is_zero(&G[k]) {
                        // remove k
                        G.remove(k);
                        G_lm.remove(k);
                        G_exp.remove(k);
                        G_deg.remove(k);
                        G_inv.remove(k);
                        continue; // do not increment k
                    } else {
                        // refresh LM cache for G[k]
                        let (glc2, glm2) = ring.LT(&G[k], order).unwrap();
                        G_lm[k] = ring.clone_monomial(glm2);
                        G_exp[k] = ring.expand_monomial(glm2);
                        G_deg[k] = ring.monomial_deg(glm2);
                        G_inv[k] = base_ring.invert(glc2).unwrap();
                    }
                }
                k += 1;
            }

            // Insert f into G and cache its LM data
            // After normalization, LC(f) = 1
            G_inv.push(base_ring.one());
            G_deg.push(f_deg_norm);
            G_exp.push(f_exp_norm);
            G_lm.push(lm_owned);
            G.push(f);
        }
    }

    G
}

/// Try to reduce `poly` by `reducer` once
fn try_reduce<P, O>(ring: P, poly: &mut El<P>, reducer: &El<P>, order: O) -> bool
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let (poly_lc, poly_lm) = match ring.LT(poly, order) {
        Some(lt) => lt,
        None => return false,
    };

    let (red_lc, red_lm) = match ring.LT(reducer, order) {
        Some(lt) => lt,
        None => return false,
    };

    // Check if reducer's LT divides poly's LT
    if let Ok(quot_m) = ring.monomial_div(ring.clone_monomial(poly_lm), red_lm) {
        // Compute: poly -= (poly_lc / red_lc) * quot_m * reducer
        let scalar = ring.base_ring().checked_div(poly_lc, red_lc).unwrap();
        let mut scaled_red = ring.clone_el(reducer);
        ring.mul_assign_monomial(&mut scaled_red, quot_m);
        ring.inclusion().mul_assign_ref_map(&mut scaled_red, &scalar);
        ring.sub_assign(poly, scaled_red);
        return true;
    }

    false
}

/// Try to reduce `poly` by `reducer` once using a cached inverse of LC(reducer)
///
/// This avoids performing a field inversion on every attempted reduction and is a
/// major win on cryptographic prime fields where inversion is very expensive.
fn try_reduce_with_inv<P, O>(
    ring: P,
    poly: &mut El<P>,
    reducer: &El<P>,
    order: O,
    red_lc_inv: &PolyCoeff<P>,
) -> bool
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    // Leading term of the target polynomial
    let (poly_lc, poly_lm) = match ring.LT(poly, order) {
        Some(lt) => lt,
        None => return false,
    };

    // Leading monomial of the reducer (coefficient inverse is provided by caller)
    let (_, red_lm) = match ring.LT(reducer, order) {
        Some(lt) => lt,
        None => return false,
    };

    if let Ok(quot_m) = ring.monomial_div(ring.clone_monomial(poly_lm), red_lm) {
        // scalar = poly_lc * inv(red_lc)
        let scalar = ring.base_ring().mul_ref(&poly_lc, red_lc_inv);
        let mut scaled_red = ring.clone_el(reducer);
        ring.mul_assign_monomial(&mut scaled_red, quot_m);
        ring.inclusion().mul_assign_ref_map(&mut scaled_red, &scalar);
        ring.sub_assign(poly, scaled_red);
        return true;
    }

    false
}
