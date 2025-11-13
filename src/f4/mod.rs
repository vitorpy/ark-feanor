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

use crate::f4::gaussian::reduce_matrix;
use crate::f4::matrix::MacaulayMatrix;
use feanor_math::algorithms::buchberger::GBAborted;
use feanor_math::divisibility::DivisibilityRingStore;
use feanor_math::field::Field;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::ring::*;
use feanor_math::rings::multivariate::*;
use std::time::Duration;

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
}

impl F4Config {
    /// Create a new default configuration
    pub fn new() -> Self {
        F4Config {
            max_degree: None,
            max_matrix_rows: None,
            s_pair_budget: None,
            time_budget: None,
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
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    f4_configured(ring, input_basis, order, F4Config::default()).unwrap()
}

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
    O: MonomialOrder + Copy,
{
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

        // Select S-polys of current degree
        let (current_spolys, remaining_spolys): (Vec<_>, Vec<_>) = spolys
            .into_iter()
            .partition(|sp| sp.lcm_degree(ring, &basis, order) == min_deg);

        spolys = remaining_spolys;

        // Check S-pair budget
        s_pairs_processed += current_spolys.len();
        if let Some(max_pairs) = config.s_pair_budget {
            if s_pairs_processed > max_pairs {
                return Err(GBAborted::SPairBudget { max_s_pairs: max_pairs });
            }
        }

        // Build Macaulay matrix
        let mut matrix = MacaulayMatrix::new(ring);

        // Add S-polynomials to matrix
        for sp in &current_spolys {
            let spoly = sp.compute(ring, &basis, order);
            if !ring.is_zero(&spoly) {
                let row = matrix.polynomial_to_row(&spoly, order);
                matrix.add_row(row);
            }
        }

        // Compute LCM monomials for S-polys (intended pivot columns)
        let pivot_monomials: Vec<_> = current_spolys
            .iter()
            .map(|sp| {
                let (_, lm_i) = ring.LT(&basis[sp.i], order).unwrap();
                let (_, lm_j) = ring.LT(&basis[sp.j], order).unwrap();
                ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j)
            })
            .collect();

        // Reorder columns: put pivot columns (LCMs) first for A|B split
        // This reduces reducer pressure and improves elimination efficiency
        matrix.reorder_columns_pivot_first(&pivot_monomials);

        // Identify pivot column indices (after reordering, these are now the leftmost columns)
        let pivot_col_indices: std::collections::HashSet<usize> = pivot_monomials
            .iter()
            .filter_map(|m| {
                let expanded = ring.expand_monomial(m);
                matrix.monomial_to_col.get(&expanded).copied()
            })
            .collect();

        // Add reducers (basis multiples) to matrix, skipping pivot columns
        add_reducers_to_matrix_with_pivot_skip(&mut matrix, &basis, order, &pivot_col_indices);

        // Check matrix size limit
        if let Some(max_rows) = config.max_matrix_rows {
            if matrix.num_rows() > max_rows {
                return Err(GBAborted::MemBudget {
                    max_bytes: max_rows * std::mem::size_of::<usize>(),
                });
            }
        }

        // Reduce the matrix via Gaussian elimination
        let new_pivot_rows = reduce_matrix(&mut matrix);

        // Extract new basis elements from reduced matrix
        let mut new_polys = Vec::new();
        for &row_idx in &new_pivot_rows {
            let poly = matrix.row_to_polynomial(&matrix.rows[row_idx]);
            if ring.is_zero(&poly) {
                continue;
            }
            if poly_in_basis(ring, &poly, &basis, order) {
                continue;
            }
            new_polys.push(poly);
        }

        // Add new polynomials to basis and generate new S-polys with Gebauer-Möller filtering
        let old_basis_len = basis.len();
        for (idx, new_poly) in new_polys.into_iter().enumerate() {
            let new_idx = old_basis_len + idx;
            let (_, new_lm) = ring.LT(&new_poly, order).unwrap();

            // Collect new pairs (i, new_idx) with product criterion
            let mut new_pairs = Vec::new();

            // Generate S-polys with OLD basis elements (apply Buchberger criterion)
            for i in 0..old_basis_len {
                let (_, lm_i) = ring.LT(&basis[i], order).unwrap();
                if !are_monomials_coprime(ring, lm_i, new_lm) {
                    new_pairs.push(SPoly::new(i, new_idx));
                }
            }

            // Generate S-polys with PREVIOUSLY ADDED NEW elements in this iteration
            for j in 0..idx {
                let (_, lm_j) = ring.LT(&basis[old_basis_len + j], order).unwrap();
                if !are_monomials_coprime(ring, lm_j, new_lm) {
                    new_pairs.push(SPoly::new(old_basis_len + j, new_idx));
                }
            }

            // Add the new element to the basis BEFORE applying GM criteria
            basis.push(new_poly);

            // Apply Gebauer-Möller criterion to OLD pairs
            // Remove old pairs (i,j) if LM(new) | LCM(i,j) and degrees of (i,new) and (j,new) ≤ deg(i,j)
            gm_filter_old_pairs(ring, &mut spolys, &basis, new_idx, order);

            // Apply Gebauer-Möller criterion to NEW pairs
            // Remove dominated pairs (those whose LCM is divisible by an earlier pair's LCM)
            gm_filter_new_pairs(ring, &mut new_pairs, &basis, order);

            // Add filtered new pairs to the global list
            spolys.extend(new_pairs);
        }
    }

    // Final interreduction
    basis = interreduce(ring, basis, order);

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

    spolys.retain(|sp| {
        // Get LCM of the old pair
        let (_, lm_i) = ring.LT(&basis[sp.i], order).unwrap();
        let (_, lm_j) = ring.LT(&basis[sp.j], order).unwrap();
        let lcm_old = ring.monomial_lcm(ring.clone_monomial(lm_i), lm_j);

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

        let mut dominated = false;
        for kept in &filtered {
            let (_, lm_ki) = ring.LT(&basis[kept.i], order).unwrap();
            let (_, lm_kj) = ring.LT(&basis[kept.j], order).unwrap();
            let lcm_kept = ring.monomial_lcm(ring.clone_monomial(lm_ki), lm_kj);
            let deg_kept = ring.monomial_deg(&lcm_kept);

            // Check if candidate LCM is divisible by kept LCM and degree is not better
            if deg_candidate >= deg_kept {
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
            index.push(LMIndexEntry {
                basis_idx: i,
                lm: ring.clone_monomial(lm),
                degree,
                exponents,
            });
        }
    }

    // Sort by degree (ascending) for faster lookup
    index.sort_by_key(|e| e.degree);

    index
}

/// Add necessary basis multiples to the matrix as reducers
///
/// Following msolve's approach: for each monomial M in the matrix, find the FIRST
/// basis element whose leading term divides M, then add (M/LT(basis)) * basis as a reducer.
///
/// Uses an LM index for O(k) lookup per column where k = # basis elements with deg ≤ deg(M).
fn add_reducers_to_matrix<P, O>(matrix: &mut MacaulayMatrix<P>, basis: &[El<P>], order: O)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    use std::collections::HashSet;
    add_reducers_to_matrix_with_pivot_skip(matrix, basis, order, &HashSet::new())
}

/// Add necessary basis multiples to the matrix as reducers, skipping pivot columns
///
/// Following msolve's approach: for each monomial M in the matrix (except pivot columns),
/// find the FIRST basis element whose leading term divides M, then add (M/LT(basis)) * basis.
///
/// Pivot columns are S-poly LCMs that will become pivots; no reducers needed for them.
fn add_reducers_to_matrix_with_pivot_skip<P, O>(
    matrix: &mut MacaulayMatrix<P>,
    basis: &[El<P>],
    order: O,
    pivot_cols: &std::collections::HashSet<usize>,
)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    use std::collections::HashSet;

    let ring = matrix.ring;

    // Build LM index for fast degree-filtered lookup
    let lm_index = build_lm_index(ring, basis, order);

    // Collect all column indices (monomials) currently in the matrix
    let mut monomial_cols: HashSet<usize> = HashSet::new();
    for row in &matrix.rows {
        for &(col, _) in &row.entries {
            monomial_cols.insert(col);
        }
    }

    // For each monomial in the matrix, try to find a reducer
    for &col in &monomial_cols {
        // Skip if this column already has a pivot (already has a reducer)
        if matrix.has_pivot(col) {
            continue;
        }

        // Skip pivot columns (S-poly LCMs) - they will become pivots
        if pivot_cols.contains(&col) {
            continue;
        }

        let monomial = &matrix.col_to_monomial[col];
        let monomial_deg = ring.monomial_deg(monomial);
        let monomial_exp = ring.expand_monomial(monomial);

        // Find the first basis element whose LT divides this monomial
        // Only scan basis elements with deg(LM) ≤ deg(monomial) using LM index
        for entry in &lm_index {
            // Degree prefilter: skip if basis LM has higher degree
            if entry.degree > monomial_deg {
                break; // Index is sorted, so we can stop here
            }

            // Variable-wise upper bound check: if any exponent in LM > monomial, skip
            let mut can_divide = true;
            for (lm_exp, mon_exp) in entry.exponents.iter().zip(monomial_exp.iter()) {
                if lm_exp > mon_exp {
                    can_divide = false;
                    break;
                }
            }

            if !can_divide {
                continue;
            }

            // Full divisibility check
            if let Ok(multiplier) = ring.monomial_div(ring.clone_monomial(monomial), &entry.lm) {
                // Compute (monomial / LT(basis)) * basis
                let mut reducer = ring.clone_el(&basis[entry.basis_idx]);
                ring.mul_assign_monomial(&mut reducer, multiplier);

                let row = matrix.polynomial_to_row(&reducer, order);
                matrix.add_row(row);

                // Only add ONE reducer per monomial
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
    let mut changed = true;
    while changed {
        changed = false;
        let mut i = 0;
        while i < basis.len() {
            // Try to reduce basis[i] by all others
            let mut poly = ring.clone_el(&basis[i]);
            let mut reduced = false;

            for j in 0..basis.len() {
                if i == j {
                    continue;
                }
                if try_reduce(ring, &mut poly, &basis[j], order) {
                    reduced = true;
                }
            }

            if ring.is_zero(&poly) {
                basis.remove(i);
                changed = true;
            } else if reduced {
                basis[i] = poly;
                changed = true;
                i += 1;
            } else {
                i += 1;
            }
        }
    }

    // Final normalization
    for f in &mut basis {
        let (lc, _) = ring.LT(f, order).unwrap();
        let lc_inv = ring.base_ring().invert(lc).unwrap();
        ring.inclusion().mul_assign_map(f, lc_inv);
    }

    basis
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

#[cfg(test)]
mod tests {
    use super::*;
    use feanor_math::homomorphism::Homomorphism;
    use feanor_math::rings::multivariate::multivariate_impl::{MultivariatePolyRingImpl, DegreeCfg};
    use feanor_math::rings::multivariate::DegRevLex;
    use feanor_math::rings::zn::zn_static;
    use std::alloc::Global;

    #[test]
    fn test_f4_simple() {
        let base = zn_static::F17;
        // Use default ring (degree 64) but with lower max_degree config
        let ring = MultivariatePolyRingImpl::new(base, 2);

        // Simpler system to avoid degree overflow: x + y - 1, x*y - 2
        let f1 = ring.from_terms([
            (base.int_hom().map(1), ring.create_monomial([1, 0])),
            (base.int_hom().map(1), ring.create_monomial([0, 1])),
            (base.int_hom().map(16), ring.create_monomial([0, 0])),
        ].into_iter());

        let f2 = ring.from_terms([
            (base.int_hom().map(1), ring.create_monomial([1, 1])),
            (base.int_hom().map(15), ring.create_monomial([0, 0])),
        ].into_iter());

        // Add time budget to prevent infinite loops
        let config = F4Config::new().with_max_degree(10).with_time_budget(std::time::Duration::from_secs(5));
        let gb = f4_configured(&ring, vec![f1, f2], DegRevLex, config).unwrap();

        // Should compute a Gröbner basis
        assert!(!gb.is_empty());
        assert!(gb.len() >= 2);
    }

    #[test]
    fn test_f4_with_degree_limit() {
        let base = zn_static::F17;
        let ring = MultivariatePolyRingImpl::new(base, 2);

        let f1 = ring.from_terms([
            (base.int_hom().map(1), ring.create_monomial([2, 0])),
            (base.int_hom().map(1), ring.create_monomial([0, 0])),
        ].into_iter());

        let f2 = ring.from_terms([
            (base.int_hom().map(1), ring.create_monomial([1, 1])),
            (base.int_hom().map(1), ring.create_monomial([0, 0])),
        ].into_iter());

        let config = F4Config::new().with_max_degree(5);

        let result = f4_configured(&ring, vec![f1, f2], DegRevLex, config);

        // Should either succeed or abort with degree exceeded
        match result {
            Ok(gb) => assert!(!gb.is_empty()),
            Err(GBAborted::DegreeExceeded { .. }) => {}, // Also acceptable
            Err(e) => panic!("Unexpected error: {:?}", e),
        }
    }

    #[test]
    fn test_gebauer_moeller_product_criterion() {
        let base = zn_static::F17;
        let ring = MultivariatePolyRingImpl::new(base, 2);

        // Test product criterion: x and y are coprime
        let x_mono = ring.create_monomial([1, 0]);
        let y_mono = ring.create_monomial([0, 1]);

        assert!(are_monomials_coprime(&ring, &x_mono, &y_mono),
                "x and y should be coprime");

        // x^2 and y are coprime
        let x2_mono = ring.create_monomial([2, 0]);
        assert!(are_monomials_coprime(&ring, &x2_mono, &y_mono),
                "x^2 and y should be coprime");

        // xy and y are NOT coprime (share y)
        let xy_mono = ring.create_monomial([1, 1]);
        assert!(!are_monomials_coprime(&ring, &xy_mono, &y_mono),
                "xy and y should not be coprime");

        // xy and x^2 are NOT coprime (share x)
        assert!(!are_monomials_coprime(&ring, &xy_mono, &x2_mono),
                "xy and x^2 should not be coprime");
    }

    #[test]
    fn test_f4_katsura_3() {
        // Test Katsura-3 system for correctness with GM criteria
        use crate::BN254_FR;

        let field = &*BN254_FR;
        let n = 3;
        let degree_cfg = DegreeCfg::new(100).with_precompute(50);
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n, degree_cfg, (0, 0), Global
        );

        // Create monomials for variables
        let mut vars = Vec::new();
        for i in 0..n {
            let mut exponents = vec![0; n];
            exponents[i] = 1;
            vars.push(poly_ring.from_terms([(
                poly_ring.base_ring().one(),
                poly_ring.create_monomial(exponents.into_iter())
            )].into_iter()));
        }

        // Build Katsura-3 system
        let mut system = Vec::new();

        // u_0 + 2*sum u_i = 1
        let mut eq0 = poly_ring.clone_el(&vars[0]);
        for i in 1..n {
            let two = poly_ring.base_ring().int_hom().map(2);
            let term = poly_ring.inclusion().map(two);
            let scaled = poly_ring.mul_ref(&term, &vars[i]);
            eq0 = poly_ring.add_ref(&eq0, &scaled);
        }
        let one = poly_ring.one();
        system.push(poly_ring.add_ref(&eq0, &poly_ring.negate(one)));

        // Katsura equations
        for k in 0..(n-1) {
            let mut poly = poly_ring.zero();
            for i in 0..n {
                for j in 0..n {
                    let sum_abs = if i >= j { i - j } else { j - i };
                    if sum_abs == k {
                        let term = poly_ring.mul_ref(&vars[i], &vars[j]);
                        poly = poly_ring.add_ref(&poly, &term);
                    }
                }
            }
            poly = poly_ring.sub_ref(&poly, &vars[k]);
            system.push(poly);
        }

        // Run F4
        let gb = f4_simple(&poly_ring, system, DegRevLex);

        // Should produce a valid basis (Katsura-3 typically gives 4 elements)
        assert!(!gb.is_empty());
        assert!(gb.len() >= 3, "Katsura-3 should produce at least 3 basis elements");
        assert!(gb.len() <= 5, "Katsura-3 should produce at most 5 basis elements");
    }
}
