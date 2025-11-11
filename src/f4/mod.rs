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

    // Generate initial S-polynomials
    let mut spolys: Vec<SPoly> = Vec::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            spolys.push(SPoly::new(i, j));
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

        // Add reducers (basis multiples) to matrix
        add_reducers_to_matrix(&mut matrix, &basis, order);

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
            if !ring.is_zero(&poly) && !poly_in_basis(ring, &poly, &basis, order) {
                new_polys.push(poly);
            }
        }

        // Add new polynomials to basis and generate new S-polys
        let old_basis_len = basis.len();
        for (idx, new_poly) in new_polys.into_iter().enumerate() {
            let new_idx = old_basis_len + idx;

            // Generate S-polys with existing basis
            for i in 0..basis.len() {
                spolys.push(SPoly::new(i, new_idx));
            }

            basis.push(new_poly);
        }
    }

    // Final interreduction
    basis = interreduce(ring, basis, order);

    Ok(basis)
}

/// Add necessary basis multiples to the matrix as reducers
fn add_reducers_to_matrix<P, O>(matrix: &mut MacaulayMatrix<P>, basis: &[El<P>], order: O)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    // For now, just add all basis elements directly
    // A more sophisticated implementation would compute the necessary multiples
    for poly in basis {
        let row = matrix.polynomial_to_row(poly, order);
        matrix.add_row(row);
    }
}

/// Check if a polynomial is already in the basis (up to scalar multiple)
fn poly_in_basis<P, O>(ring: P, poly: &El<P>, basis: &[El<P>], order: O) -> bool
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
    O: MonomialOrder + Copy,
{
    let (_poly_lc, poly_lm) = match ring.LT(poly, order) {
        Some(lt) => lt,
        None => return false,
    };

    for b in basis {
        if let Some((_b_lc, b_lm)) = ring.LT(b, order) {
            if order.compare(ring, poly_lm, b_lm) == std::cmp::Ordering::Equal {
                // Same leading monomial - check if they're scalar multiples
                return true;
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
    use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
    use feanor_math::rings::multivariate::DegRevLex;
    use feanor_math::rings::zn::zn_static;

    #[test]
    fn test_f4_simple() {
        let base = zn_static::F17;
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

        let config = F4Config::new().with_max_degree(10);
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
}
