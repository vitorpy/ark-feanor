//! Matrix construction for F4 algorithm
//!
//! This module provides the Macaulay matrix construction that is central to the F4 algorithm.
//! Unlike the F4-style Buchberger which reduces polynomials individually, true F4 builds
//! a large sparse matrix from S-polynomials and basis multiples, then performs a single
//! Gaussian elimination pass.

use feanor_math::ring::*;
use feanor_math::rings::multivariate::*;
use std::collections::HashMap;

/// A sparse row in the Macaulay matrix
///
/// Each row represents a polynomial, stored as (column_index, coefficient) pairs.
/// Columns are sorted by index for efficient operations.
#[derive(Debug)]
pub struct SparseRow<F> {
    /// (column_index, coefficient) pairs, sorted by column_index
    pub entries: Vec<(usize, F)>,
}

impl<F> SparseRow<F> {
    /// Create a new empty sparse row
    pub fn new() -> Self {
        SparseRow { entries: Vec::new() }
    }

    /// Create a sparse row with given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        SparseRow {
            entries: Vec::with_capacity(capacity),
        }
    }

    /// Add an entry to the row (must maintain sorted order)
    pub fn push(&mut self, col: usize, coeff: F) {
        self.entries.push((col, coeff));
    }

    /// Sort entries by column index
    pub fn sort(&mut self) {
        self.entries.sort_by_key(|(col, _)| *col);
    }

    /// Find the pivot (leftmost non-zero column)
    pub fn pivot(&self) -> Option<usize> {
        self.entries.first().map(|(col, _)| *col)
    }

    /// Check if row is zero (empty)
    pub fn is_zero(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Macaulay matrix for F4 algorithm
///
/// The matrix is stored in sparse row format. Each row corresponds to either:
/// - An S-polynomial that needs to be reduced
/// - A multiple of a basis element (used as a reducer)
///
/// Columns correspond to monomials in a fixed order.
pub struct MacaulayMatrix<P>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    /// The polynomial ring
    pub ring: P,

    /// Sparse rows of the matrix
    pub rows: Vec<SparseRow<PolyCoeff<P>>>,

    /// Mapping from column index to monomial
    pub col_to_monomial: Vec<PolyMonomial<P>>,

    /// Mapping from monomial to column index
    pub monomial_to_col: HashMap<Vec<usize>, usize>,

    /// Pivot columns (columns with known reducers)
    /// pivot_rows[col] = Some(row_index) if column `col` has a pivot in `row_index`
    pub pivot_rows: Vec<Option<usize>>,
}

impl<P> MacaulayMatrix<P>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
{
    /// Create a new empty Macaulay matrix
    pub fn new(ring: P) -> Self {
        MacaulayMatrix {
            ring,
            rows: Vec::new(),
            col_to_monomial: Vec::new(),
            monomial_to_col: HashMap::new(),
            pivot_rows: Vec::new(),
        }
    }

    /// Get or create a column index for a given monomial
    pub fn get_or_create_column(&mut self, monomial: &PolyMonomial<P>) -> usize {
        let expanded: Vec<usize> = self.ring.expand_monomial(monomial);

        if let Some(&col) = self.monomial_to_col.get(&expanded) {
            return col;
        }

        let col = self.col_to_monomial.len();
        self.col_to_monomial.push(self.ring.clone_monomial(monomial));
        self.monomial_to_col.insert(expanded, col);
        self.pivot_rows.push(None);
        col
    }

    /// Convert a polynomial to a sparse row
    ///
    /// Each term (coeff, monomial) becomes an entry (column_index, coeff) in the row.
    pub fn polynomial_to_row<O>(&mut self, poly: &El<P>, _order: O) -> SparseRow<PolyCoeff<P>>
    where
        O: MonomialOrder + Copy,
    {
        // Collect terms first to avoid borrow checker issues
        let terms: Vec<_> = self.ring.terms(poly)
            .map(|(c, m)| (self.ring.base_ring().clone_el(c), self.ring.clone_monomial(m)))
            .collect();

        let mut row = SparseRow::with_capacity(terms.len());

        for (coeff, monomial) in terms {
            let col = self.get_or_create_column(&monomial);
            row.push(col, coeff);
        }

        // Sort by column index (monomial order)
        row.sort();
        row
    }

    /// Add a row to the matrix
    pub fn add_row(&mut self, row: SparseRow<PolyCoeff<P>>) -> usize {
        let row_idx = self.rows.len();
        self.rows.push(row);
        row_idx
    }

    /// Convert a sparse row back to a polynomial
    pub fn row_to_polynomial(&self, row: &SparseRow<PolyCoeff<P>>) -> El<P> {
        let terms: Vec<_> = row
            .entries
            .iter()
            .map(|(col, coeff)| {
                (
                    self.ring.base_ring().clone_el(coeff),
                    self.ring.clone_monomial(&self.col_to_monomial[*col]),
                )
            })
            .collect();

        self.ring.from_terms(terms.into_iter())
    }

    /// Get the number of rows
    pub fn num_rows(&self) -> usize {
        self.rows.len()
    }

    /// Get the number of columns
    pub fn num_cols(&self) -> usize {
        self.col_to_monomial.len()
    }

    /// Mark a column as having a pivot in a given row
    pub fn mark_pivot(&mut self, col: usize, row: usize) {
        if col < self.pivot_rows.len() {
            self.pivot_rows[col] = Some(row);
        }
    }

    /// Check if a column has a pivot
    pub fn has_pivot(&self, col: usize) -> bool {
        col < self.pivot_rows.len() && self.pivot_rows[col].is_some()
    }

    /// Get the row index of a pivot column
    pub fn get_pivot_row(&self, col: usize) -> Option<usize> {
        self.pivot_rows.get(col).and_then(|&r| r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feanor_math::homomorphism::Homomorphism;
    use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
    use feanor_math::rings::multivariate::DegRevLex;
    use feanor_math::rings::zn::zn_static;

    #[test]
    fn test_sparse_row_basic() {
        let mut row = SparseRow::<i32>::new();
        row.push(5, 3);
        row.push(2, 7);
        row.push(10, -1);

        row.sort();

        assert_eq!(row.entries, vec![(2, 7), (5, 3), (10, -1)]);
        assert_eq!(row.pivot(), Some(2));
        assert!(!row.is_zero());
    }

    #[test]
    fn test_polynomial_to_row() {
        let base = zn_static::F17;
        let ring = MultivariatePolyRingImpl::new(base, 3);

        let mut matrix = MacaulayMatrix::new(&ring);

        // Create polynomial: 2x² + 3xy + 5y²
        let poly = ring.from_terms([
            (base.int_hom().map(2), ring.create_monomial([2, 0, 0])),
            (base.int_hom().map(3), ring.create_monomial([1, 1, 0])),
            (base.int_hom().map(5), ring.create_monomial([0, 2, 0])),
        ].into_iter());

        let row = matrix.polynomial_to_row(&poly, DegRevLex);

        assert_eq!(row.entries.len(), 3);
        assert!(!row.is_zero());

        // Convert back and check equality
        let poly2 = matrix.row_to_polynomial(&row);
        let diff = ring.sub(ring.clone_el(&poly), poly2);
        assert!(ring.is_zero(&diff));
    }

    #[test]
    fn test_column_management() {
        let base = zn_static::F17;
        let ring = MultivariatePolyRingImpl::new(base, 2);

        let mut matrix = MacaulayMatrix::new(&ring);

        let m1 = ring.create_monomial([1, 0]);
        let m2 = ring.create_monomial([0, 1]);
        let m3 = ring.create_monomial([1, 0]); // Same as m1

        let col1 = matrix.get_or_create_column(&m1);
        let col2 = matrix.get_or_create_column(&m2);
        let col3 = matrix.get_or_create_column(&m3);

        assert_eq!(col1, col3); // Same monomial -> same column
        assert_ne!(col1, col2); // Different monomials -> different columns
        assert_eq!(matrix.num_cols(), 2);
    }
}
