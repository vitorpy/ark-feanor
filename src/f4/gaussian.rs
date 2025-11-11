//! Sparse Gaussian elimination for F4 algorithm
//!
//! This module implements row-echelon form reduction for sparse matrices over fields.
//! Optimized for cryptographically large prime fields (BN254, BLS12-381) where field
//! operations are expensive.

use super::matrix::{MacaulayMatrix, SparseRow};
use feanor_math::divisibility::DivisibilityRingStore;
use feanor_math::field::Field;
use feanor_math::ring::*;
use feanor_math::rings::multivariate::*;

/// Reduce the Macaulay matrix to row echelon form
///
/// This performs Gaussian elimination on the matrix, reducing all rows.
/// The matrix is modified in place.
///
/// Returns the indices of rows with new pivots (i.e., non-zero rows after reduction)
pub fn reduce_matrix<P>(matrix: &mut MacaulayMatrix<P>) -> Vec<usize>
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
{
    let num_rows = matrix.num_rows();
    let mut new_pivot_rows = Vec::new();

    // Process each row
    for row_idx in 0..num_rows {
        // Reduce the current row by all known pivots
        reduce_row_by_pivots(matrix, row_idx);

        // Check if the row has a new pivot
        if let Some(pivot_col) = matrix.rows[row_idx].pivot() {
            // Normalize the row so the pivot coefficient is 1
            normalize_row(&matrix.ring, &mut matrix.rows[row_idx], pivot_col);

            // If this column doesn't have a pivot yet, mark it
            if !matrix.has_pivot(pivot_col) {
                matrix.mark_pivot(pivot_col, row_idx);
                new_pivot_rows.push(row_idx);
            }
        }
    }

    new_pivot_rows
}

/// Reduce a single row by all known pivot rows
fn reduce_row_by_pivots<P>(matrix: &mut MacaulayMatrix<P>, row_idx: usize)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
{
    let base_ring = matrix.ring.base_ring();

    loop {
        let pivot_col = match matrix.rows[row_idx].pivot() {
            Some(col) => col,
            None => break, // Row is zero
        };

        // Check if there's a pivot row for this column
        let pivot_row_idx = match matrix.get_pivot_row(pivot_col) {
            Some(idx) if idx != row_idx => idx,
            _ => break, // No reducer available or we'd be reducing by ourselves
        };

        // Clone the pivot row entries to avoid borrow checker issues
        let pivot_row_entries: Vec<(usize, PolyCoeff<P>)> = matrix.rows[pivot_row_idx]
            .entries
            .iter()
            .map(|(col, coeff)| (*col, base_ring.clone_el(coeff)))
            .collect();
        let pivot_row: SparseRow<PolyCoeff<P>> = SparseRow { entries: pivot_row_entries };

        // Reduce: row -= (row[pivot_col] / pivot_row[pivot_col]) * pivot_row
        let row_leading_coeff = base_ring.clone_el(&matrix.rows[row_idx].entries[0].1);
        let pivot_leading_coeff = base_ring.clone_el(&pivot_row.entries[0].1);

        // Compute multiplier: row_leading_coeff / pivot_leading_coeff
        let multiplier = base_ring.checked_div(&row_leading_coeff, &pivot_leading_coeff).unwrap();

        // Perform: row -= multiplier * pivot_row
        subtract_scaled_row(
            &base_ring,
            &mut matrix.rows[row_idx],
            &pivot_row,
            &multiplier,
        );
    }
}

/// Normalize a row so that the coefficient at `pivot_col` is 1
fn normalize_row<P>(ring: &P, row: &mut SparseRow<PolyCoeff<P>>, pivot_col: usize)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
{
    let base_ring = ring.base_ring();

    // Find the pivot coefficient
    let pivot_entry = row.entries.iter().find(|(col, _)| *col == pivot_col);
    if let Some((_, pivot_coeff)) = pivot_entry {
        let inv = base_ring.invert(pivot_coeff).unwrap();

        // Multiply all coefficients by the inverse
        for (_, coeff) in &mut row.entries {
            base_ring.mul_assign(coeff, base_ring.clone_el(&inv));
        }
    }
}

/// Compute: target -= multiplier * source
///
/// This is the core operation in Gaussian elimination. We subtract a scaled version
/// of the source row from the target row, canceling the leading term.
fn subtract_scaled_row<R>(
    base_ring: &R,
    target: &mut SparseRow<El<R>>,
    source: &SparseRow<El<R>>,
    multiplier: &El<R>,
) where
    R: RingStore + Copy,
    R::Type: Field,
{
    // Merge source into target
    let mut result = Vec::new();
    let mut target_iter = target.entries.iter();
    let mut source_iter = source.entries.iter();

    let mut target_next = target_iter.next();
    let mut source_next = source_iter.next();

    while target_next.is_some() || source_next.is_some() {
        match (target_next, source_next) {
            (Some((tc, tv)), Some((sc, sv))) => {
                if tc == sc {
                    // Same column: subtract
                    let mut new_val = base_ring.clone_el(tv);
                    let scaled = base_ring.mul_ref(multiplier, sv);
                    base_ring.sub_assign(&mut new_val, scaled);

                    if !base_ring.is_zero(&new_val) {
                        result.push((*tc, new_val));
                    }

                    target_next = target_iter.next();
                    source_next = source_iter.next();
                } else if tc < sc {
                    // Target column comes first
                    result.push((*tc, base_ring.clone_el(tv)));
                    target_next = target_iter.next();
                } else {
                    // Source column comes first
                    let mut new_val = base_ring.mul_ref(multiplier, sv);
                    base_ring.negate_inplace(&mut new_val);
                    if !base_ring.is_zero(&new_val) {
                        result.push((*sc, new_val));
                    }
                    source_next = source_iter.next();
                }
            }
            (Some((tc, tv)), None) => {
                result.push((*tc, base_ring.clone_el(tv)));
                target_next = target_iter.next();
            }
            (None, Some((sc, sv))) => {
                let mut new_val = base_ring.mul_ref(multiplier, sv);
                base_ring.negate_inplace(&mut new_val);
                if !base_ring.is_zero(&new_val) {
                    result.push((*sc, new_val));
                }
                source_next = source_iter.next();
            }
            (None, None) => unreachable!(),
        }
    }

    target.entries = result;
}

#[cfg(test)]
mod tests {
    use super::*;
    use feanor_math::homomorphism::Homomorphism;
    use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
    use feanor_math::rings::multivariate::DegRevLex;
    use feanor_math::rings::zn::zn_static;

    #[test]
    fn test_subtract_scaled_row() {
        let base = zn_static::F17;

        // target = [(0, 5), (2, 3), (5, 7)]
        let mut target = SparseRow {
            entries: vec![
                (0, base.int_hom().map(5)),
                (2, base.int_hom().map(3)),
                (5, base.int_hom().map(7)),
            ],
        };

        // source = [(0, 1), (3, 2), (5, 1)]
        let source = SparseRow {
            entries: vec![
                (0, base.int_hom().map(1)),
                (3, base.int_hom().map(2)),
                (5, base.int_hom().map(1)),
            ],
        };

        // multiplier = 5 (so we should cancel column 0)
        let multiplier = base.int_hom().map(5);

        subtract_scaled_row(&base, &mut target, &source, &multiplier);

        // Expected: column 0 cancels (5 - 5*1 = 0)
        //           column 2 stays: 3
        //           column 3 appears: -5*2 = -10 = 7 (mod 17)
        //           column 5: 7 - 5*1 = 2

        assert_eq!(target.entries.len(), 3);
        assert_eq!(target.entries[0], (2, base.int_hom().map(3)));
        assert_eq!(target.entries[1], (3, base.int_hom().map(7))); // -10 mod 17 = 7
        assert_eq!(target.entries[2], (5, base.int_hom().map(2)));
    }

    #[test]
    fn test_matrix_reduction() {
        let base = zn_static::F17;
        let ring = MultivariatePolyRingImpl::new(base, 2);

        let mut matrix = MacaulayMatrix::new(&ring);

        // Create two polynomials
        // p1: x + 2y
        let p1 = ring.from_terms([
            (base.int_hom().map(1), ring.create_monomial([1, 0])),
            (base.int_hom().map(2), ring.create_monomial([0, 1])),
        ].into_iter());

        // p2: 3x + 6y (which is 3*p1)
        let p2 = ring.from_terms([
            (base.int_hom().map(3), ring.create_monomial([1, 0])),
            (base.int_hom().map(6), ring.create_monomial([0, 1])),
        ].into_iter());

        let row1 = matrix.polynomial_to_row(&p1, DegRevLex);
        let row2 = matrix.polynomial_to_row(&p2, DegRevLex);

        matrix.add_row(row1);
        matrix.add_row(row2);

        // Reduce the matrix
        let new_pivots = reduce_matrix(&mut matrix);

        // Should have one pivot (row 0), row 1 should reduce to zero
        assert_eq!(new_pivots.len(), 1);
        assert!(matrix.rows[1].is_zero());
    }
}
