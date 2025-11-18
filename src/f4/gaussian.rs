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

/// AB phase reduction: reduce tr rows by rr rows using A-block columns only (columns < ncl)
///
/// This implements the AB phase of msolve's sparse linear algebra:
/// - Build rr_map: mapping from A-block columns to rr row indices
/// - For each tr row, reduce it by rr rows while its pivot is in the A-block
/// - Does NOT modify matrix.pivot_rows (rr rows are not global pivots)
/// - After this phase, tr rows should have leftmost column >= ncl or be zero
///
/// # Arguments
/// * `matrix` - The Macaulay matrix
/// * `rr_range` - Range of rr row indices (reducer rows)
/// * `tr_range` - Range of tr row indices (to-reduce rows)
/// * `ncl` - A-block width (columns 0..ncl are A-block, >= ncl are B/D blocks)
pub fn reduce_tr_by_rr_on_A<P>(
    matrix: &mut MacaulayMatrix<P>,
    rr_range: std::ops::Range<usize>,
    tr_range: std::ops::Range<usize>,
    ncl: usize,
)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
{
    let base_ring = matrix.ring.base_ring();

    // Step 1: Build rr_map: A-column -> rr_row_idx
    // Map each A-block column to the rr row whose leftmost entry is that column
    // Also cache the inverse of each rr row's leading coefficient for fast multiplier computation
    let mut rr_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut rr_inv_cache: std::collections::HashMap<usize, PolyCoeff<P>> = std::collections::HashMap::new();

    for rr_idx in rr_range {
        if let Some(pivot_col) = matrix.rows[rr_idx].pivot() {
            if pivot_col < ncl {
                // Only map A-block columns
                rr_map.insert(pivot_col, rr_idx);

                // Cache the inverse of the leading coefficient
                let rr_leading_coeff = &matrix.rows[rr_idx].entries[0].1;
                let rr_inv = base_ring.invert(rr_leading_coeff).unwrap();
                rr_inv_cache.insert(rr_idx, rr_inv);
            }
        }
    }

    // Step 2: For each tr row, reduce by rr_map while pivot < ncl
    for tr_idx in tr_range {
        let mut iterations = 0;
        let max_iterations = ncl + 10; // Safety bound: at most one reduction per A-column plus slack

        loop {
            // Safety guard: prevent infinite loops
            iterations += 1;
            if iterations > max_iterations {
                eprintln!("WARNING: AB reduction infinite loop detected for tr row {}", tr_idx);
                eprintln!("  Pivot column stuck, breaking. This may indicate a bug.");
                break;
            }

            // Get the current pivot of the tr row
            let pivot_col = match matrix.rows[tr_idx].pivot() {
                Some(col) => col,
                None => break, // Row is zero, done
            };

            // If pivot is outside A-block, we're done with this tr row
            if pivot_col >= ncl {
                break;
            }

            // Look up the rr row for this A-column
            let rr_idx = match rr_map.get(&pivot_col) {
                Some(&idx) => idx,
                None => break, // No rr reducer for this column, done
            };

            // Safety guard: prevent self-reduction (should never happen with correct ranges)
            if rr_idx == tr_idx {
                eprintln!("ERROR: AB reduction attempted self-reduction! rr_idx = tr_idx = {}", tr_idx);
                eprintln!("  This indicates rr_range incorrectly includes tr rows.");
                break;
            }

            // Compute multiplier: (tr_row leading coeff) * (cached inverse of rr_row leading coeff)
            // This is much faster than division as we pre-computed the inverse
            let tr_leading_coeff = base_ring.clone_el(&matrix.rows[tr_idx].entries[0].1);
            let rr_inv = rr_inv_cache.get(&rr_idx).unwrap();
            let multiplier = base_ring.mul_ref(&tr_leading_coeff, rr_inv);

            // Temporarily take out the tr row to allow immutable access to rr row
            let mut tr_row = std::mem::replace(&mut matrix.rows[tr_idx], SparseRow::new());

            // Reduce: tr_row -= multiplier * rr_row
            subtract_scaled_row(
                &base_ring,
                &mut tr_row,
                &matrix.rows[rr_idx],
                &multiplier,
            );

            // Put the tr row back
            matrix.rows[tr_idx] = tr_row;
        }

        // After AB reduction, tr row should have pivot >= ncl or be zero
        // This is enforced by the break conditions above
    }
}

/// CD phase reduction: reduce tr rows and discover pivots only in columns >= ncl (B/D blocks)
///
/// This is a simple one-pass CD phase that reduces tr rows by existing pivots
/// and discovers new pivots in the B/D blocks (columns >= ncl).
///
/// Returns: (new_pivot_rows, tr_new_pivots_count, anomaly_a_pivots)
pub fn reduce_matrix_tr<P>(
    matrix: &mut MacaulayMatrix<P>,
    tr_range: std::ops::Range<usize>,
    ncl: usize,
) -> (Vec<usize>, usize, usize)
where
    P: RingStore + Copy,
    P::Type: MultivariatePolyRing,
    <<P::Type as RingExtension>::BaseRing as RingStore>::Type: Field,
{
    let mut new_pivot_rows = Vec::new();
    let mut anomaly_a_pivots = 0;

    // Cache for pivot row leading coefficient inverses (speeds up multiplier computation)
    let mut pivot_inv_cache: std::collections::HashMap<usize, PolyCoeff<P>> = std::collections::HashMap::new();

    // Process each tr row
    for row_idx in tr_range {
        // Reduce the current row by all known pivots (including rr pivots in A-block)
        reduce_row_by_pivots(matrix, row_idx, &mut pivot_inv_cache);

        // Check if the row has a pivot after reduction
        if let Some(pivot_col) = matrix.rows[row_idx].pivot() {
            // Only discover pivots in B/D blocks (columns >= ncl)
            if pivot_col >= ncl {
                // Normalize the row so the pivot coefficient is 1
                normalize_row(&matrix.ring, &mut matrix.rows[row_idx], pivot_col);

                // If this column doesn't have a pivot yet, mark it
                if !matrix.has_pivot(pivot_col) {
                    matrix.mark_pivot(pivot_col, row_idx);
                    new_pivot_rows.push(row_idx);
                }
            } else {
                // If pivot_col < ncl, it means AB reduction didn't fully reduce this row
                // This is an anomaly
                anomaly_a_pivots += 1;
            }
        }
    }

    let tr_new_pivots_count = new_pivot_rows.len();
    (new_pivot_rows, tr_new_pivots_count, anomaly_a_pivots)
}

/// CD phase elimination: reduce tr rows and discover pivots only in columns >= ncl (B/D blocks)
///
/// This implements the CD phase of msolve's sparse linear algebra:
/// - Operates only on tr rows (ignores rr rows for pivot discovery)
/// - Discovers pivots only in columns >= ncl (B/D blocks)
/// - Uses existing pivots (from rr rows in A-block) for reduction
/// - Returns indices of tr rows with new pivots and counts for diagnostics
///
/// # Arguments
/// * `matrix` - The Macaulay matrix
/// * `tr_range` - Range of tr row indices (to-reduce rows)
/// * `ncl` - A-block width (skip pivot discovery for columns < ncl)
///
/// # Returns
/// Tuple of (pivot_rows, tr_new_pivots_count, anomaly_a_pivots_count)
/// - pivot_rows: Vector of tr row indices that have new pivots in B/D blocks
/// - tr_new_pivots_count: Count of new pivots discovered in tr rows (should equal pivot_rows.len())
/// - anomaly_a_pivots_count: Count of tr rows that still have pivot < ncl (should be 0)
/// Reduce a single row by all known pivot rows
fn reduce_row_by_pivots<P>(
    matrix: &mut MacaulayMatrix<P>,
    row_idx: usize,
    pivot_inv_cache: &mut std::collections::HashMap<usize, PolyCoeff<P>>,
)
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

        // Compute multiplier using cached inverse or compute and cache it
        let row_leading_coeff = base_ring.clone_el(&matrix.rows[row_idx].entries[0].1);
        let pivot_inv = pivot_inv_cache.entry(pivot_row_idx).or_insert_with(|| {
            let pivot_leading_coeff = &matrix.rows[pivot_row_idx].entries[0].1;
            base_ring.invert(pivot_leading_coeff).unwrap()
        });
        let multiplier = base_ring.mul_ref(&row_leading_coeff, pivot_inv);

        // Temporarily replace the target row with empty to allow immutable access to pivot row
        let mut row = std::mem::replace(&mut matrix.rows[row_idx], SparseRow::new());

        // Reduce: row -= multiplier * pivot_row
        subtract_scaled_row(
            &base_ring,
            &mut row,
            &matrix.rows[pivot_row_idx],
            &multiplier,
        );

        // Put the row back
        matrix.rows[row_idx] = row;
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
    // Reserve capacity to avoid reallocations during merge
    let mut result = Vec::with_capacity(target.entries.len() + source.entries.len());
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
}
