//! Matrix construction for F4 algorithm
//!
//! This module provides the Macaulay matrix construction that is central to the F4 algorithm.
//! Unlike the F4-style Buchberger which reduces polynomials individually, true F4 builds
//! a large sparse matrix from S-polynomials and basis multiples, then performs a single
//! Gaussian elimination pass.

use feanor_math::ring::*;
use feanor_math::rings::multivariate::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Thread-safe pivot rows using atomic operations
///
/// This structure wraps a Vec<AtomicUsize> to enable lock-free concurrent pivot marking.
/// Encoding: 0 = None (no pivot), row_idx+1 = Some(row_idx)
#[derive(Debug)]
pub struct AtomicPivotRows {
    /// Internal storage: 0 = None, row+1 = Some(row)
    inner: Vec<AtomicUsize>,
}

impl AtomicPivotRows {
    /// Create a new AtomicPivotRows with given size, all initialized to None
    pub fn new(size: usize) -> Self {
        AtomicPivotRows {
            inner: (0..size).map(|_| AtomicUsize::new(0)).collect(),
        }
    }

    /// Check if a column has a pivot
    pub fn has_pivot(&self, col: usize) -> bool {
        self.inner[col].load(Ordering::Acquire) != 0
    }

    /// Get the pivot row for a column (if it exists)
    pub fn get_pivot_row(&self, col: usize) -> Option<usize> {
        let val = self.inner[col].load(Ordering::Acquire);
        if val == 0 {
            None
        } else {
            Some(val - 1)
        }
    }

    /// Try to mark a column as having a pivot at the given row
    ///
    /// Uses compare-and-swap to atomically mark the pivot.
    /// Returns true if successful, false if another thread already marked it.
    pub fn try_mark_pivot(&self, col: usize, row: usize) -> bool {
        self.inner[col]
            .compare_exchange(0, row + 1, Ordering::AcqRel, Ordering::Acquire)
            .is_ok()
    }

    /// Mark a pivot unconditionally (for single-threaded code)
    pub fn mark_pivot(&self, col: usize, row: usize) {
        self.inner[col].store(row + 1, Ordering::Release);
    }

    /// Get the number of columns
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Extend with more columns (initialized to None)
    pub fn extend(&mut self, count: usize) {
        self.inner.extend((0..count).map(|_| AtomicUsize::new(0)));
    }

    /// Convert to Vec<Option<usize>> for compatibility with single-threaded code
    pub fn to_vec(&self) -> Vec<Option<usize>> {
        self.inner
            .iter()
            .map(|atomic| {
                let val = atomic.load(Ordering::Acquire);
                if val == 0 {
                    None
                } else {
                    Some(val - 1)
                }
            })
            .collect()
    }

    /// Create from Vec<Option<usize>> (for migration/testing)
    pub fn from_vec(vec: Vec<Option<usize>>) -> Self {
        AtomicPivotRows {
            inner: vec
                .into_iter()
                .map(|opt| AtomicUsize::new(opt.map_or(0, |r| r + 1)))
                .collect(),
        }
    }

    /// Remap pivot rows according to a column mapping
    pub fn remap(&self, old_to_new: &[usize]) -> Self {
        let mut new_inner: Vec<AtomicUsize> = (0..self.len())
            .map(|_| AtomicUsize::new(0))
            .collect();
        for (old_col, &new_col) in old_to_new.iter().enumerate() {
            let val = self.inner[old_col].load(Ordering::Acquire);
            new_inner[new_col].store(val, Ordering::Release);
        }
        AtomicPivotRows { inner: new_inner }
    }
}

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

/// Type of row in the Macaulay matrix for AB/CD semantics
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowType {
    /// Reducer row (rr): first generator per LCM group, forms A-block reducers
    Rr,
    /// To-reduce row (tr): remaining generators per LCM group, reduced by rr rows
    Tr,
    /// Dedicated reducer: added explicitly for pivot columns
    Dedicated,
    /// Closure reducer: added during symbolic closure
    Closure,
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

    /// Cached expanded exponents for each column (same indexing as `col_to_monomial`)
    pub col_exponents: Vec<Vec<usize>>,

    /// Mapping from monomial to column index
    pub monomial_to_col: HashMap<Vec<usize>, usize>,

    /// Pivot columns (columns with known reducers)
    /// pivot_rows[col] = Some(row_index) if column `col` has a pivot in `row_index`
    pub pivot_rows: AtomicPivotRows,

    /// Row types for AB/CD semantics (tracks which rows are rr vs tr)
    pub row_types: Vec<RowType>,

    /// A-block width (number of pivot LCM columns) after SHT reordering
    /// Set after reorder_columns_with_sht; used for AB/CD phase boundary
    pub ncl: Option<usize>,
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
            col_exponents: Vec::new(),
            monomial_to_col: HashMap::new(),
            pivot_rows: AtomicPivotRows::new(0),
            row_types: Vec::new(),
            ncl: None,
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
        self.col_exponents.push(expanded.clone());
        self.monomial_to_col.insert(expanded, col);
        self.pivot_rows.extend(1);
        col
    }

    /// Get cached expanded exponents for a column
    #[inline]
    pub fn col_exp(&self, col: usize) -> &Vec<usize> {
        &self.col_exponents[col]
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

    /// Add a row to the matrix with specified row type
    pub fn add_row(&mut self, row: SparseRow<PolyCoeff<P>>, row_type: RowType) -> usize {
        let row_idx = self.rows.len();
        self.rows.push(row);
        self.row_types.push(row_type);
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
            self.pivot_rows.mark_pivot(col, row);
        }
    }

    /// Check if a column has a pivot
    pub fn has_pivot(&self, col: usize) -> bool {
        col < self.pivot_rows.len() && self.pivot_rows.has_pivot(col)
    }

    /// Get the row index of a pivot column
    pub fn get_pivot_row(&self, col: usize) -> Option<usize> {
        if col < self.pivot_rows.len() {
            self.pivot_rows.get_pivot_row(col)
        } else {
            None
        }
    }

    /// Reorder columns using an SHT-like index: entries with idx==2 (pivots)
    /// come first, then idx==1, then the remaining columns in their current order.
    /// The `sht_idx` map uses expanded monomials as keys.
    pub fn reorder_columns_with_sht(&mut self, sht_idx: &HashMap<Vec<usize>, u8>) {
        let ncols = self.col_to_monomial.len();
        if ncols == 0 || sht_idx.is_empty() {
            return;
        }

        // Build lists of columns by idx class
        let mut piv: Vec<usize> = Vec::new();
        let mut red: Vec<usize> = Vec::new();
        let mut seen: Vec<bool> = vec![false; ncols];

        for (exp, idx) in sht_idx {
            if let Some(&col) = self.monomial_to_col.get(exp) {
                match *idx {
                    2 => { piv.push(col); seen[col] = true; },
                    1 => { red.push(col); seen[col] = true; },
                    _ => {},
                }
            }
        }

        // Sort pivot and reducer columns deterministically to eliminate HashMap iteration order effects
        // This ensures reproducible column ordering across runs, which is critical for debugging
        // Sort by (degree, expanded monomial) to match msolve's deterministic ordering
        piv.sort_unstable_by(|&a, &b| {
            let mono_a = &self.col_to_monomial[a];
            let mono_b = &self.col_to_monomial[b];
            let deg_a = self.ring.monomial_deg(mono_a);
            let deg_b = self.ring.monomial_deg(mono_b);
            deg_a.cmp(&deg_b).then_with(|| {
                self.col_exponents[a].cmp(&self.col_exponents[b])
            })
        });

        red.sort_unstable_by(|&a, &b| {
            let mono_a = &self.col_to_monomial[a];
            let mono_b = &self.col_to_monomial[b];
            let deg_a = self.ring.monomial_deg(mono_a);
            let deg_b = self.ring.monomial_deg(mono_b);
            deg_a.cmp(&deg_b).then_with(|| {
                self.col_exponents[a].cmp(&self.col_exponents[b])
            })
        });

        let mut tails: Vec<usize> = Vec::new();
        for col in 0..ncols {
            if !seen[col] { tails.push(col); }
        }

        // Build new order: pivots, reducers, then tails
        let mut new_order: Vec<usize> = Vec::with_capacity(ncols);
        new_order.extend_from_slice(&piv);
        new_order.extend_from_slice(&red);
        new_order.extend_from_slice(&tails);

        // Map old -> new
        let mut old_to_new = vec![0usize; ncols];
        for (new_col, &old_col) in new_order.iter().enumerate() {
            old_to_new[old_col] = new_col;
        }

        // Remap rows and re-sort by new column indices
        // Note: sorting is necessary because column remapping is not order-preserving
        // (pivots→first, reducers→middle, tails→last reordering changes the order)
        for row in &mut self.rows {
            for (col, _) in &mut row.entries {
                *col = old_to_new[*col];
            }
            // Use sort_unstable_by_key for better performance (stability not needed)
            row.entries.sort_unstable_by_key(|(col, _)| *col);
        }

        // Reorder col_to_monomial/col_exponents and rebuild monomial_to_col
        let mut new_cols: Vec<PolyMonomial<P>> = Vec::with_capacity(ncols);
        let mut new_exps: Vec<Vec<usize>> = Vec::with_capacity(ncols);
        for &old_col in &new_order {
            new_cols.push(self.ring.clone_monomial(&self.col_to_monomial[old_col]));
            new_exps.push(self.col_exponents[old_col].clone());
        }
        self.col_to_monomial = new_cols;
        self.col_exponents = new_exps;

        self.monomial_to_col.clear();
        self.monomial_to_col.reserve(ncols);
        for (new_col, em) in self.col_exponents.iter().enumerate() {
            self.monomial_to_col.insert(em.clone(), new_col);
        }

        // Remap pivot_rows
        self.pivot_rows = self.pivot_rows.remap(&old_to_new);
    }

    /// Reorder columns so pivot monomials come first
    ///
    /// This implements the A|B split from msolve: columns corresponding to S-poly LCMs
    /// (intended pivot columns) are moved to the left, reducing reducer pressure and
    /// improving elimination efficiency.
    pub fn reorder_columns_pivot_first(&mut self, pivot_monomials: &[PolyMonomial<P>]) {
        let num_cols = self.num_cols();
        if num_cols == 0 {
            return;
        }

        // Find column indices for pivot monomials
        let mut pivot_cols = Vec::new();
        for pivot_mono in pivot_monomials {
            let expanded = self.ring.expand_monomial(pivot_mono);
            if let Some(&col) = self.monomial_to_col.get(&expanded) {
                pivot_cols.push(col);
            }
        }

        // Remove duplicates and sort for stable ordering
        pivot_cols.sort_unstable();
        pivot_cols.dedup();

        // Build new column order: pivot columns first, then non-pivot columns
        // Use HashSet for O(1) lookup instead of Vec::contains which is O(n)
        let pivot_set: std::collections::HashSet<usize> = pivot_cols.iter().copied().collect();
        let mut new_order = pivot_cols.clone();
        for col in 0..num_cols {
            if !pivot_set.contains(&col) {
                new_order.push(col);
            }
        }

        // Build old_to_new mapping
        let mut old_to_new = vec![0; num_cols];
        for (new_col, &old_col) in new_order.iter().enumerate() {
            old_to_new[old_col] = new_col;
        }

        // Remap all rows
        for row in &mut self.rows {
            for entry in &mut row.entries {
                entry.0 = old_to_new[entry.0];
            }
            // Re-sort after remapping
            row.sort();
        }

        // Reorder col_to_monomial
        let old_col_to_monomial = std::mem::take(&mut self.col_to_monomial);
        self.col_to_monomial = new_order
            .iter()
            .map(|&old_col| self.ring.clone_monomial(&old_col_to_monomial[old_col]))
            .collect();

        // Rebuild monomial_to_col using cached exponents; also keep col_exponents aligned
        let old_exps = std::mem::take(&mut self.col_exponents);
        self.col_exponents = new_order.iter().map(|&old_col| old_exps[old_col].clone()).collect();
        self.monomial_to_col.clear();
        for (new_col, em) in self.col_exponents.iter().enumerate() {
            self.monomial_to_col.insert(em.clone(), new_col);
        }

        // Reorder pivot_rows
        self.pivot_rows = self.pivot_rows.remap(&old_to_new);
    }

    /// Reorder columns so that the provided pivot monomials come first.
    ///
    /// This mimics the A|B split used by optimized F4 implementations: columns
    /// corresponding to LCMs of selected S-polynomials are placed first, which
    /// improves reducer effectiveness and elimination locality.
    pub fn reorder_columns_with_pivots(&mut self, pivot_monomials: &[PolyMonomial<P>]) {
        if self.col_to_monomial.is_empty() || pivot_monomials.is_empty() {
            return;
        }

        // Build a set of expanded pivot monomials for quick membership tests
        let mut pivot_set: HashSet<Vec<usize>> = HashSet::with_capacity(pivot_monomials.len());
        for m in pivot_monomials {
            pivot_set.insert(self.ring.expand_monomial(m));
        }

        // Determine new column order: first pivot columns, then the rest
        let ncols = self.col_to_monomial.len();
        let mut pivots: Vec<usize> = Vec::new();
        let mut tails: Vec<usize> = Vec::new();
        for col in 0..ncols {
            let em = &self.col_exponents[col];
            if pivot_set.contains(em) {
                pivots.push(col);
            } else {
                tails.push(col);
            }
        }

        // If nothing to reorder, exit early
        if pivots.is_empty() {
            return;
        }

        let mut old_to_new = vec![0usize; ncols];
        let mut new_order: Vec<usize> = Vec::with_capacity(ncols);
        for (i, &c) in pivots.iter().enumerate() {
            old_to_new[c] = i;
            new_order.push(c);
        }
        for (k, &c) in tails.iter().enumerate() {
            let idx = pivots.len() + k;
            old_to_new[c] = idx;
            new_order.push(c);
        }

        // Remap rows and re-sort by new column indices
        for row in &mut self.rows {
            for (col, _) in &mut row.entries {
                *col = old_to_new[*col];
            }
            // Use sort_unstable_by_key for better performance (stability not needed)
            row.entries.sort_unstable_by_key(|(col, _)| *col);
        }

        // Reorder col_to_monomial/col_exponents and rebuild monomial_to_col
        let mut new_cols: Vec<PolyMonomial<P>> = Vec::with_capacity(ncols);
        let mut new_exps: Vec<Vec<usize>> = Vec::with_capacity(ncols);
        for &old_col in &new_order {
            new_cols.push(self.ring.clone_monomial(&self.col_to_monomial[old_col]));
            new_exps.push(self.col_exponents[old_col].clone());
        }
        self.col_to_monomial = new_cols;
        self.col_exponents = new_exps;

        self.monomial_to_col.clear();
        self.monomial_to_col.reserve(ncols);
        for (new_col, em) in self.col_exponents.iter().enumerate() {
            self.monomial_to_col.insert(em.clone(), new_col);
        }

        // Remap pivot_rows to the new order
        self.pivot_rows = self.pivot_rows.remap(&old_to_new);
    }

    /// Reorder columns so that the provided pivot column indices come first.
    ///
    /// This variant avoids monomial expansion and works directly with known
    /// column indices of pivot (LCM) monomials.
    pub fn reorder_columns_with_pivot_indices(&mut self, pivot_indices: &[usize]) {
        let ncols = self.col_to_monomial.len();
        if ncols == 0 || pivot_indices.is_empty() {
            return;
        }

        // Build a boolean marker array and an ordered unique list of pivots
        let mut is_pivot = vec![false; ncols];
        let mut piv: Vec<usize> = Vec::new();
        piv.reserve(pivot_indices.len());
        for &c in pivot_indices {
            if c < ncols && !is_pivot[c] {
                is_pivot[c] = true;
                piv.push(c);
            }
        }
        if piv.is_empty() {
            return;
        }

        // Build new column order: pivots (in given order) followed by tails
        let mut new_order: Vec<usize> = Vec::with_capacity(ncols);
        new_order.extend_from_slice(&piv);
        for col in 0..ncols {
            if !is_pivot[col] {
                new_order.push(col);
            }
        }

        // Map old -> new
        let mut old_to_new = vec![0usize; ncols];
        for (new_col, &old_col) in new_order.iter().enumerate() {
            old_to_new[old_col] = new_col;
        }

        // Remap rows
        for row in &mut self.rows {
            for (col, _) in &mut row.entries {
                *col = old_to_new[*col];
            }
            row.sort();
        }

        // Reorder col_to_monomial/col_exponents and rebuild monomial_to_col
        let mut new_cols: Vec<PolyMonomial<P>> = Vec::with_capacity(ncols);
        let mut new_exps: Vec<Vec<usize>> = Vec::with_capacity(ncols);
        for &old_col in &new_order {
            new_cols.push(self.ring.clone_monomial(&self.col_to_monomial[old_col]));
            new_exps.push(self.col_exponents[old_col].clone());
        }
        self.col_to_monomial = new_cols;
        self.col_exponents = new_exps;

        self.monomial_to_col.clear();
        self.monomial_to_col.reserve(ncols);
        for (new_col, em) in self.col_exponents.iter().enumerate() {
            self.monomial_to_col.insert(em.clone(), new_col);
        }

        // Remap pivot_rows
        self.pivot_rows = self.pivot_rows.remap(&old_to_new);
    }

    /// Sort rows in the given range by pivot column (ascending) then density (ascending).
    ///
    /// This improves cache locality and reduction efficiency during Gaussian elimination.
    /// Rows are sorted by:
    /// 1. Pivot column (leftmost non-zero column) - ascending
    /// 2. Density (number of entries) - ascending
    ///
    /// Zero rows (empty rows) are placed at the end of the range.
    /// The corresponding row_types are permuted along with the rows.
    pub fn sort_rows_by_pivot_and_density(&mut self, range: std::ops::Range<usize>) {
        if range.start >= range.end || range.end > self.rows.len() {
            return; // Empty or invalid range
        }

        // Build index array for stable sorting
        let len = range.end - range.start;
        let mut indices: Vec<usize> = (0..len).collect();

        // Compute sorting keys once: (pivot_col, density) for each row
        let keys: Vec<(Option<usize>, usize)> = (range.start..range.end)
            .map(|i| {
                let pivot = self.rows[i].pivot();
                let density = self.rows[i].entries.len();
                (pivot, density)
            })
            .collect();

        // Sort indices by keys: None pivots (zero rows) last, then by (pivot, density)
        indices.sort_by_key(|&idx| {
            let (pivot, density) = keys[idx];
            match pivot {
                Some(p) => (0, p, density), // Non-zero rows: sort by pivot then density
                None => (1, usize::MAX, density), // Zero rows: place at end
            }
        });

        // Check if already sorted (optimization to avoid unnecessary work)
        if indices.iter().enumerate().all(|(i, &idx)| i == idx) {
            return;
        }

        // Apply permutation in-place using cycle decomposition
        let mut done = vec![false; len];
        for start in 0..len {
            if done[start] {
                continue;
            }

            let mut current = start;
            while !done[current] {
                done[current] = true;
                let next = indices[current];
                if next != start && !done[next] {
                    // Swap rows and row_types
                    self.rows.swap(range.start + current, range.start + next);
                    self.row_types.swap(range.start + current, range.start + next);
                }
                current = next;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feanor_math::homomorphism::Homomorphism;
    use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
    use feanor_math::rings::multivariate::DegRevLex;
    use crate::BN254_FR;

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
        let base = &*BN254_FR;
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
        let base = &*BN254_FR;
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
