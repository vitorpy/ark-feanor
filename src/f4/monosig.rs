//! Monomial signature (divmask) for fast divisibility prefiltering.
//!
//! A MonSig is a lightweight summary of a monomial that allows us to quickly reject
//! impossible divisibility candidates before doing expensive per-variable checks and
//! monomial_div operations.
//!
//! The prefilter is based on three conservative tests:
//! 1. Degree bound: deg(LM) <= deg(target)
//! 2. Presence subset: vars(LM) ⊆ vars(target)
//! 3. Bucket upper bound: for each var, min(exp_LM, 15) <= min(exp_target, 15)
//!
//! If any test fails, we can immediately reject the candidate without calling monomial_div.

/// Monomial signature for fast divisibility prefiltering.
#[derive(Debug, Clone)]
pub struct MonSig {
    /// Total degree of the monomial
    pub deg: u16,
    /// Bitset marking which variables have non-zero exponents
    pub presence: PresenceMask,
    /// Optional nibble buckets: 4-bit caps per variable (min(exp, 15))
    pub buckets: Option<Vec<u64>>,
}

/// Bitset for tracking which variables appear in a monomial.
/// Uses compact u64 representation for small variable counts.
#[derive(Debug, Clone)]
pub enum PresenceMask {
    /// For polynomials with <= 64 variables
    U64(u64),
    /// For polynomials with > 64 variables
    Chunks(Vec<u64>),
}

impl MonSig {
    /// Create a MonSig from an exponent vector.
    ///
    /// # Arguments
    /// * `exp` - Exponent vector where exp[i] is the exponent of variable i
    ///
    /// # Returns
    /// A new MonSig with degree, presence mask, and nibble buckets computed
    pub fn from_exponents(exp: &[usize]) -> Self {
        let deg = exp.iter().sum::<usize>() as u16;
        let presence = PresenceMask::from_exponents(exp);
        let buckets = Some(pack_nibbles(exp));

        MonSig { deg, presence, buckets }
    }

    /// Check if this signature's presence mask is a subset of another's.
    /// Returns true if every variable present in self is also present in other.
    ///
    /// This is a necessary (but not sufficient) condition for divisibility:
    /// if LM | target, then vars(LM) ⊆ vars(target).
    #[inline]
    pub fn presence_subset(&self, other: &Self) -> bool {
        self.presence.is_subset(&other.presence)
    }

    /// Check if this signature's bucket values are all <= the other's.
    /// Returns true if for every variable i, bucket_self[i] <= bucket_other[i].
    ///
    /// This is a conservative approximation: if the real exponent self.exp[i] > 15,
    /// we cap it at 15, so we may not reject some candidates. But we never incorrectly
    /// reject a true divisor.
    #[inline]
    pub fn buckets_leq(&self, other: &Self) -> bool {
        match (&self.buckets, &other.buckets) {
            (Some(b1), Some(b2)) => compare_buckets(b1, b2),
            _ => true, // If buckets not computed, conservatively pass
        }
    }
}

impl PresenceMask {
    /// Create a PresenceMask from an exponent vector.
    fn from_exponents(exp: &[usize]) -> Self {
        let nvars = exp.len();

        if nvars <= 64 {
            let mut mask = 0u64;
            for (i, &e) in exp.iter().enumerate() {
                if e > 0 {
                    mask |= 1u64 << i;
                }
            }
            PresenceMask::U64(mask)
        } else {
            let nchunks = (nvars + 63) / 64;
            let mut chunks = vec![0u64; nchunks];
            for (i, &e) in exp.iter().enumerate() {
                if e > 0 {
                    let chunk_idx = i / 64;
                    let bit_idx = i % 64;
                    chunks[chunk_idx] |= 1u64 << bit_idx;
                }
            }
            PresenceMask::Chunks(chunks)
        }
    }

    /// Check if this presence mask is a subset of another.
    /// Returns true if every bit set in self is also set in other.
    #[inline]
    fn is_subset(&self, other: &Self) -> bool {
        match (self, other) {
            (PresenceMask::U64(a), PresenceMask::U64(b)) => {
                // a is subset of b iff (a & !b) == 0
                (a & !b) == 0
            }
            (PresenceMask::Chunks(a), PresenceMask::Chunks(b)) => {
                // Check each chunk
                for (chunk_a, chunk_b) in a.iter().zip(b.iter()) {
                    if (chunk_a & !chunk_b) != 0 {
                        return false;
                    }
                }
                true
            }
            _ => {
                // Mismatched types (shouldn't happen in practice)
                // Conservatively return true to avoid rejecting valid candidates
                true
            }
        }
    }
}

/// Pack exponents into nibble buckets.
/// Each u64 holds 16 variables, with 4 bits per variable storing min(exp, 15).
fn pack_nibbles(exp: &[usize]) -> Vec<u64> {
    let nvars = exp.len();
    let nwords = (nvars + 15) / 16;
    let mut buckets = vec![0u64; nwords];

    for (i, &e) in exp.iter().enumerate() {
        let word_idx = i / 16;
        let nibble_idx = i % 16;
        let capped = std::cmp::min(e, 15) as u64;
        buckets[word_idx] |= capped << (nibble_idx * 4);
    }

    buckets
}

/// Compare two bucket vectors.
/// Returns true if for every variable i, bucket_a[i] <= bucket_b[i].
#[inline]
fn compare_buckets(buckets_a: &[u64], buckets_b: &[u64]) -> bool {
    for (word_a, word_b) in buckets_a.iter().zip(buckets_b.iter()) {
        // Extract and compare all 16 nibbles in this word
        // For efficiency, we can compare the whole word at once in some cases
        if word_a == word_b {
            continue; // All nibbles equal, move to next word
        }

        // Compare nibble by nibble
        for shift in 0..16 {
            let nibble_a = (word_a >> (shift * 4)) & 0xF;
            let nibble_b = (word_b >> (shift * 4)) & 0xF;
            if nibble_a > nibble_b {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_presence_mask_u64() {
        let exp1 = vec![1, 0, 2, 0]; // x^1 * z^2
        let exp2 = vec![2, 1, 3, 0]; // x^2 * y * z^3
        let exp3 = vec![0, 1, 0, 0]; // y

        let sig1 = MonSig::from_exponents(&exp1);
        let sig2 = MonSig::from_exponents(&exp2);
        let sig3 = MonSig::from_exponents(&exp3);

        // sig1 has vars {0, 2}, sig2 has vars {0, 1, 2}, sig3 has vars {1}
        assert!(sig1.presence_subset(&sig2)); // {0,2} ⊆ {0,1,2}
        assert!(sig3.presence_subset(&sig2)); // {1} ⊆ {0,1,2}
        assert!(!sig2.presence_subset(&sig1)); // {0,1,2} ⊄ {0,2}
        assert!(!sig2.presence_subset(&sig3)); // {0,1,2} ⊄ {1}
    }

    #[test]
    fn test_degree_check() {
        let exp1 = vec![1, 2]; // degree 3
        let exp2 = vec![2, 3]; // degree 5

        let sig1 = MonSig::from_exponents(&exp1);
        let sig2 = MonSig::from_exponents(&exp2);

        assert_eq!(sig1.deg, 3);
        assert_eq!(sig2.deg, 5);
    }

    #[test]
    fn test_buckets_leq() {
        let exp1 = vec![1, 2, 3];    // all < 15
        let exp2 = vec![2, 3, 4];    // all < 15, all >= exp1
        let exp3 = vec![0, 3, 2];    // mixed
        let exp4 = vec![20, 30, 10]; // some > 15 (capped)

        let sig1 = MonSig::from_exponents(&exp1);
        let sig2 = MonSig::from_exponents(&exp2);
        let sig3 = MonSig::from_exponents(&exp3);
        let sig4 = MonSig::from_exponents(&exp4);

        assert!(sig1.buckets_leq(&sig2)); // [1,2,3] <= [2,3,4]
        assert!(!sig2.buckets_leq(&sig1)); // [2,3,4] > [1,2,3]
        assert!(!sig1.buckets_leq(&sig3)); // [1,2,3] vs [0,3,2]: 1>0 fails

        // Test with capped values
        assert!(sig1.buckets_leq(&sig4)); // [1,2,3] <= [15,15,10]
    }

    #[test]
    fn test_presence_chunks() {
        // Test with > 64 variables
        let mut exp1 = vec![0; 100];
        exp1[0] = 1;
        exp1[65] = 1;

        let mut exp2 = vec![0; 100];
        exp2[0] = 1;
        exp2[32] = 1;
        exp2[65] = 1;

        let sig1 = MonSig::from_exponents(&exp1);
        let sig2 = MonSig::from_exponents(&exp2);

        assert!(sig1.presence_subset(&sig2)); // {0,65} ⊆ {0,32,65}
        assert!(!sig2.presence_subset(&sig1)); // {0,32,65} ⊄ {0,65}
    }

    #[test]
    fn test_divisibility_prefilter() {
        // Simulate a divisibility check: does x^2*y divide x^3*y^2*z?
        let lm_exp = vec![2, 1, 0];    // x^2*y, degree 3
        let tgt_exp = vec![3, 2, 1];   // x^3*y^2*z, degree 6

        let lm_sig = MonSig::from_exponents(&lm_exp);
        let tgt_sig = MonSig::from_exponents(&tgt_exp);

        // Should pass all prefilter checks
        assert!(lm_sig.deg <= tgt_sig.deg);
        assert!(lm_sig.presence_subset(&tgt_sig));
        assert!(lm_sig.buckets_leq(&tgt_sig));
    }

    #[test]
    fn test_divisibility_prefilter_reject() {
        // x^2*y*z should NOT divide x*y (degree too high)
        let lm_exp = vec![2, 1, 1];    // x^2*y*z, degree 4
        let tgt_exp = vec![1, 1, 0];   // x*y, degree 2

        let lm_sig = MonSig::from_exponents(&lm_exp);
        let tgt_sig = MonSig::from_exponents(&tgt_exp);

        // Should fail degree check
        assert!(lm_sig.deg > tgt_sig.deg);

        // Also fails presence check (z present in lm but not target)
        assert!(!lm_sig.presence_subset(&tgt_sig));
    }
}
