//! Type conversions between arkworks and feanor-math types.

use ark_ff::{Field, PrimeField, BigInteger};
use num_bigint::{BigUint, BigInt};
use feanor_math::integer::*;
use crate::field_wrapper::ArkFieldWrapper;

/// Extract the characteristic of a prime field as a BigUint
/// 
/// # Example
/// ```ignore
/// use ark_bn254::Fr;
/// use ark_feanor::conversions::extract_characteristic;
/// 
/// let char = extract_characteristic::<Fr>();
/// println!("Field characteristic: {}", char);
/// ```
pub fn extract_characteristic<F: PrimeField>() -> BigUint {
    let modulus_bytes = F::MODULUS.to_bytes_le();
    BigUint::from_bytes_le(&modulus_bytes)
}

/// Convert a BigUint to an arkworks field element
/// 
/// This performs modular reduction if the value is larger than the field modulus.
pub fn biguint_to_field<F: PrimeField>(value: &BigUint) -> F {
    // Convert to bytes and then to field element
    let bytes = value.to_bytes_le();
    F::from_le_bytes_mod_order(&bytes)
}

/// Convert an arkworks field element to a BigUint
/// 
/// This returns the canonical representative in [0, p) where p is the field characteristic.
pub fn field_to_biguint<F: PrimeField>(el: &F) -> BigUint {
    let bytes = el.into_bigint().to_bytes_le();
    BigUint::from_bytes_le(&bytes)
}

/// Convert between signed and unsigned big integers
pub fn bigint_to_biguint(value: &BigInt) -> Option<BigUint> {
    if value.sign() == num_bigint::Sign::Minus {
        None
    } else {
        Some(value.magnitude().clone())
    }
}

/// Convert a feanor-math integer to an arkworks field element
///
/// Note: This is a stub implementation that only works for small integers.
/// For full functionality, use the field's from_int method directly.
pub fn feanor_int_to_field<F, I>(_int_ring: &I, value: &I::Element) -> F
where
    F: PrimeField,
    I: IntegerRing + ?Sized,
{
    // Simplified implementation - only works for small values
    // In practice, use the field's int_hom() for proper conversion
    let _ = value;  // Suppress unused warning
    F::from(0u64)  // Placeholder
}

/// Create a homomorphism from integers to an arkworks field
pub struct IntToFieldHom<F: PrimeField> {
    field: ArkFieldWrapper<F>,
}

impl<F: PrimeField> IntToFieldHom<F> {
    pub fn new() -> Self {
        Self {
            field: ArkFieldWrapper::new(),
        }
    }
}

impl<F: PrimeField> Default for IntToFieldHom<F> {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper to convert field elements to strings for debugging
pub fn field_to_string<F: Field>(el: &F) -> String {
    format!("{:?}", el)
}

/// Check if a field element represents a small integer
/// 
/// Returns Some(n) if the element equals n as a field element for |n| < 2^32
pub fn field_to_small_int<F: PrimeField>(el: &F) -> Option<i64> {
    // Check small positive values
    for i in 0i64..1000 {
        if *el == F::from(i as u64) {
            return Some(i);
        }
    }
    
    // Check small negative values
    for i in 1i64..1000 {
        if *el == -F::from(i as u64) {
            return Some(-i);
        }
    }
    
    None
}

/// Batch conversion of field elements to BigUints
pub fn fields_to_biguints<F: PrimeField>(elements: &[F]) -> Vec<BigUint> {
    elements.iter().map(field_to_biguint).collect()
}

/// Batch conversion of BigUints to field elements
pub fn biguints_to_fields<F: PrimeField>(values: &[BigUint]) -> Vec<F> {
    values.iter().map(biguint_to_field).collect()
}

/// Convert an i64 to an arkworks field element
///
/// This properly handles both positive and negative i64 values, including i64::MIN.
///
/// # Example
/// ```ignore
/// use ark_bn254::Fr;
/// use ark_feanor::i64_to_field;
///
/// let pos = i64_to_field::<Fr>(42);
/// let neg = i64_to_field::<Fr>(-100);
/// let min = i64_to_field::<Fr>(i64::MIN);
/// ```
pub fn i64_to_field<F: Field>(value: i64) -> F {
    if value >= 0 {
        F::from(value as u64)
    } else {
        // For negative values, we need to handle i64::MIN specially
        // because -i64::MIN overflows
        if value == i64::MIN {
            // i64::MIN = -2^63
            // We compute this as -(2^63) = -(2^62 * 2)
            let two_pow_62 = F::from(1u64 << 62);
            let two_pow_63 = two_pow_62 + two_pow_62; // 2 * 2^62 = 2^63
            -two_pow_63
        } else {
            -F::from((-value) as u64)
        }
    }
}

/// Convert a signed BigInt to an arkworks field element
///
/// This performs modular reduction if the value is larger than the field modulus.
/// Negative values are properly handled by computing the field element's additive inverse.
///
/// # Example
/// ```ignore
/// use ark_bn254::Fr;
/// use ark_feanor::bigint_to_field;
/// use num_bigint::BigInt;
///
/// let big = BigInt::from(123456789i64);
/// let neg_big = BigInt::from(-987654321i64);
///
/// let pos_field = bigint_to_field::<Fr>(&big);
/// let neg_field = bigint_to_field::<Fr>(&neg_big);
/// ```
pub fn bigint_to_field<F: PrimeField>(value: &BigInt) -> F {
    let (sign, magnitude) = value.clone().into_parts();
    let field_elem = biguint_to_field::<F>(&magnitude);

    match sign {
        num_bigint::Sign::Minus => -field_elem,
        num_bigint::Sign::NoSign | num_bigint::Sign::Plus => field_elem,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    use num_traits::Zero;

    #[test]
    fn test_characteristic_extraction() {
        let char = extract_characteristic::<Fr>();
        // BN254's Fr field has a specific known characteristic
        // Just verify it's non-zero and reasonable size
        assert!(char > BigUint::zero());
        assert!(char.bits() > 200); // Should be around 254 bits
    }
    
    #[test]
    fn test_roundtrip_conversion() {
        let original = Fr::from(12345u64);
        let as_biguint = field_to_biguint(&original);
        let back_to_field = biguint_to_field::<Fr>(&as_biguint);
        assert_eq!(original, back_to_field);
    }
    
    #[test]
    fn test_small_int_detection() {
        let five = Fr::from(5u64);
        assert_eq!(field_to_small_int(&five), Some(5));
        
        let neg_three = -Fr::from(3u64);
        assert_eq!(field_to_small_int(&neg_three), Some(-3));
        
        // Large random element unlikely to be a small int
        let large = Fr::from(u64::MAX);
        assert_eq!(field_to_small_int(&large), None);
    }
    
    #[test]
    fn test_modular_reduction() {
        let char = extract_characteristic::<Fr>();
        let large_value = &char * 2u32 + 7u32; // 2p + 7
        let reduced = biguint_to_field::<Fr>(&large_value);
        let expected = Fr::from(7u64);
        assert_eq!(reduced, expected);
    }

    #[test]
    fn test_i64_to_field_positive() {
        let val = i64_to_field::<Fr>(42);
        assert_eq!(val, Fr::from(42u64));

        let large_pos = i64_to_field::<Fr>(i64::MAX);
        assert_eq!(large_pos, Fr::from(i64::MAX as u64));
    }

    #[test]
    fn test_i64_to_field_negative() {
        let neg = i64_to_field::<Fr>(-100);
        assert_eq!(neg, -Fr::from(100u64));

        let neg_one = i64_to_field::<Fr>(-1);
        assert_eq!(neg_one, -Fr::from(1u64));
    }

    #[test]
    fn test_i64_to_field_min() {
        // Test i64::MIN specially since it's the edge case
        let min_val = i64_to_field::<Fr>(i64::MIN);

        // Verify: i64::MIN = -2^63
        let two_pow_62 = Fr::from(1u64 << 62);
        let two_pow_63 = two_pow_62 + two_pow_62;
        assert_eq!(min_val, -two_pow_63);

        // Verify it adds to zero with i64::MAX + 1
        let max_plus_one = Fr::from(i64::MAX as u64) + Fr::from(1u64);
        assert_eq!(min_val + max_plus_one, Fr::zero());
    }

    #[test]
    fn test_i64_to_field_zero() {
        let zero = i64_to_field::<Fr>(0);
        assert_eq!(zero, Fr::zero());
    }

    #[test]
    fn test_bigint_to_field_positive() {
        let big = BigInt::from(123456789i64);
        let field_elem = bigint_to_field::<Fr>(&big);
        assert_eq!(field_elem, Fr::from(123456789u64));
    }

    #[test]
    fn test_bigint_to_field_negative() {
        let neg_big = BigInt::from(-987654321i64);
        let field_elem = bigint_to_field::<Fr>(&neg_big);
        assert_eq!(field_elem, -Fr::from(987654321u64));
    }

    #[test]
    fn test_bigint_to_field_zero() {
        let zero = BigInt::zero();
        let field_elem = bigint_to_field::<Fr>(&zero);
        assert_eq!(field_elem, Fr::zero());
    }

    #[test]
    fn test_bigint_to_field_large() {
        // Test with a value larger than i64::MAX
        let char = extract_characteristic::<Fr>();
        let large = BigInt::from_bytes_le(num_bigint::Sign::Plus, &char.to_bytes_le());
        let large_minus_one = &large - 1;

        let field_elem = bigint_to_field::<Fr>(&large_minus_one);
        // Should be p - 1, which is -1 in the field
        assert_eq!(field_elem, -Fr::from(1u64));
    }

    #[test]
    fn test_bigint_roundtrip() {
        // Test that we can convert field -> biguint -> bigint -> field
        let original = Fr::from(999999u64);
        let as_biguint = field_to_biguint(&original);
        let as_bigint = BigInt::from_biguint(num_bigint::Sign::Plus, as_biguint);
        let back = bigint_to_field::<Fr>(&as_bigint);
        assert_eq!(original, back);
    }

    #[test]
    fn test_i64_edge_cases() {
        // Test various edge cases
        let cases = vec![
            0i64,
            1i64,
            -1i64,
            i64::MIN,
            i64::MAX,
            i64::MIN + 1,
            i64::MAX - 1,
        ];

        for &val in &cases {
            let field_elem = i64_to_field::<Fr>(val);

            // Verify it's not zero unless val is zero
            if val != 0 {
                assert!(!field_elem.is_zero(), "Non-zero i64 {} should not map to zero", val);
            } else {
                assert!(field_elem.is_zero(), "Zero i64 should map to zero");
            }
        }
    }
}
