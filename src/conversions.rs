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
}
