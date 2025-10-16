//! Prime field specializations that enable division and field operations.

use ark_ff::{Field, PrimeField, BigInteger};
use feanor_math::divisibility::*;
use feanor_math::pid::*;
use crate::field_wrapper::ArkFieldWrapper;

/// Marker trait to indicate that an arkworks field implements PrimeField
pub trait IsPrimeField: Field {
    fn field_characteristic() -> num_bigint::BigUint;
}

impl<F: PrimeField> IsPrimeField for F {
    fn field_characteristic() -> num_bigint::BigUint {
        // Convert arkworks BigInt to num_bigint::BigUint
        let modulus_bytes = F::MODULUS.to_bytes_le();
        num_bigint::BigUint::from_bytes_le(&modulus_bytes)
    }
}

// Implement DivisibilityRing for prime fields (allows division by non-zero elements)
impl<F: PrimeField> DivisibilityRing for ArkFieldWrapper<F> {
    fn checked_left_div(&self, lhs: &Self::Element, rhs: &Self::Element) -> Option<Self::Element> {
        if rhs.is_zero() {
            None
        } else {
            rhs.inverse().map(|inv| *lhs * inv)
        }
    }

    fn is_unit(&self, el: &Self::Element) -> bool {
        !el.is_zero()
    }
}

// Implement Domain (no zero divisors)
impl<F: PrimeField> Domain for ArkFieldWrapper<F> {}

// Implement PrincipalIdealRing (every ideal is principal)
impl<F: PrimeField> PrincipalIdealRing for ArkFieldWrapper<F> {
    fn ideal_gen(&self, lhs: &Self::Element, rhs: &Self::Element) -> Self::Element {
        // In a field, the GCD is either 0 (if both are zero) or 1
        if lhs.is_zero() && rhs.is_zero() {
            F::zero()
        } else {
            // At least one is non-zero, so GCD is 1
            F::one()
        }
    }

    fn extended_ideal_gen(&self, lhs: &Self::Element, rhs: &Self::Element)
        -> (Self::Element, Self::Element, Self::Element) {
        // In a field, the GCD is either 0 (if both are zero) or 1
        if lhs.is_zero() && rhs.is_zero() {
            (F::zero(), F::one(), F::zero())
        } else if lhs.is_zero() {
            (*rhs, F::zero(), F::one())
        } else if rhs.is_zero() {
            (*lhs, F::one(), F::zero())
        } else {
            // Both non-zero, GCD is 1
            // We need to find a, b such that a*lhs + b*rhs = 1
            // Since we're in a field: 1 = lhs * (1/lhs) + rhs * 0
            let lhs_inv = lhs.inverse().unwrap();
            (F::one(), lhs_inv, F::zero())
        }
    }

    fn checked_div_min(&self, lhs: &Self::Element, rhs: &Self::Element) -> Option<Self::Element> {
        self.checked_left_div(lhs, rhs)
    }
}

// Implement EuclideanRing
impl<F: PrimeField> EuclideanRing for ArkFieldWrapper<F> {
    fn euclidean_div_rem(&self, lhs: Self::Element, rhs: &Self::Element)
        -> (Self::Element, Self::Element) {
        if rhs.is_zero() {
            panic!("Division by zero");
        }
        let quotient = lhs * rhs.inverse().unwrap();
        (quotient, F::zero())
    }

    fn euclidean_deg(&self, el: &Self::Element) -> Option<usize> {
        if el.is_zero() {
            None
        } else {
            Some(0)
        }
    }
}

// Finally, implement Field trait from feanor-math
impl<F: PrimeField> feanor_math::field::Field for ArkFieldWrapper<F> {
    fn div(&self, lhs: &Self::Element, rhs: &Self::Element) -> Self::Element {
        if rhs.is_zero() {
            panic!("Division by zero");
        }
        *lhs * rhs.inverse().expect("Division by zero")
    }
}

// Helper trait to extract field properties
pub trait FieldProperties {
    /// Get the characteristic of the field as a BigUint
    fn field_characteristic_biguint(&self) -> num_bigint::BigUint;

    /// Check if the field is a prime field
    fn is_prime_field(&self) -> bool;

    /// Get the field size (number of elements)
    fn field_size(&self) -> num_bigint::BigUint;
}

impl<F: PrimeField> FieldProperties for ArkFieldWrapper<F> {
    fn field_characteristic_biguint(&self) -> num_bigint::BigUint {
        F::field_characteristic()
    }

    fn is_prime_field(&self) -> bool {
        true
    }

    fn field_size(&self) -> num_bigint::BigUint {
        // For prime fields, size = characteristic
        self.field_characteristic_biguint()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;

    #[test]
    fn test_division() {
        let field = ArkFieldWrapper::<Fr>::new();
        
        let a = field.from_int(10);
        let b = field.from_int(5);
        let c = field.div(&a, &b);
        let expected = field.from_int(2);
        assert!(field.eq_el(&c, &expected));
    }

    #[test]
    fn test_inverse() {
        let field = ArkFieldWrapper::<Fr>::new();
        
        let a = field.from_int(7);
        let a_inv = field.checked_left_div(&field.one(), &a).unwrap();
        let product = field.mul_ref(&a, &a_inv);
        assert!(field.is_one(&product));
    }

    #[test]
    fn test_is_unit() {
        let field = ArkFieldWrapper::<Fr>::new();
        
        assert!(!field.is_unit(&field.zero()));
        assert!(field.is_unit(&field.one()));
        assert!(field.is_unit(&field.from_int(42)));
    }

    #[test]
    #[should_panic(expected = "Division by zero")]
    fn test_division_by_zero() {
        let field = ArkFieldWrapper::<Fr>::new();
        let a = field.from_int(5);
        let _ = field.div(&a, &field.zero());
    }
}
