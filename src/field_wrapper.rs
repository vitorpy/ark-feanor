//! Core field wrapper that implements feanor-math's RingBase trait for arkworks fields.

use ark_ff::Field;
use feanor_math::ring::*;
use feanor_math::integer::IntegerRing;
use std::marker::PhantomData;
use std::fmt::{self, Debug, Display};

/// Wrapper struct that allows arkworks `Field` types to be used as feanor-math rings.
/// 
/// This is a zero-sized type that acts as a bridge between the two type systems.
/// All arkworks field operations are mapped to their feanor-math equivalents.
#[derive(Clone, Copy)]
pub struct ArkFieldWrapper<F: Field> {
    _phantom: PhantomData<F>,
}

impl<F: Field> ArkFieldWrapper<F> {
    /// Create a new instance of the field wrapper.
    /// 
    /// # Example
    /// ```ignore
    /// use ark_bn254::Fr;
    /// use ark_feanor::ArkFieldWrapper;
    /// 
    /// let field = ArkFieldWrapper::<Fr>::new();
    /// ```
    pub const fn new() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<F: Field> Default for ArkFieldWrapper<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Field> Debug for ArkFieldWrapper<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ArkFieldWrapper<{}>", std::any::type_name::<F>())
    }
}

impl<F: Field> Display for ArkFieldWrapper<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Field({})", std::any::type_name::<F>())
    }
}

// All instances of the same field type are considered equal
impl<F: Field> PartialEq for ArkFieldWrapper<F> {
    fn eq(&self, _other: &Self) -> bool {
        true // All instances of the same field type are the same ring
    }
}

impl<F: Field> Eq for ArkFieldWrapper<F> {}

impl<F: Field> RingBase for ArkFieldWrapper<F> {
    type Element = F;

    fn clone_el(&self, el: &Self::Element) -> Self::Element {
        *el
    }

    fn add_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        *lhs += rhs;
    }

    fn add_ref(&self, lhs: &Self::Element, rhs: &Self::Element) -> Self::Element {
        *lhs + *rhs
    }

    fn negate_inplace(&self, el: &mut Self::Element) {
        *el = -*el;
    }

    fn mul_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        *lhs *= rhs;
    }

    fn mul_ref(&self, lhs: &Self::Element, rhs: &Self::Element) -> Self::Element {
        *lhs * *rhs
    }

    fn from_int(&self, value: i32) -> Self::Element {
        if value >= 0 {
            F::from(value as u64)
        } else {
            -F::from((-value) as u64)
        }
    }

    fn eq_el(&self, lhs: &Self::Element, rhs: &Self::Element) -> bool {
        lhs == rhs
    }

    fn is_zero(&self, el: &Self::Element) -> bool {
        el.is_zero()
    }

    fn is_one(&self, el: &Self::Element) -> bool {
        el.is_one()
    }

    fn is_neg_one(&self, el: &Self::Element) -> bool {
        *el == -F::one()
    }

    fn zero(&self) -> Self::Element {
        F::zero()
    }

    fn one(&self) -> Self::Element {
        F::one()
    }

    fn neg_one(&self) -> Self::Element {
        -F::one()
    }

    fn is_commutative(&self) -> bool {
        true // All fields are commutative
    }

    fn is_noetherian(&self) -> bool {
        true // All fields are noetherian
    }

    fn is_approximate(&self) -> bool {
        false // Finite fields are exact
    }

    fn dbg<'a>(&self, el: &Self::Element, out: &mut std::fmt::Formatter<'a>) -> std::fmt::Result {
        write!(out, "{:?}", el)
    }

    fn dbg_within<'a>(&self, el: &Self::Element, out: &mut std::fmt::Formatter<'a>, _env: EnvBindingStrength) -> std::fmt::Result {
        write!(out, "{:?}", el)
    }

    fn characteristic<I>(&self, _int_ring: I) -> Option<<I::Type as RingBase>::Element>
    where
        I: RingStore + Copy,
        I::Type: IntegerRing,
    {
        // For fields with known characteristic, we would return it here
        // For now, return None to indicate unknown/zero characteristic
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;

    #[test]
    fn test_basic_operations() {
        let field = ArkFieldWrapper::<Fr>::new();
        
        // Test zero and one
        assert!(field.is_zero(&field.zero()));
        assert!(field.is_one(&field.one()));
        
        // Test basic arithmetic
        let a = field.from_int(5);
        let b = field.from_int(3);
        let c = field.add_ref(&a, &b);
        let expected = field.from_int(8);
        assert!(field.eq_el(&c, &expected));
        
        // Test multiplication
        let d = field.mul_ref(&a, &b);
        let expected = field.from_int(15);
        assert!(field.eq_el(&d, &expected));
        
        // Test negation
        let neg_a = field.negate(field.clone_el(&a));
        let zero = field.add_ref(&a, &neg_a);
        assert!(field.is_zero(&zero));
    }
}
