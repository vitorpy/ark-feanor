//! Common field type aliases for easy access to popular cryptographic fields.
//! 
//! This module provides pre-configured wrappers for commonly used fields from
//! various elliptic curves, making it easy to work with them in feanor-math.

use ark_bn254::{Fr as BnFr, Fq as BnFq};
use ark_bls12_381::{Fr as BlsFr, Fq as BlsFq};
use feanor_math::ring::*;
use crate::field_wrapper::ArkFieldWrapper;

// Type aliases for BN254 curve fields
/// The scalar field of BN254 (Fr)
pub type BN254ScalarField = ArkFieldWrapper<BnFr>;
/// The base field of BN254 (Fq)
pub type BN254BaseField = ArkFieldWrapper<BnFq>;

// Type aliases for BLS12-381 curve fields
/// The scalar field of BLS12-381 (Fr)
pub type BLS12_381ScalarField = ArkFieldWrapper<BlsFr>;
/// The base field of BLS12-381 (Fq)
pub type BLS12_381BaseField = ArkFieldWrapper<BlsFq>;

// RingValue wrappers for convenient usage
/// BN254 Fr field as a RingValue
pub type BnFrRing = RingValue<ArkFieldWrapper<BnFr>>;
/// BN254 Fq field as a RingValue
pub type BnFqRing = RingValue<ArkFieldWrapper<BnFq>>;
/// BLS12-381 Fr field as a RingValue
pub type BlsFrRing = RingValue<ArkFieldWrapper<BlsFr>>;
/// BLS12-381 Fq field as a RingValue
pub type BlsFqRing = RingValue<ArkFieldWrapper<BlsFq>>;

// Const instances for direct use
lazy_static::lazy_static! {
    /// Pre-initialized BN254 scalar field (Fr)
    pub static ref BN254_FR: BnFrRing = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    
    /// Pre-initialized BN254 base field (Fq)
    pub static ref BN254_FQ: BnFqRing = RingValue::from(ArkFieldWrapper::<BnFq>::new());
    
    /// Pre-initialized BLS12-381 scalar field (Fr)
    pub static ref BLS12_381_FR: BlsFrRing = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    
    /// Pre-initialized BLS12-381 base field (Fq)
    pub static ref BLS12_381_FQ: BlsFqRing = RingValue::from(ArkFieldWrapper::<BlsFq>::new());
}

/// Create a new BN254 scalar field instance
#[inline]
pub fn bn254_scalar_field() -> BN254ScalarField {
    ArkFieldWrapper::new()
}

/// Create a new BN254 base field instance
#[inline]
pub fn bn254_base_field() -> BN254BaseField {
    ArkFieldWrapper::new()
}

/// Create a new BLS12-381 scalar field instance
#[inline]
pub fn bls12_381_scalar_field() -> BLS12_381ScalarField {
    ArkFieldWrapper::new()
}

/// Create a new BLS12-381 base field instance
#[inline]
pub fn bls12_381_base_field() -> BLS12_381BaseField {
    ArkFieldWrapper::new()
}

/// Helper trait to identify which curve a field belongs to
pub trait CurveField {
    /// The name of the curve this field belongs to
    fn curve_name() -> &'static str;
    
    /// Whether this is the scalar field (Fr) or base field (Fq)
    fn field_type() -> FieldType;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldType {
    ScalarField,
    BaseField,
}

impl CurveField for BnFr {
    fn curve_name() -> &'static str {
        "BN254"
    }
    
    fn field_type() -> FieldType {
        FieldType::ScalarField
    }
}

impl CurveField for BnFq {
    fn curve_name() -> &'static str {
        "BN254"
    }
    
    fn field_type() -> FieldType {
        FieldType::BaseField
    }
}

impl CurveField for BlsFr {
    fn curve_name() -> &'static str {
        "BLS12-381"
    }
    
    fn field_type() -> FieldType {
        FieldType::ScalarField
    }
}

impl CurveField for BlsFq {
    fn curve_name() -> &'static str {
        "BLS12-381"
    }
    
    fn field_type() -> FieldType {
        FieldType::BaseField
    }
}

/// Get field information as a string
pub fn field_info<F: CurveField>() -> String {
    format!(
        "{} {:?}",
        F::curve_name(),
        F::field_type()
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_field_creation() {
        let _bn_fr = bn254_scalar_field();
        let _bn_fq = bn254_base_field();
        let _bls_fr = bls12_381_scalar_field();
        let _bls_fq = bls12_381_base_field();
    }
    
    #[test]
    fn test_field_info() {
        assert_eq!(field_info::<BnFr>(), "BN254 ScalarField");
        assert_eq!(field_info::<BlsFq>(), "BLS12-381 BaseField");
    }
    
    #[test]
    fn test_static_fields() {
        // Access the lazy static fields
        let field = &*BN254_FR;
        let one = field.one();
        let two = field.int_hom().map(2);
        let three = field.add(one, two);
        assert_eq!(field.int_hom().map(3), three);
    }
}
