//! Tests for field properties and polynomial system basics over cryptographic fields
//!
//! Note: Gröbner basis computation in feanor-math v3 requires unstable APIs.
//! These tests focus on the field wrapper functionality and basic polynomial operations.

use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use ark_bls12_381::Fr as BlsFr;
use ark_feanor::prime_field::FieldProperties;
use feanor_math::homomorphism::Homomorphism;

#[test]
fn test_field_properties_bn254() {
    let field = ArkFieldWrapper::<BnFr>::new();

    // Test that the field is indeed a prime field
    assert!(field.is_prime_field());

    // Get characteristic
    let char = field.field_characteristic_biguint();
    assert!(char > num_bigint::BigUint::from(0u32));

    // Field size should equal characteristic for prime fields
    let size = field.field_size();
    assert_eq!(size, char);

    println!("BN254 Fr characteristic: {} bits", char.bits());
}

#[test]
fn test_field_properties_bls12_381() {
    let field = ArkFieldWrapper::<BlsFr>::new();

    // Test that the field is indeed a prime field
    assert!(field.is_prime_field());

    // Get characteristic
    let char = field.field_characteristic_biguint();
    assert!(char > num_bigint::BigUint::from(0u32));

    // Verify field properties
    assert!(field.is_prime_field());
    println!("BLS12-381 Fr characteristic: {} bits", char.bits());
}

#[test]
fn test_different_field_characteristics() {
    let bn_field = ArkFieldWrapper::<BnFr>::new();
    let bls_field = ArkFieldWrapper::<BlsFr>::new();

    let bn_char = bn_field.field_characteristic_biguint();
    let bls_char = bls_field.field_characteristic_biguint();

    // Different curves should have different characteristics
    assert_ne!(bn_char, bls_char);

    println!("BN254 Fr: {} bits", bn_char.bits());
    println!("BLS12-381 Fr: {} bits", bls_char.bits());
}

#[test]
fn test_field_element_properties() {
    use feanor_math::field::Field;

    let field = ArkFieldWrapper::<BnFr>::new();

    // Every non-zero element should have an inverse
    for i in 1..10 {
        let el = field.from_int(i);
        let inv = field.div(&field.one(), &el);
        let product = field.mul_ref(&el, &inv);
        assert!(field.is_one(&product), "Failed for i={}", i);
    }
}

#[test]
fn test_pre_configured_fields() {
    // Test that pre-configured fields work
    let bn_field = &*BN254_FR;
    let bls_field = &*BLS12_381_FR;

    // Basic arithmetic
    let a = bn_field.int_hom().map(100);
    let b = bn_field.int_hom().map(50);
    let sum = bn_field.add_ref(&a, &b);
    let expected = bn_field.int_hom().map(150);
    assert!(bn_field.eq_el(&sum, &expected));

    // Division with BLS field
    let x = bls_field.int_hom().map(100);
    let y = bls_field.int_hom().map(50);
    let quotient = bls_field.checked_left_div(&x, &y).expect("Division failed");
    let reconstructed = bls_field.mul_ref(&quotient, &y);
    assert!(bls_field.eq_el(&reconstructed, &x));
}

#[test]
fn test_field_arithmetic_consistency() {
    let field = ArkFieldWrapper::<BlsFr>::new();

    // Test that (a + b) * c = a*c + b*c
    let a = field.from_int(7);
    let b = field.from_int(13);
    let c = field.from_int(19);

    let left = field.mul_ref(&field.add_ref(&a, &b), &c);
    let right = field.add_ref(
        &field.mul_ref(&a, &c),
        &field.mul_ref(&b, &c),
    );

    assert!(field.eq_el(&left, &right));
}

#[test]
fn test_characteristic_field_operations() {
    use feanor_math::divisibility::DivisibilityRing;

    let field = ArkFieldWrapper::<BnFr>::new();

    // Test division
    let a = field.from_int(42);
    let b = field.from_int(7);

    let quotient = field.checked_left_div(&a, &b).expect("Division should succeed");
    let expected = field.from_int(6);
    assert!(field.eq_el(&quotient, &expected));

    // Test that zero is not a unit
    assert!(!field.is_unit(&field.zero()));

    // Test that all non-zero elements are units
    for i in 1..20 {
        let el = field.from_int(i);
        assert!(field.is_unit(&el), "Element {} should be a unit", i);
    }
}
