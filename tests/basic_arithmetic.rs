//! Tests for basic arithmetic operations and ring axioms

use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use ark_bls12_381::{Fr as BlsFr, Fq as BlsFq};
use feanor_math::field::Field;
use ark_ff::{Field as ArkField, One, Zero};

#[test]
fn test_ring_axioms_bn254() {
    let field = ArkFieldWrapper::<BnFr>::new();
    
    // Test additive identity
    let a = field.from_int(42);
    let zero = field.zero();
    assert!(field.eq_el(&field.add_ref(&a, &zero), &a));
    assert!(field.eq_el(&field.add_ref(&zero, &a), &a));
    
    // Test multiplicative identity
    let one = field.one();
    assert!(field.eq_el(&field.mul_ref(&a, &one), &a));
    assert!(field.eq_el(&field.mul_ref(&one, &a), &a));
    
    // Test additive inverse
    let neg_a = field.negate(field.clone_el(&a));
    let sum = field.add_ref(&a, &neg_a);
    assert!(field.is_zero(&sum));
    
    // Test commutativity
    let b = field.from_int(17);
    assert!(field.eq_el(
        &field.add_ref(&a, &b),
        &field.add_ref(&b, &a)
    ));
    assert!(field.eq_el(
        &field.mul_ref(&a, &b),
        &field.mul_ref(&b, &a)
    ));
    
    // Test associativity
    let c = field.from_int(23);
    let left = field.add_ref(&field.add_ref(&a, &b), &c);
    let right = field.add_ref(&a, &field.add_ref(&b, &c));
    assert!(field.eq_el(&left, &right));
    
    let left = field.mul_ref(&field.mul_ref(&a, &b), &c);
    let right = field.mul_ref(&a, &field.mul_ref(&b, &c));
    assert!(field.eq_el(&left, &right));
    
    // Test distributivity
    let left = field.mul_ref(&a, &field.add_ref(&b, &c));
    let right = field.add_ref(
        &field.mul_ref(&a, &b),
        &field.mul_ref(&a, &c)
    );
    assert!(field.eq_el(&left, &right));
}

#[test]
fn test_division_bls12_381() {
    let field = ArkFieldWrapper::<BlsFr>::new();
    
    // Test simple division
    let a = field.from_int(20);
    let b = field.from_int(4);
    let c = field.div(&a, &b);
    let expected = field.from_int(5);
    assert!(field.eq_el(&c, &expected));
    
    // Test inverse
    let x = field.from_int(7);
    let x_inv = field.div(&field.one(), &x);
    let product = field.mul_ref(&x, &x_inv);
    assert!(field.is_one(&product));
    
    // Test division cancellation: (a * b) / b = a
    let a = field.from_int(13);
    let b = field.from_int(11);
    let product = field.mul_ref(&a, &b);
    let quotient = field.div(&product, &b);
    assert!(field.eq_el(&quotient, &a));
}

#[test]
fn test_powers() {
    let field_base = ArkFieldWrapper::<BlsFq>::new();

    let base = BlsFq::from(3u64);

    // Test x^0 = 1
    let result = base.pow([0u64]);
    assert!(result.is_one());

    // Test x^1 = x
    let result = base.pow([1u64]);
    assert_eq!(result, base);

    // Test x^2
    let result = base.pow([2u64]);
    let expected = base * base;
    assert_eq!(result, expected);

    // Test x^5
    let result = base.pow([5u64]);
    let expected = BlsFq::from(243u64); // 3^5 = 243
    assert_eq!(result, expected);
}

#[test]
fn test_characteristic() {
    use ark_feanor::prime_field::FieldProperties;

    let bn_field = ArkFieldWrapper::<BnFr>::new();
    let bls_field = ArkFieldWrapper::<BlsFr>::new();

    // Get characteristics
    let bn_char = bn_field.field_characteristic_biguint();
    let bls_char = bls_field.field_characteristic_biguint();

    // Verify they're different and non-zero
    assert!(bn_char != bls_char);
    assert!(bn_char > num_bigint::BigUint::from(0u32));
    assert!(bls_char > num_bigint::BigUint::from(0u32));

    // Verify they're prime field
    assert!(bn_field.is_prime_field());
    assert!(bls_field.is_prime_field());
}

#[test]
fn test_from_int_range() {
    let field = ArkFieldWrapper::<BnFr>::new();
    
    // Test positive values
    for i in 0..100 {
        let el = field.from_int(i);
        if i == 0 {
            assert!(field.is_zero(&el));
        } else if i == 1 {
            assert!(field.is_one(&el));
        }
    }
    
    // Test negative values
    for i in -100..0 {
        let el = field.from_int(i);
        let pos = field.from_int(-i);
        let sum = field.add_ref(&el, &pos);
        assert!(field.is_zero(&sum));
    }
}

#[test]
fn test_field_arithmetic_properties() {
    let field = ArkFieldWrapper::<BlsFr>::new();
    
    // Test that every non-zero element has an inverse
    for i in 1..20 {
        let el = field.from_int(i);
        let inv = field.div(&field.one(), &el);
        let product = field.mul_ref(&el, &inv);
        assert!(field.is_one(&product));
    }
    
    // Test field properties
    assert!(field.is_commutative());
    assert!(field.is_noetherian());
    assert!(!field.is_approximate());
}

#[test]
fn test_zero_and_one() {
    let field = ArkFieldWrapper::<BnFr>::new();
    
    let zero = field.zero();
    let one = field.one();
    
    // Basic checks
    assert!(field.is_zero(&zero));
    assert!(!field.is_zero(&one));
    assert!(field.is_one(&one));
    assert!(!field.is_one(&zero));
    
    // Arithmetic properties
    let sum = field.add_ref(&zero, &zero);
    assert!(field.is_zero(&sum));
    
    let product = field.mul_ref(&one, &one);
    assert!(field.is_one(&product));
    
    let product = field.mul_ref(&zero, &one);
    assert!(field.is_zero(&product));
}

#[test]
fn test_subtract() {
    let field = ArkFieldWrapper::<BlsFq>::new();
    
    let a = field.from_int(10);
    let b = field.from_int(3);
    
    // a - b = a + (-b)
    let neg_b = field.negate(field.clone_el(&b));
    let diff = field.add_ref(&a, &neg_b);
    let expected = field.from_int(7);
    assert!(field.eq_el(&diff, &expected));
}

#[test]
#[should_panic(expected = "Division by zero")]
fn test_division_by_zero_panics() {
    let field = ArkFieldWrapper::<BnFr>::new();
    let a = field.from_int(5);
    let _ = field.div(&a, &field.zero());
}

#[test]
fn test_is_unit() {
    let field = ArkFieldWrapper::<BlsFr>::new();
    
    // In a field, every non-zero element is a unit
    assert!(!field.is_unit(&field.zero()));
    assert!(field.is_unit(&field.one()));
    assert!(field.is_unit(&field.from_int(2)));
    assert!(field.is_unit(&field.from_int(-1)));
    assert!(field.is_unit(&field.from_int(12345)));
}
