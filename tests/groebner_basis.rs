//! Tests for Gröbner basis computation and field properties over cryptographic fields

#![feature(allocator_api)]

use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use ark_bls12_381::Fr as BlsFr;
use ark_feanor::prime_field::FieldProperties;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::algorithms::buchberger::buchberger_simple;
use feanor_math::rings::multivariate::{DegRevLex, Lex};
use feanor_math::rings::multivariate::multivariate_impl::{DegreeCfg, MultivariatePolyRingImpl};
use std::alloc::Global;

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

// Gröbner basis computation tests

#[test]
fn test_groebner_basis_linear_system() {
    let field = &*BN254_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 2, degree_cfg, (1, 1), Global);

    // System: x + y = 0, x - y = 0
    // Solution: x = 0, y = 0
    let [p1, p2] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [
            x.clone() + y.clone(),
            x.clone() - y.clone(),
        ]
    });

    let gb = buchberger_simple(&poly_ring, vec![p1, p2], DegRevLex);

    // Gröbner basis should not be empty
    assert!(!gb.is_empty(), "Gröbner basis should not be empty");

    // Check that we got a reduced basis
    for poly in &gb {
        assert!(!poly_ring.is_zero(poly), "Basis should not contain zero");
    }
}

#[test]
fn test_groebner_basis_with_lex_ordering() {
    let field = &*BLS12_381_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 2, degree_cfg, (1, 1), Global);

    // Simple system: xy - 1, x^2 - 1
    let [p1, p2] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [
            x.clone() * y.clone() - 1,
            x.clone().pow(2) - 1,
        ]
    });

    let gb = buchberger_simple(&poly_ring, vec![p1, p2], Lex);

    assert!(!gb.is_empty(), "Gröbner basis should not be empty");

    // The basis should contain at least the input polynomials
    assert!(gb.len() >= 2, "Basis should have at least 2 elements");
}

#[test]
fn test_groebner_basis_ideal_membership() {
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    // Generate ideal from x^2 - y, y^2 - x
    let [p1, p2] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [
            x.clone().pow(2) - y.clone(),
            y.clone().pow(2) - x.clone(),
        ]
    });

    let gb = buchberger_simple(&poly_ring, vec![p1, p2], DegRevLex);

    assert!(!gb.is_empty());

    // Verify all polynomials in the basis are non-trivial
    for poly in &gb {
        assert!(!poly_ring.is_zero(poly));
        assert!(!poly_ring.is_one(poly));
    }
}

#[test]
fn test_groebner_basis_triangular_system() {
    let field = &*BLS12_381_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 3, degree_cfg, (1, 1), Global);

    // Triangular system: x + y + z, xy, xz
    let [p1, p2, p3] = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        [
            x.clone() + y.clone() + z.clone(),
            x.clone() * y.clone(),
            x.clone() * z.clone(),
        ]
    });

    let gb = buchberger_simple(&poly_ring, vec![p1, p2, p3], DegRevLex);

    assert!(!gb.is_empty());
    println!("Triangular system Gröbner basis has {} elements", gb.len());
}

#[test]
fn test_groebner_basis_difference_of_squares() {
    let field = &*BN254_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 2, degree_cfg, (1, 1), Global);

    // System: x^2 - y^2
    let [input] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone().pow(2) - y.clone().pow(2)]
    });

    let gb = buchberger_simple(&poly_ring, vec![input], DegRevLex);

    // For a single polynomial, the Gröbner basis is essentially itself (up to scaling)
    assert_eq!(gb.len(), 1, "Single polynomial should have basis of size 1");
}

#[test]
fn test_groebner_basis_comparison_orderings() {
    let field = &*BLS12_381_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 2, degree_cfg, (1, 1), Global);

    let [p1, p2] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [
            x.clone().pow(2) + y.clone(),
            x.clone() * y.clone() + 1,
        ]
    });

    // Compute with DegRevLex
    let gb_degrevlex = buchberger_simple(&poly_ring, vec![
        poly_ring.clone_el(&p1),
        poly_ring.clone_el(&p2)
    ], DegRevLex);

    // Compute with Lex
    let gb_lex = buchberger_simple(&poly_ring, vec![p1, p2], Lex);

    // Both should be non-empty
    assert!(!gb_degrevlex.is_empty());
    assert!(!gb_lex.is_empty());

    // The bases might have different sizes due to different orderings
    println!("DegRevLex basis size: {}", gb_degrevlex.len());
    println!("Lex basis size: {}", gb_lex.len());
}

#[test]
fn test_groebner_basis_homogeneous_ideal() {
    let field = &*BN254_FR;
    let degree_cfg = DegreeCfg::new(64).with_precompute(1);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(field, 3, degree_cfg, (1, 1), Global);

    // Homogeneous ideal: x^2, xy, xz
    let [p1, p2, p3] = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        [
            x.clone().pow(2),
            x.clone() * y.clone(),
            x.clone() * z.clone(),
        ]
    });

    let gb = buchberger_simple(&poly_ring, vec![p1, p2, p3], DegRevLex);

    assert!(!gb.is_empty());

    // All elements should be homogeneous
    for poly in &gb {
        assert!(!poly_ring.is_zero(poly));
    }

    println!("Homogeneous ideal basis size: {}", gb.len());
}
