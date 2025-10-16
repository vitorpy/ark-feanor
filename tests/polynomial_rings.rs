//! Tests for polynomial ring creation and operations

use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use ark_bls12_381::Fr as BlsFr;
use feanor_math::rings::poly::dense_poly::DensePolyRing;
use feanor_math::rings::poly::PolyRingStore;
use feanor_math::homomorphism::Homomorphism;

#[test]
fn test_univariate_polynomial_creation() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = DensePolyRing::new(&field, "x");

    // Create polynomial: 3x^2 + 2x + 1
    let coeffs = vec![
        field.int_hom().map(1),
        field.int_hom().map(2),
        field.int_hom().map(3),
    ];
    let poly = poly_ring.from_terms(coeffs.iter().enumerate().map(|(i, c)| (c.clone(), i)));

    // Check degree
    assert_eq!(poly_ring.degree(&poly).unwrap(), 2);

    // Evaluate at x = 2: 3(4) + 2(2) + 1 = 12 + 4 + 1 = 17
    let x = field.int_hom().map(2);
    let hom = feanor_math::homomorphism::identity(&field);
    let result = poly_ring.evaluate(&poly, &x, &hom);
    let expected = field.int_hom().map(17);
    assert!(field.eq_el(&result, &expected));
}

#[test]
fn test_polynomial_arithmetic() {
    let field = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    let poly_ring = DensePolyRing::new(&field, "x");

    // p(x) = x + 1
    let p = poly_ring.from_terms([
        (field.one(), 0),
        (field.one(), 1),
    ].iter().cloned());

    // q(x) = x - 1
    let q = poly_ring.from_terms([
        (field.neg_one(), 0),
        (field.one(), 1),
    ].iter().cloned());

    // Test addition: (x + 1) + (x - 1) = 2x
    let sum = poly_ring.add_ref(&p, &q);
    let expected = poly_ring.from_terms([(field.int_hom().map(2), 1)].iter().cloned());
    assert!(poly_ring.eq_el(&sum, &expected));

    // Test multiplication: (x + 1)(x - 1) = x^2 - 1
    let product = poly_ring.mul_ref(&p, &q);
    let expected = poly_ring.from_terms([
        (field.neg_one(), 0),
        (field.one(), 2),
    ].iter().cloned());
    assert!(poly_ring.eq_el(&product, &expected));
}

#[test]
fn test_zero_polynomial() {
    let field = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    let poly_ring = DensePolyRing::new(&field, "x");

    let zero = poly_ring.zero();
    assert!(poly_ring.is_zero(&zero));

    // Adding zero doesn't change polynomial
    let p = poly_ring.from_terms([(field.one(), 1)].iter().cloned());
    let sum = poly_ring.add_ref(&p, &zero);
    assert!(poly_ring.eq_el(&sum, &p));
}

#[test]
fn test_polynomial_degree() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = DensePolyRing::new(&field, "x");

    // Create polynomials of different degrees
    let p1 = poly_ring.from_terms([(field.one(), 0)].iter().cloned());
    assert_eq!(poly_ring.degree(&p1).unwrap(), 0);

    let p2 = poly_ring.from_terms([(field.one(), 3)].iter().cloned());
    assert_eq!(poly_ring.degree(&p2).unwrap(), 3);

    let zero = poly_ring.zero();
    assert!(poly_ring.degree(&zero).is_none());
}

#[test]
fn test_common_fields() {
    // Test pre-configured fields
    let bn_field = &*BN254_FR;
    let bls_field = &*BLS12_381_FR;

    // Basic operations
    let a = bn_field.int_hom().map(5);
    let b = bn_field.int_hom().map(3);
    let c = bn_field.add_ref(&a, &b);
    let expected = bn_field.int_hom().map(8);
    assert!(bn_field.eq_el(&c, &expected));

    // Same for BLS field
    let x = bls_field.int_hom().map(10);
    let y = bls_field.int_hom().map(20);
    let z = bls_field.mul_ref(&x, &y);
    let expected = bls_field.int_hom().map(200);
    assert!(bls_field.eq_el(&z, &expected));
}
