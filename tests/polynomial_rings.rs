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

// Multivariate polynomial tests

#[test]
fn test_multivariate_polynomial_creation() {
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    // Create polynomial: x^2 + xy + y^2
    let [poly] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone().pow(2) + x.clone() * y.clone() + y.clone().pow(2)]
    });

    assert!(!poly_ring.is_zero(&poly));
    // Verify this is a non-trivial polynomial
    assert!(!poly_ring.is_one(&poly));
}

#[test]
fn test_multivariate_polynomial_arithmetic() {
    let field = &*BLS12_381_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    let [p, q] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [
            x.clone() + y.clone(),  // p = x + y
            x.clone() - y.clone(),  // q = x - y
        ]
    });

    // Test addition: (x + y) + (x - y) = 2x
    let sum = poly_ring.add_ref(&p, &q);
    let [expected] = poly_ring.with_wrapped_indeterminates(|[x, _y]| {
        [x.clone() + x.clone()]
    });
    assert!(poly_ring.eq_el(&sum, &expected));

    // Test multiplication: (x + y)(x - y) = x^2 - y^2
    let product = poly_ring.mul_ref(&p, &q);
    let [expected_product] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone().pow(2) - y.clone().pow(2)]
    });
    assert!(poly_ring.eq_el(&product, &expected_product));
}

#[test]
fn test_multivariate_three_variables() {
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 3);

    // Create polynomial: xyz + x^2 + y^2 + z^2
    let [poly] = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        [x.clone() * y.clone() * z.clone()
            + x.clone().pow(2)
            + y.clone().pow(2)
            + z.clone().pow(2)]
    });

    assert!(!poly_ring.is_zero(&poly));
    assert!(!poly_ring.is_one(&poly));

    // Test that adding zero doesn't change the polynomial
    let zero = poly_ring.zero();
    let sum = poly_ring.add_ref(&poly, &zero);
    assert!(poly_ring.eq_el(&sum, &poly));
}

#[test]
fn test_multivariate_constant_polynomial() {
    let field = &*BLS12_381_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    // Test polynomial identity: (x + y) * 1 = x + y
    let [linear] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone() + y.clone()]
    });

    let one = poly_ring.one();
    let product = poly_ring.mul_ref(&one, &linear);
    assert!(poly_ring.eq_el(&product, &linear));

    // Test zero: (x + y) * 0 = 0
    let zero = poly_ring.zero();
    let product_with_zero = poly_ring.mul_ref(&zero, &linear);
    assert!(poly_ring.is_zero(&product_with_zero));
}

#[test]
fn test_multivariate_homogeneous_polynomial() {
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 3);

    // Create homogeneous polynomial of degree 2: x^2 + xy + xz + y^2 + yz + z^2
    let [poly] = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        [x.clone().pow(2)
            + x.clone() * y.clone()
            + x.clone() * z.clone()
            + y.clone().pow(2)
            + y.clone() * z.clone()
            + z.clone().pow(2)]
    });

    assert!(!poly_ring.is_zero(&poly));
    assert!(!poly_ring.is_one(&poly));

    // Verify it's a non-trivial homogeneous polynomial by checking it has multiple terms
    let term_count = poly_ring.terms(&poly).count();
    assert!(term_count > 1, "Should have multiple terms");
}

#[test]
fn test_multivariate_zero_and_one() {
    let field = &*BLS12_381_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    let zero = poly_ring.zero();
    let one = poly_ring.one();

    assert!(poly_ring.is_zero(&zero));
    assert!(poly_ring.is_one(&one));
    assert!(!poly_ring.is_zero(&one));
    assert!(!poly_ring.is_one(&zero));

    // Test that 1 * p = p
    let [p] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone() * y.clone().pow(2)]
    });

    let product = poly_ring.mul_ref(&one, &p);
    assert!(poly_ring.eq_el(&product, &p));
}

#[test]
fn test_multivariate_negation() {
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);

    // Create polynomial: x + y
    let [p] = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        [x.clone() + y.clone()]
    });

    // Negate it: -(x + y) = -x - y
    let neg_p = poly_ring.negate(poly_ring.clone_el(&p));

    // Add them: should get zero
    let sum = poly_ring.add_ref(&p, &neg_p);
    assert!(poly_ring.is_zero(&sum));
}
