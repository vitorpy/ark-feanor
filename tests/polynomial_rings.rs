//! Tests for polynomial ring creation and operations

use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use ark_bls12_381::Fr as BlsFr;
use feanor_math::rings::poly::dense_poly::DensePolyRing;
use feanor_math::rings::poly::PolyRingStore;
use feanor_math::rings::multivariate::*;

#[test]
fn test_univariate_polynomial_creation() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = DensePolyRing::new(field, "x");
    
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
    let result = poly_ring.evaluate(&poly, &x, &field.identity());
    let expected = field.int_hom().map(17);
    assert!(field.eq_el(&result, &expected));
}

#[test]
fn test_polynomial_arithmetic() {
    let field = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    let poly_ring = DensePolyRing::new(field, "x");
    
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
fn test_multivariate_polynomial_creation() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = MultivariatePolyRingImpl::new(&field, 2);
    
    // Create polynomial: x^2 + y^2 + xy
    let x_squared = poly_ring.from_terms([
        (field.one(), Monomial::new([2, 0]))
    ].iter().cloned());
    
    let y_squared = poly_ring.from_terms([
        (field.one(), Monomial::new([0, 2]))
    ].iter().cloned());
    
    let xy = poly_ring.from_terms([
        (field.one(), Monomial::new([1, 1]))
    ].iter().cloned());
    
    let poly = poly_ring.add(poly_ring.add(x_squared, y_squared), xy);
    
    // Check number of terms
    let terms: Vec<_> = poly_ring.terms(&poly).collect();
    assert_eq!(terms.len(), 3);
}

#[test]
fn test_multivariate_with_indeterminates() {
    let field = &*BLS12_381_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 3);
    
    // Use the convenient wrapped indeterminates
    let polys = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        vec![
            // x + y + z
            x.clone() + y.clone() + z.clone(),
            // x*y + y*z + x*z
            x.clone() * y.clone() + y.clone() * z.clone() + x.clone() * z.clone(),
            // x*y*z
            x.clone() * y.clone() * z.clone(),
        ]
    });
    
    assert_eq!(polys.len(), 3);
    
    // Check the first polynomial (x + y + z) has 3 terms
    let terms: Vec<_> = poly_ring.terms(&polys[0]).collect();
    assert_eq!(terms.len(), 3);
}

#[test]
fn test_polynomial_evaluation() {
    let field = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    let poly_ring = MultivariatePolyRingImpl::new(&field, 2);
    
    // Create polynomial: x^2 + 2xy + y^2 = (x + y)^2
    let poly = poly_ring.from_terms([
        (field.one(), Monomial::new([2, 0])),      // x^2
        (field.int_hom().map(2), Monomial::new([1, 1])),  // 2xy
        (field.one(), Monomial::new([0, 2])),      // y^2
    ].iter().cloned());
    
    // Evaluate at x = 3, y = 2
    let values = [field.int_hom().map(3), field.int_hom().map(2)];
    let result = poly_ring.evaluate(&poly, values.iter().cloned(), &field.identity());
    
    // (3 + 2)^2 = 25
    let expected = field.int_hom().map(25);
    assert!(field.eq_el(&result, &expected));
}

#[test]
fn test_polynomial_degree() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = MultivariatePolyRingImpl::new(&field, 2);
    
    // Create polynomial: x^3*y^2 + x^2*y + x + 1
    let poly = poly_ring.from_terms([
        (field.one(), Monomial::new([3, 2])),  // x^3*y^2
        (field.one(), Monomial::new([2, 1])),  // x^2*y
        (field.one(), Monomial::new([1, 0])),  // x
        (field.one(), Monomial::new([0, 0])),  // 1
    ].iter().cloned());
    
    // Total degree should be 5 (from x^3*y^2)
    let terms: Vec<_> = poly_ring.terms(&poly).collect();
    let max_degree = terms.iter()
        .map(|(_, m)| m.total_degree())
        .max()
        .unwrap();
    assert_eq!(max_degree, 5);
}

#[test]
fn test_zero_polynomial() {
    let field = RingValue::from(ArkFieldWrapper::<BlsFr>::new());
    let poly_ring = DensePolyRing::new(field, "x");
    
    let zero = poly_ring.zero();
    assert!(poly_ring.is_zero(&zero));
    
    // Adding zero doesn't change polynomial
    let p = poly_ring.from_terms([(field.one(), 1)].iter().cloned());
    let sum = poly_ring.add_ref(&p, &zero);
    assert!(poly_ring.eq_el(&sum, &p));
}

#[test]
fn test_polynomial_subtraction() {
    let field = RingValue::from(ArkFieldWrapper::<BnFr>::new());
    let poly_ring = MultivariatePolyRingImpl::new(&field, 2);
    
    let polys = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        let p = x.clone() + y.clone();  // x + y
        let q = x.clone() - y.clone();  // x - y
        
        // (x + y) - (x - y) = 2y
        let diff = p - q;
        diff
    });
    
    // Check result is 2y
    let terms: Vec<_> = poly_ring.terms(&polys).collect();
    assert_eq!(terms.len(), 1);
    assert_eq!(terms[0].1.total_degree(), 1);
    assert!(field.eq_el(&terms[0].0, &field.int_hom().map(2)));
}
