//! Example of field operations over cryptographic fields
//!
//! Note: Gröbner basis computation in feanor-math v3 requires unstable APIs.
//! This example demonstrates basic field and polynomial operations instead.

use ark_feanor::*;
use feanor_math::rings::poly::dense_poly::DensePolyRing;
use feanor_math::rings::poly::PolyRingStore;
use feanor_math::homomorphism::Homomorphism;

fn main() {
    println!("=== Field Operations over Cryptographic Fields ===\n");

    // Use the BLS12-381 scalar field
    let field = &*BLS12_381_FR;

    println!("Working with BLS12-381 Fr field");
    println!("Field characteristic has {} bits\n", {
        use ark_feanor::prime_field::FieldProperties;
        let wrapper = ArkFieldWrapper::<ark_bls12_381::Fr>::new();
        wrapper.field_characteristic_biguint().bits()
    });

    // Basic field arithmetic
    println!("Basic field arithmetic:");
    let a = field.int_hom().map(17);
    let b = field.int_hom().map(23);
    let sum = field.add_ref(&a, &b);
    let product = field.mul_ref(&a, &b);

    println!("  a = 17");
    println!("  b = 23");
    println!("  a + b = {} (checking: {})", field.format(&sum), field.eq_el(&sum, &field.int_hom().map(40)));
    println!("  a * b = {} (checking: {})", field.format(&product), field.eq_el(&product, &field.int_hom().map(391)));
    println!();

    // Univariate polynomials
    println!("Univariate polynomial operations:");
    let poly_ring = DensePolyRing::new(field, "x");

    // Create polynomial: x^2 + 2x + 1 = (x + 1)^2
    let poly = poly_ring.from_terms([
        (field.one(), 0),
        (field.int_hom().map(2), 1),
        (field.one(), 2),
    ].iter().cloned());

    println!("  p(x) = x^2 + 2x + 1");
    println!("  degree = {}", poly_ring.degree(&poly).unwrap());

    // Polynomial multiplication
    let linear = poly_ring.from_terms([
        (field.one(), 0),
        (field.one(), 1),
    ].iter().cloned());

    let product_poly = poly_ring.mul_ref(&poly, &linear);
    println!("  p(x) * (x + 1) has degree {}", poly_ring.degree(&product_poly).unwrap());
    println!();

    // Field division
    println!("Field division:");
    use feanor_math::field::Field as FeanorField;
    let field_base = ArkFieldWrapper::<ark_bls12_381::Fr>::new();
    let x = field_base.from_int(100);
    let y = field_base.from_int(7);
    let quotient = field_base.div(&x, &y);
    let reconstructed = field_base.mul_ref(&quotient, &y);

    println!("  100 / 7 = quotient");
    println!("  quotient * 7 ≈ 100 (in field arithmetic)");
    println!("  Verification: {}", field_base.eq_el(&reconstructed, &x));
    println!();

    println!("Example completed successfully!");
}
