//! Example of polynomial operations over cryptographic fields
//!
//! Note: Gröbner basis computation in feanor-math v3 requires unstable APIs.
//! This example demonstrates univariate polynomial operations instead.

use ark_feanor::*;
use feanor_math::rings::poly::dense_poly::DensePolyRing;
use feanor_math::rings::poly::PolyRingStore;
use feanor_math::homomorphism::Homomorphism;

fn main() {
    println!("=== Polynomial Operations over Cryptographic Fields ===\n");

    example_polynomial_division();
    println!("\n{}\n", "=".repeat(60));
    example_polynomial_evaluation();
    println!("\n{}\n", "=".repeat(60));
    example_polynomial_composition();
}

/// Divide polynomials and verify factorizations
fn example_polynomial_division() {
    println!("Example 1: Polynomial Division over BN254 Fr");
    println!("---------------------------------------------");

    let field = &*BN254_FR;
    let poly_ring = DensePolyRing::new(field, "x");

    // Create p(x) = x² - 5x + 6 = (x - 2)(x - 3)
    let p = poly_ring.from_terms([
        (field.int_hom().map(6), 0),
        (field.negate(field.int_hom().map(5)), 1),
        (field.one(), 2),
    ].iter().cloned());

    // Create factor q(x) = x - 2
    let q = poly_ring.from_terms([
        (field.negate(field.int_hom().map(2)), 0),
        (field.one(), 1),
    ].iter().cloned());

    println!("p(x) = x² - 5x + 6 = (x - 2)(x - 3)");
    println!("q(x) = x - 2");
    println!();

    // Verify division
    use feanor_math::divisibility::DivisibilityRingStore;
    if let Some(quotient) = poly_ring.checked_div(&p, &q) {
        println!("✓ q(x) divides p(x)");
        println!("Quotient degree: {}", poly_ring.degree(&quotient).unwrap_or(0));
        println!("Expected: degree 1 (should be x - 3)");

        // Verify by multiplication
        let product = poly_ring.mul_ref(&q, &quotient);
        if poly_ring.eq_el(&product, &p) {
            println!("✓ Verified: q(x) * quotient = p(x)");
        }
    } else {
        println!("✗ Division failed");
    }
}

/// Evaluate polynomials at specific points
fn example_polynomial_evaluation() {
    println!("Example 2: Polynomial Evaluation over BLS12-381 Fr");
    println!("--------------------------------------------------");

    let field = &*BLS12_381_FR;
    let poly_ring = DensePolyRing::new(field, "x");

    // Create p(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
    let p = poly_ring.from_terms([
        (field.negate(field.int_hom().map(6)), 0),
        (field.int_hom().map(11), 1),
        (field.negate(field.int_hom().map(6)), 2),
        (field.one(), 3),
    ].iter().cloned());

    println!("p(x) = x³ - 6x² + 11x - 6");
    println!("Roots: x = 1, 2, 3");
    println!();

    // Evaluate at the roots
    println!("Evaluating polynomial:");
    for i in 1..=5 {
        let x = field.int_hom().map(i);
        let result = poly_ring.evaluate(&p, &x, &field.identity());

        let is_root = field.is_zero(&result);
        let status = if is_root { "✓ ROOT" } else { "  " };

        println!("  p({}) = {} {}", i, field.format(&result), status);
    }
    println!();

    // Show degree properties
    let degree = poly_ring.degree(&p).unwrap();
    println!("Properties:");
    println!("  Degree: {}", degree);
    println!("  Leading coefficient: {}", field.format(poly_ring.lc(&p).unwrap()));
}

/// Compose polynomials
fn example_polynomial_composition() {
    println!("Example 3: Polynomial Composition over BN254 Fr");
    println!("-----------------------------------------------");

    let field = &*BN254_FR;
    let poly_ring = DensePolyRing::new(field, "x");

    // f(x) = x² + 1
    let f = poly_ring.from_terms([
        (field.one(), 0),
        (field.one(), 2),
    ].iter().cloned());

    // g(x) = x + 2
    let g = poly_ring.from_terms([
        (field.int_hom().map(2), 0),
        (field.one(), 1),
    ].iter().cloned());

    println!("f(x) = x² + 1");
    println!("g(x) = x + 2");
    println!();

    // Compute f(g(x)) = (x + 2)² + 1 = x² + 4x + 5
    let composition = poly_ring.evaluate(&f, &g, &poly_ring.inclusion());

    println!("Computing f(g(x)):");
    println!("  f(g(x)) = f(x + 2) = (x + 2)² + 1");
    println!("  Expected: x² + 4x + 5");
    println!();

    let degree = poly_ring.degree(&composition).unwrap();
    println!("Result degree: {}", degree);

    // Verify by evaluating at a point
    let test_point = field.int_hom().map(3);
    let direct = poly_ring.evaluate(&composition, &test_point, &field.identity());

    // f(g(3)) = f(5) = 25 + 1 = 26
    let g_at_3 = poly_ring.evaluate(&g, &test_point, &field.identity());
    let expected = poly_ring.evaluate(&f, &g_at_3, &field.identity());

    println!("Verification at x = 3:");
    println!("  f(g(3)) = {} (direct evaluation)", field.format(&direct));
    println!("  f(g(3)) = {} (step-by-step)", field.format(&expected));

    if field.eq_el(&direct, &expected) {
        println!("  ✓ Composition correct");
    }

    println!();
    println!("Polynomial operations completed successfully!");
}
