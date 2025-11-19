//! F4 scaling test
//! Tests F4 performance on progressively larger systems

#![feature(allocator_api)]

use ark_feanor::*;
use ark_feanor::f4::f4_simple;
use ark_bn254::Fr as BnFr;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::rings::multivariate::{DegRevLex};
use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
use feanor_math::ring::RingStore;
use std::alloc::Global;
use std::time::Instant;

type BN254PolyRing = MultivariatePolyRingImpl<&'static RingValue<ArkFieldWrapper<BnFr>>, Global>;

/// Create Katsura-n polynomial system
fn create_katsura_system(poly_ring: &BN254PolyRing, n: usize) -> Vec<El<BN254PolyRing>> {
    let mut system = Vec::new();

    // Create monomials for variables
    let mut vars = Vec::new();
    for i in 0..n {
        let mut exponents = vec![0; n];
        exponents[i] = 1;
        vars.push(poly_ring.from_terms([(
            poly_ring.base_ring().one(),
            poly_ring.create_monomial(exponents.into_iter())
        )].into_iter()));
    }

    // u_0 + 2*sum u_i = 1
    let mut eq0 = poly_ring.clone_el(&vars[0]);
    for i in 1..n {
        let two = poly_ring.base_ring().int_hom().map(2);
        let term = poly_ring.inclusion().map(two);
        let scaled = poly_ring.mul_ref(&term, &vars[i]);
        eq0 = poly_ring.add_ref(&eq0, &scaled);
    }
    let one = poly_ring.one();
    system.push(poly_ring.add_ref(&eq0, &poly_ring.negate(one)));

    // Katsura equations
    for k in 0..(n-1) {
        let mut poly = poly_ring.zero();

        for i in 0..n {
            for j in 0..n {
                let sum_abs = if i >= j { i - j } else { j - i };
                if sum_abs == k {
                    let term = poly_ring.mul_ref(&vars[i], &vars[j]);
                    poly = poly_ring.add_ref(&poly, &term);
                }
            }
        }

        poly = poly_ring.sub_ref(&poly, &vars[k]);
        system.push(poly);
    }

    system
}

fn run_test(n: usize) {
    let field = &*BN254_FR;

    let poly_ring = MultivariatePolyRingImpl::new(field, n);

    println!("==== Katsura-{} ====", n);

    // Run F4
    println!("Running F4...");
    let system_f4 = create_katsura_system(&poly_ring, n);
    let start = Instant::now();
    let gb_f4 = f4_simple::<_, DegRevLex>(
        &poly_ring,
        system_f4,
        DegRevLex
    );
    let time_f4 = start.elapsed();

    println!("  Time: {:?}", time_f4);
    println!("  Basis size: {}", gb_f4.len());
    println!();
}

fn main() {
    println!("F4 Scaling Test");
    println!();

    // Test on progressively larger systems
    run_test(5);  // medium
    run_test(6);  // larger
    run_test(7);  // large
    run_test(9);  // very large

    println!("Test complete!");
}
