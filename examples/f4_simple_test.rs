//! Simplest possible F4 test: single variable

#![feature(allocator_api)]

use ark_feanor::*;
use ark_feanor::f4::f4_simple;
use ark_bn254::Fr as BnFr;
use feanor_math::rings::multivariate::{DegRevLex};
use feanor_math::rings::multivariate::multivariate_impl::MultivariatePolyRingImpl;
use feanor_math::ring::RingStore;
use std::alloc::Global;
use std::time::Instant;

type BN254PolyRing = MultivariatePolyRingImpl<&'static RingValue<ArkFieldWrapper<BnFr>>, Global>;

/// Create trivial system: x^2 - 1, which should give basis {x^2 - 1} or {x - 1, x + 1}
fn create_trivial_system(poly_ring: &BN254PolyRing) -> Vec<El<BN254PolyRing>> {
    let mut system = Vec::new();

    // x^2 - 1
    let f1 = poly_ring.from_terms([
        (poly_ring.base_ring().one(), poly_ring.create_monomial([2])),
        (poly_ring.base_ring().negate(poly_ring.base_ring().one()), poly_ring.create_monomial([0])),
    ].into_iter());

    system.push(f1);
    system
}

fn main() {
    let field = &*BN254_FR;
    let n_vars = 1;

    let poly_ring = MultivariatePolyRingImpl::new(field, n_vars);

    println!("Testing F4 on trivial 1-variable system");
    println!("System: x^2 - 1");
    println!();

    // Run F4
    println!("Running F4...");
    let system_f4 = create_trivial_system(&poly_ring);
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

    println!("Result:");
    println!("  Gröbner basis has {} element(s)", gb_f4.len());

    // Print the basis elements
    for (i, _poly) in gb_f4.iter().enumerate() {
        println!("  G[{}]", i);
    }
}
