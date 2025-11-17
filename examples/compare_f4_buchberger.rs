//! Quick comparison: F4 vs Buchberger
//! Tests correctness and timing on a tiny system

#![feature(allocator_api)]

use ark_feanor::*;
use ark_feanor::f4::f4_simple;
use ark_bn254::Fr as BnFr;
use feanor_math::algorithms::buchberger::buchberger_simple;
use feanor_math::rings::multivariate::{DegRevLex};
use feanor_math::rings::multivariate::multivariate_impl::{DegreeCfg, MultivariatePolyRingImpl};
use feanor_math::ring::RingStore;
use std::alloc::Global;
use std::time::Instant;

type BN254PolyRing = MultivariatePolyRingImpl<&'static RingValue<ArkFieldWrapper<BnFr>>, Global>;

/// Create simple 2-variable system: x^2 + y^2 - 1, xy - 1
fn create_simple_system(poly_ring: &BN254PolyRing) -> Vec<El<BN254PolyRing>> {
    let mut system = Vec::new();

    // x^2 + y^2 - 1
    let f1 = poly_ring.from_terms([
        (poly_ring.base_ring().one(), poly_ring.create_monomial([2, 0])),
        (poly_ring.base_ring().one(), poly_ring.create_monomial([0, 2])),
        (poly_ring.base_ring().negate(poly_ring.base_ring().one()), poly_ring.create_monomial([0, 0])),
    ].into_iter());

    // xy - 1
    let f2 = poly_ring.from_terms([
        (poly_ring.base_ring().one(), poly_ring.create_monomial([1, 1])),
        (poly_ring.base_ring().negate(poly_ring.base_ring().one()), poly_ring.create_monomial([0, 0])),
    ].into_iter());

    system.push(f1);
    system.push(f2);
    system
}

fn main() {
    let field = &*BN254_FR;
    let n_vars = 2;

    // Create ring with high degree limit for F4
    let degree_cfg = DegreeCfg::new(100).with_precompute(50);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
        field, n_vars, degree_cfg, (0, 0), Global
    );

    println!("Testing F4 vs Buchberger on simple 2-variable system");
    println!("System: x^2 + y^2 - 1, xy - 1");
    println!();

    // Run Buchberger
    println!("Running Buchberger...");
    let system_buch = create_simple_system(&poly_ring);
    let start = Instant::now();
    let gb_buchberger = buchberger_simple::<_, DegRevLex>(
        &poly_ring,
        system_buch,
        DegRevLex
    );
    let time_buchberger = start.elapsed();

    println!("  Time: {:?}", time_buchberger);
    println!("  Basis size: {}", gb_buchberger.len());
    println!();

    // Run F4
    println!("Running F4...");
    let system_f4 = create_simple_system(&poly_ring);
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

    // Compare results
    println!("Comparison:");
    println!("  Buchberger basis size: {}", gb_buchberger.len());
    println!("  F4 basis size: {}", gb_f4.len());

    if gb_buchberger.len() == gb_f4.len() {
        println!("  ✓ Basis sizes match");
    } else {
        println!("  ✗ BASIS SIZES DIFFER!");
    }

    println!();
    println!("Timing:");
    println!("  Buchberger: {:?}", time_buchberger);
    println!("  F4:         {:?}", time_f4);

    let ratio = time_f4.as_secs_f64() / time_buchberger.as_secs_f64();
    if ratio < 1.0 {
        println!("  F4 is {:.2}x FASTER", 1.0 / ratio);
    } else {
        println!("  F4 is {:.2}x SLOWER", ratio);
    }
}
