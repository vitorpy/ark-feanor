//! F4 vs Buchberger scaling comparison
//! Tests on progressively larger systems

#![feature(allocator_api)]

use ark_feanor::*;
use ark_feanor::f4::f4_simple;
use ark_bn254::Fr as BnFr;
use feanor_math::algorithms::buchberger::buchberger_simple;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::rings::multivariate::{DegRevLex};
use feanor_math::rings::multivariate::multivariate_impl::{DegreeCfg, MultivariatePolyRingImpl};
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

fn run_comparison(n: usize, skip_buchberger: bool) {
    let field = &*BN254_FR;

    // Adjust degree limits based on problem size
    let max_deg = match n {
        3 => 50,
        4 => 100,
        5 => 200,
        6 => 500,
        7 => 800,
        _ => 1200,
    };

    let degree_cfg = DegreeCfg::new(max_deg).with_precompute(50);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
        field, n, degree_cfg, (0, 0), Global
    );

    println!("==== Katsura-{} ====", n);

    // Run Buchberger (unless skipped for large problems)
    let (gb_buchberger, time_buchberger) = if !skip_buchberger {
        println!("Running Buchberger...");
        let system_buch = create_katsura_system(&poly_ring, n);
        let start = Instant::now();
        let gb = buchberger_simple::<_, DegRevLex>(
            &poly_ring,
            system_buch,
            DegRevLex
        );
        let time = start.elapsed();

        println!("  Time: {:?}", time);
        println!("  Basis size: {}", gb.len());
        println!();
        (Some(gb), Some(time))
    } else {
        println!("Skipping Buchberger (too slow for this size)...\n");
        (None, None)
    };

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

    // Compare
    println!("Results:");
    if let Some(ref gb_buch) = gb_buchberger {
        println!("  Basis sizes: Buchberger={}, F4={}",
                 gb_buch.len(), gb_f4.len());
    } else {
        println!("  Basis size: F4={}",
                 gb_f4.len());
    }

    // Timing comparison
    if let Some(time_buch) = time_buchberger {
        let ratio = time_f4.as_secs_f64() / time_buch.as_secs_f64();

        println!("\n  Buchberger vs F4:");
        if ratio < 1.0 {
            println!("    F4 is {:.2}x FASTER", 1.0 / ratio);
        } else {
            println!("    F4 is {:.2}x SLOWER", ratio);
        }
    }
    println!();
}

fn main() {
    println!("F4 vs Buchberger Scaling Test");
    println!();

    // Size 5 (medium)
    //run_comparison(5, false);

    // Size 6 (larger)
    //run_comparison(6, false);

    // Size 7 (large)
    run_comparison(9, false);

    println!("Test complete!");
}
