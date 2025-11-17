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

fn run_comparison(n: usize) {
    let field = &*BN254_FR;

    // Adjust degree limits based on problem size
    let max_deg = match n {
        3 => 50,
        4 => 100,
        5 => 200,
        6 => 400,
        _ => 800,
    };

    let degree_cfg = DegreeCfg::new(max_deg).with_precompute(50);
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
        field, n, degree_cfg, (0, 0), Global
    );

    println!("==== Katsura-{} ====", n);

    // Run Buchberger
    println!("Running Buchberger...");
    let system_buch = create_katsura_system(&poly_ring, n);
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
    println!("  Basis sizes: Buchberger={}, F4={}", gb_buchberger.len(), gb_f4.len());

    if gb_buchberger.len() == gb_f4.len() {
        println!("  ✓ Basis sizes match");
    } else {
        println!("  ⚠ Basis sizes differ (this is OK if they generate the same ideal)");
    }

    // Check correctness: both bases should generate the same ideal
    // Test: every Buchberger element should reduce to zero by F4 basis
    let mut all_reduce = true;
    for poly in &gb_buchberger {
        let mut remainder = poly_ring.clone_el(poly);
        loop {
            let lt = match poly_ring.LT(&remainder, DegRevLex) {
                Some(lt) => lt,
                None => break, // Reduced to zero - good!
            };

            let mut reduced = false;
            for basis_elem in &gb_f4 {
                if let Some((basis_lc, basis_lm)) = poly_ring.LT(basis_elem, DegRevLex) {
                    if let Ok(mult) = poly_ring.monomial_div(poly_ring.clone_monomial(lt.1), basis_lm) {
                        let coeff_mult = poly_ring.base_ring().checked_div(&lt.0, &basis_lc).unwrap();
                        let mut scaled = poly_ring.clone_el(basis_elem);
                        poly_ring.mul_assign_monomial(&mut scaled, mult);
                        poly_ring.inclusion().mul_assign_map(&mut scaled, coeff_mult);
                        remainder = poly_ring.sub_ref(&remainder, &scaled);
                        reduced = true;
                        break;
                    }
                }
            }
            if !reduced {
                all_reduce = false;
                break;
            }
        }
        if !all_reduce {
            break;
        }
    }

    if all_reduce {
        println!("  ✓ Correctness verified: Buchberger basis ⊆ F4 ideal");
    } else {
        println!("  ✗ CORRECTNESS FAILED: Bases generate different ideals!");
    }

    let ratio = time_f4.as_secs_f64() / time_buchberger.as_secs_f64();
    if ratio < 1.0 {
        println!("  F4 is {:.2}x FASTER", 1.0 / ratio);
    } else {
        println!("  F4 is {:.2}x SLOWER", ratio);
    }
    println!();
}

fn main() {
    println!("F4 vs Buchberger Scaling Test");
    println!("Testing Katsura systems of increasing size\n");

    // Start with size 3 (known to work)
    run_comparison(3);

    // Size 4
    run_comparison(4);

    // Size 5 (getting larger)
    run_comparison(5);

    println!("Test complete!");
}
