//! Benchmark: Buchberger Algorithm Baseline (No AVX512)
//!
//! This benchmark tests the traditional Buchberger algorithm for computing Gröbner bases
//! over cryptographic prime fields (BN254, BLS12-381) using standard arkworks operations
//! WITHOUT AVX512 optimizations.
//!
//! Test systems:
//! - Cyclic-7: 7 variables, homogeneous system, degree 7
//! - Katsura-9: 9 variables, dense system from Katsura problem
//! - Custom large systems for stress testing

#![feature(allocator_api)]

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::algorithms::buchberger::{buchberger_simple, BuchbergerConfig};
use feanor_math::rings::multivariate::{DegRevLex, Lex};
use feanor_math::rings::multivariate::multivariate_impl::{DegreeCfg, MultivariatePolyRingImpl};
use feanor_math::ring::RingStore;
use std::alloc::Global;
use std::time::Duration;

type BN254PolyRing = MultivariatePolyRingImpl<&'static RingValue<ArkFieldWrapper<BnFr>>, Global>;

/// Create cyclic-n polynomial system
/// Cyclic-n: {sum(x_i) = 0, sum(x_i * x_{i+1}) = 0, ..., prod(x_i) - 1 = 0}
fn create_cyclic_system(poly_ring: &BN254PolyRing, n: usize) -> Vec<El<BN254PolyRing>> {
    assert!(n >= 3, "Cyclic system requires at least 3 variables");

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

    // Create cyclic equations: sum of products of k consecutive variables
    for k in 1..n {
        let mut poly = poly_ring.zero();
        for i in 0..n {
            let mut term = poly_ring.one();
            for j in 0..k {
                let idx = (i + j) % n;
                term = poly_ring.mul_ref(&term, &vars[idx]);
            }
            poly = poly_ring.add_ref(&poly, &term);
        }
        system.push(poly);
    }

    // Final equation: product of all variables minus 1
    let mut prod = poly_ring.one();
    for var in &vars {
        prod = poly_ring.mul_ref(&prod, var);
    }
    let neg_one = poly_ring.negate(poly_ring.one());
    system.push(poly_ring.add_ref(&prod, &neg_one));

    system
}

/// Create Katsura-n polynomial system
/// Katsura-n: System arising from magnetism problems
fn create_katsura_system(poly_ring: &BN254PolyRing, n: usize) -> Vec<El<BN254PolyRing>> {
    assert!(n >= 3, "Katsura system requires at least 3 variables");

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

    // u_i = u_{-i} (symmetry), u_0 + 2*sum_{i=1}^{n-1} u_i = 1
    let mut eq0 = poly_ring.clone_el(&vars[0]);
    for i in 1..n {
        let two = poly_ring.base_ring().int_hom().map(2);
        let term = poly_ring.inclusion().map(two);
        let scaled = poly_ring.mul_ref(&term, &vars[i]);
        eq0 = poly_ring.add_ref(&eq0, &scaled);
    }
    let one = poly_ring.one();
    system.push(poly_ring.add_ref(&eq0, &poly_ring.negate(one)));

    // sum_{j=-n+1}^{n-1} u_i u_j = u_k for each k
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

        // Subtract u_k
        poly = poly_ring.sub_ref(&poly, &vars[k]);
        system.push(poly);
    }

    system
}

/// Benchmark Buchberger on Cyclic-7 system
fn bench_buchberger_cyclic7(c: &mut Criterion) {
    let mut group = c.benchmark_group("buchberger_cyclic7");
    group.sample_size(10); // Reduce sample size for long-running tests
    group.measurement_time(Duration::from_secs(120)); // 2 minutes max per benchmark

    let field = &*BN254_FR;
    let n_vars = 7;

    // Use degree configuration - precompute up to max degree for binomial tables
    let degree_cfg = DegreeCfg::new(20).with_precompute(20);

    group.bench_function("cyclic7_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_cyclic_system(&poly_ring, n_vars);
            let result = buchberger_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark Buchberger on Katsura-7 system (reduced from 9 for practical runtime)
fn bench_buchberger_katsura7(c: &mut Criterion) {
    let mut group = c.benchmark_group("buchberger_katsura7");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(120));

    let field = &*BN254_FR;
    let n_vars = 7;

    let degree_cfg = DegreeCfg::new(15).with_precompute(15);

    group.bench_function("katsura7_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_katsura_system(&poly_ring, n_vars);
            let result = buchberger_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark Buchberger on larger Cyclic systems
fn bench_buchberger_cyclic_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("buchberger_cyclic_scaling");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(180));

    let field = &*BN254_FR;

    for n_vars in [5, 6, 7] {
        let degree_cfg = DegreeCfg::new(20).with_precompute(20);

        group.bench_with_input(
            BenchmarkId::from_parameter(n_vars),
            &n_vars,
            |b, &n| {
                let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
                    field, n, degree_cfg, (0, 0), Global
                );

                b.iter(|| {
                    let system = create_cyclic_system(&poly_ring, n);
                    let result = buchberger_simple::<_, DegRevLex>(
                        &poly_ring,
                        system,
                        DegRevLex
                    );
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_buchberger_cyclic7,
    bench_buchberger_katsura7,
    bench_buchberger_cyclic_scaling,
);

criterion_main!(benches);
