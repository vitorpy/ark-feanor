//! Benchmark: F4 Algorithm on Tiny Systems (No AVX512)
//!
//! Tests the F4 algorithm on very small polynomial systems to measure
//! F4's overhead on tiny problems where batch processing may not help.
//!
//! Test systems:
//! - Cyclic-3: 3 variables (fast, ~100ms)
//! - Cyclic-4: 4 variables (moderate, ~500ms)
//! - Katsura-3: 3 variables (fast, ~50ms)
//! - Katsura-4: 4 variables (moderate, ~200ms)

#![feature(allocator_api)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_feanor::*;
use ark_feanor::f4::f4_simple;
use ark_bn254::Fr as BnFr;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::rings::multivariate::{DegRevLex};
use feanor_math::rings::multivariate::multivariate_impl::{DegreeCfg, MultivariatePolyRingImpl};
use feanor_math::ring::RingStore;
use std::alloc::Global;
use std::time::Duration;

type BN254PolyRing = MultivariatePolyRingImpl<&'static RingValue<ArkFieldWrapper<BnFr>>, Global>;

/// Create cyclic-n polynomial system
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

    // Create cyclic equations
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

    // Final equation: product - 1
    let mut prod = poly_ring.one();
    for var in &vars {
        prod = poly_ring.mul_ref(&prod, var);
    }
    let neg_one = poly_ring.negate(poly_ring.one());
    system.push(poly_ring.add_ref(&prod, &neg_one));

    system
}

/// Create Katsura-n polynomial system
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

/// Benchmark F4 on Cyclic-3 (very fast)
fn bench_f4_cyclic3(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic3");
    group.sample_size(50);
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let n_vars = 3;
    let degree_cfg = DegreeCfg::new(100).with_precompute(50);

    group.bench_function("cyclic3_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_cyclic_system(&poly_ring, n_vars);
            let result = f4_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark F4 on Cyclic-4
fn bench_f4_cyclic4(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic4");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(15));

    let field = &*BN254_FR;
    let n_vars = 4;
    let degree_cfg = DegreeCfg::new(200).with_precompute(50);

    group.bench_function("cyclic4_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_cyclic_system(&poly_ring, n_vars);
            let result = f4_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark F4 on Katsura-3 (very fast)
fn bench_f4_katsura3(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_katsura3");
    group.sample_size(50);
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let n_vars = 3;
    let degree_cfg = DegreeCfg::new(50).with_precompute(30);

    group.bench_function("katsura3_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_katsura_system(&poly_ring, n_vars);
            let result = f4_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark F4 on Katsura-4
fn bench_f4_katsura4(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_katsura4");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(15));

    let field = &*BN254_FR;
    let n_vars = 4;
    let degree_cfg = DegreeCfg::new(100).with_precompute(30);

    group.bench_function("katsura4_degrevlex", |b| {
        let poly_ring = MultivariatePolyRingImpl::new_with_mult_table_ex(
            field, n_vars, degree_cfg, (0, 0), Global
        );

        b.iter(|| {
            let system = create_katsura_system(&poly_ring, n_vars);
            let result = f4_simple::<_, DegRevLex>(
                &poly_ring,
                system,
                DegRevLex
            );
            black_box(result)
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_f4_cyclic3,
    bench_f4_cyclic4,
    bench_f4_katsura3,
    bench_f4_katsura4,
);

criterion_main!(benches);
