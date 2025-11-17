//! Microbenchmark: Raw Field Operations (No AVX512)
//!
//! Tests pure field arithmetic performance to isolate computational overhead
//! from algorithmic complexity. This provides a baseline for comparison with AVX512.
//!
//! Operations tested:
//! - Field multiplication (1M operations)
//! - Field addition (1M operations)
//! - Field squaring (1M operations)
//! - Field inversion (10K operations - expensive)

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use ark_feanor::*;
use ark_bn254::Fr as BnFr;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::ring::RingStore;
use std::time::Duration;

const ITER_MUL: usize = 1_000_000;
const ITER_ADD: usize = 1_000_000;
const ITER_SQR: usize = 1_000_000;
const ITER_INV: usize = 10_000;

/// Benchmark field multiplication
fn bench_field_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_mul");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let a = field.int_hom().map(12345);
    let b = field.int_hom().map(67890);

    group.bench_function("bn254_mul_1M", |bench| {
        bench.iter(|| {
            let mut result = field.clone_el(&a);
            for _ in 0..ITER_MUL {
                result = field.mul_ref(&result, black_box(&b));
            }
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark field addition
fn bench_field_add(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_add");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let a = field.int_hom().map(12345);
    let b = field.int_hom().map(67890);

    group.bench_function("bn254_add_1M", |bench| {
        bench.iter(|| {
            let mut result = field.clone_el(&a);
            for _ in 0..ITER_ADD {
                result = field.add_ref(&result, black_box(&b));
            }
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark field squaring
fn bench_field_sqr(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_sqr");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let a = field.int_hom().map(12345);

    group.bench_function("bn254_sqr_1M", |bench| {
        bench.iter(|| {
            let mut result = field.clone_el(&a);
            for _ in 0..ITER_SQR {
                result = field.mul_ref(black_box(&result), black_box(&result));
            }
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark field inversion (expensive operation)
fn bench_field_inv(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_inv");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let a = field.int_hom().map(12345);

    group.bench_function("bn254_inv_10K", |bench| {
        bench.iter(|| {
            let mut result = field.clone_el(&a);
            for _ in 0..ITER_INV {
                result = field.invert(black_box(&result)).unwrap();
            }
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark batch operations (simulate potential AVX512 usage pattern)
fn bench_field_batch_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_batch");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;

    // Create arrays of 8 elements (AVX512 batch size)
    let a: Vec<_> = (0..8).map(|i| field.int_hom().map(12345 + i)).collect();
    let b: Vec<_> = (0..8).map(|i| field.int_hom().map(67890 + i)).collect();

    group.bench_function("bn254_mul_batch8_100K", |bench| {
        bench.iter(|| {
            let mut results = Vec::with_capacity(8);
            for _ in 0..100_000 {
                for i in 0..8 {
                    results.push(field.mul_ref(black_box(&a[i]), black_box(&b[i])));
                }
                results.clear();
            }
            black_box(results)
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_field_mul,
    bench_field_add,
    bench_field_sqr,
    bench_field_inv,
    bench_field_batch_ops,
);

criterion_main!(benches);
