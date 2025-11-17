//! Microbenchmark: Raw Field Operations WITH AVX512-IFMA Batch Operations
//!
//! Tests AVX512-IFMA batch operations directly to measure their performance potential.
//! Compares sequential operations vs batched AVX512-IFMA operations.
//!
//! **IMPORTANT**: This benchmark demonstrates what AVX512 *could* achieve if integrated
//! into the Gröbner basis algorithms. Currently, F4/Buchberger don't use these ops.
//!
//! Operations tested:
//! - Sequential field ops (baseline)
//! - AVX512 batch-8 multiplication (mont_mul_batch_8)

#![feature(allocator_api)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_bn254::Fr as BnFr;
use ark_feanor::*;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::ring::RingStore;
use std::time::Duration;

#[cfg(all(feature = "avx512-ifma", target_arch = "x86_64"))]
use ark_ff::fields::models::fp::avx512_backend;

/// Benchmark sequential field multiplication (baseline)
fn bench_field_mul_sequential(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_mul_sequential");
    group.measurement_time(Duration::from_secs(10));

    let field = &*BN254_FR;
    let a: Vec<_> = (0..8).map(|i| field.int_hom().map(12345 + i)).collect();
    let b: Vec<_> = (0..8).map(|i| field.int_hom().map(67890 + i)).collect();

    group.bench_function("bn254_mul_8elem_100K_sequential", |bench| {
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

/// Benchmark AVX512-IFMA batch multiplication
#[cfg(all(feature = "avx512-ifma", target_arch = "x86_64"))]
fn bench_field_mul_avx512_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_mul_avx512");
    group.measurement_time(Duration::from_secs(10));

    // Create arrays of BnFr elements
    let a: [BnFr; 8] = [
        BnFr::from(12345u64),
        BnFr::from(12346u64),
        BnFr::from(12347u64),
        BnFr::from(12348u64),
        BnFr::from(12349u64),
        BnFr::from(12350u64),
        BnFr::from(12351u64),
        BnFr::from(12352u64),
    ];

    let b: [BnFr; 8] = [
        BnFr::from(67890u64),
        BnFr::from(67891u64),
        BnFr::from(67892u64),
        BnFr::from(67893u64),
        BnFr::from(67894u64),
        BnFr::from(67895u64),
        BnFr::from(67896u64),
        BnFr::from(67897u64),
    ];

    group.bench_function("bn254_mul_8elem_100K_avx512_batch", |bench| {
        bench.iter(|| {
            let mut result = [BnFr::from(0u64); 8];
            for _ in 0..100_000 {
                avx512_backend::mont_mul_batch_8(black_box(&a), black_box(&b), &mut result);
            }
            black_box(result)
        });
    });

    group.finish();
}

/// Fallback when AVX512-IFMA not available
#[cfg(not(all(feature = "avx512-ifma", target_arch = "x86_64")))]
fn bench_field_mul_avx512_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_mul_avx512");

    group.bench_function("avx512_not_available", |bench| {
        bench.iter(|| {
            black_box("AVX512-IFMA not enabled or not on x86_64")
        });
    });

    group.finish();
}

/// Benchmark to verify AVX512 feature detection
fn bench_avx512_detection(c: &mut Criterion) {
    c.bench_function("avx512_feature_check", |b| {
        b.iter(|| {
            #[cfg(all(feature = "avx512-ifma", target_arch = "x86_64"))]
            {
                black_box("AVX512-IFMA: ENABLED")
            }

            #[cfg(not(all(feature = "avx512-ifma", target_arch = "x86_64")))]
            {
                black_box("AVX512-IFMA: NOT ENABLED")
            }
        });
    });
}

criterion_group!(
    benches,
    bench_avx512_detection,
    bench_field_mul_sequential,
    bench_field_mul_avx512_batch,
);

criterion_main!(benches);
