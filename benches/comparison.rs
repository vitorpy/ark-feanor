//! Benchmarks comparing native arkworks operations vs wrapped operations

use ark_bn254::Fr;
use ark_feanor::ArkFieldWrapper;
use ark_ff::{Field, One};
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use feanor_math::ring::*;
use feanor_math::homomorphism::Homomorphism;
use feanor_math::field::Field as FeanorField;

fn bench_arkworks_addition(c: &mut Criterion) {
    let mut group = c.benchmark_group("addition");

    // Native arkworks
    group.bench_function("arkworks_native", |bencher| {
        let a = Fr::from(12345u64);
        let b = Fr::from(67890u64);
        bencher.iter(|| {
            let c = black_box(a) + black_box(b);
            black_box(c)
        })
    });

    // Through wrapper
    group.bench_function("ark_feanor_wrapper", |bencher| {
        let field = ArkFieldWrapper::<Fr>::new();
        let a = field.from_int(12345);
        let b = field.from_int(67890);
        bencher.iter(|| {
            let c = field.add_ref(&black_box(a), &black_box(b));
            black_box(c)
        })
    });

    group.finish();
}

fn bench_arkworks_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("multiplication");
    
    // Native arkworks
    group.bench_function("arkworks_native", |b| {
        let a = Fr::from(12345u64);
        let x = Fr::from(67890u64);
        b.iter(|| {
            let c = black_box(a) * black_box(x);
            black_box(c)
        })
    });
    
    // Through wrapper
    group.bench_function("ark_feanor_wrapper", |b| {
        let field = ArkFieldWrapper::<Fr>::new();
        let a = field.from_int(12345);
        let x = field.from_int(67890);
        b.iter(|| {
            let c = field.mul_ref(&black_box(a), &black_box(x));
            black_box(c)
        })
    });
    
    group.finish();
}

fn bench_arkworks_inversion(c: &mut Criterion) {
    let mut group = c.benchmark_group("inversion");
    
    // Native arkworks
    group.bench_function("arkworks_native", |b| {
        let a = Fr::from(12345u64);
        b.iter(|| {
            let inv = black_box(a).inverse().unwrap();
            black_box(inv)
        })
    });
    
    // Through wrapper
    group.bench_function("ark_feanor_wrapper", |b| {
        let field = ArkFieldWrapper::<Fr>::new();
        let a = field.from_int(12345);
        b.iter(|| {
            let inv = field.div(&field.one(), &black_box(a));
            black_box(inv)
        })
    });
    
    group.finish();
}

fn bench_batch_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_operations");
    
    for size in [10, 100, 1000].iter() {
        // Native arkworks
        group.bench_with_input(
            BenchmarkId::new("arkworks_native", size),
            size,
            |b, &size| {
                let elements: Vec<Fr> = (0..size)
                    .map(|i| Fr::from((i * 1337) as u64))
                    .collect();
                
                b.iter(|| {
                    let mut acc = Fr::one();
                    for el in &elements {
                        acc = acc * el + el;
                    }
                    black_box(acc)
                })
            }
        );
        
        // Through wrapper
        group.bench_with_input(
            BenchmarkId::new("ark_feanor_wrapper", size),
            size,
            |b, &size| {
                let field = ArkFieldWrapper::<Fr>::new();
                let elements: Vec<Fr> = (0..size)
                    .map(|i| field.from_int((i * 1337) as i32))
                    .collect();
                
                b.iter(|| {
                    let mut acc = field.one();
                    for el in &elements {
                        acc = field.add_ref(&field.mul_ref(&acc, el), el);
                    }
                    black_box(acc)
                })
            }
        );
    }
    
    group.finish();
}

fn bench_power_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("power");
    
    for exp in [10u64, 100, 1000].iter() {
        // Native arkworks
        group.bench_with_input(
            BenchmarkId::new("arkworks_native", exp),
            exp,
            |b, &exp| {
                let base = Fr::from(3u64);
                b.iter(|| {
                    let result = black_box(base).pow(&[exp]);
                    black_box(result)
                })
            }
        );
        
        // Through wrapper
        group.bench_with_input(
            BenchmarkId::new("ark_feanor_wrapper", exp),
            exp,
            |b, &exp| {
                let field = RingValue::from(ArkFieldWrapper::<Fr>::new());
                let base = field.int_hom().map(3);
                b.iter(|| {
                    let result = field.pow(black_box(base), exp as usize);
                    black_box(result)
                })
            }
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_arkworks_addition,
    bench_arkworks_multiplication,
    bench_arkworks_inversion,
    bench_batch_operations,
    bench_power_operations,
);

criterion_main!(benches);
