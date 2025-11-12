# Benchmark Instructions: Buchberger vs F4 with AVX512 IFMA

This document provides instructions for running and comparing the Buchberger baseline with the F4 algorithm using AVX512 IFMA optimizations.

## Prerequisites

### Hardware Requirements
- **CPU**: AMD Ryzen 5 7640U (Zen 4) or any CPU with AVX512F, AVX512DQ, and AVX512IFMA support
- **Verification**: Check CPU features with `grep avx512_ifma /proc/cpuinfo`

### Software Setup
- Rust toolchain (edition 2021+)
- Forked ark-algebra with AVX512 IFMA support at `../ark-algebra`
- feanor-math at `vendor/feanor-math`

## Build Configuration

The project is configured to use AVX512 IFMA via `.cargo/config.toml`:

```toml
[build]
rustflags = ["-C", "target-cpu=native"]
```

This enables all CPU features including AVX512 IFMA on supported hardware.

## Running Benchmarks

### 1. Buchberger Baseline (No AVX512 modifications)

```bash
# Run Buchberger baseline benchmarks
cargo bench --bench buchberger_baseline

# View results
ls -lh target/criterion/buchberger_*
```

**Test Systems:**
- `cyclic7_degrevlex`: Cyclic-7 system (7 variables, degree 7)
- `katsura9_degrevlex`: Katsura-9 system (9 variables, dense)
- `cyclic_scaling`: Cyclic-5, Cyclic-6, Cyclic-7 for scaling analysis

### 2. F4 with AVX512 IFMA

```bash
# Run F4 benchmarks with AVX512 IFMA enabled
cargo bench --bench f4_avx512 --features avx512,avx512-ifma

# View results
ls -lh target/criterion/f4_*
```

**Test Systems:** (Same as Buchberger for direct comparison)
- `cyclic7_degrevlex_avx512`: Cyclic-7 with AVX512
- `katsura9_degrevlex_avx512`: Katsura-9 with AVX512
- `cyclic_scaling`: Cyclic-5, 6, 7 with AVX512

### 3. Run Both Sequentially

```bash
# Run all benchmarks in sequence
cargo bench --bench buchberger_baseline && \
cargo bench --bench f4_avx512 --features avx512,avx512-ifma
```

## Benchmark Configuration

### Sample Size & Duration
- **Sample size**: 10 iterations (reduced for long-running Gröbner basis computations)
- **Measurement time**: 120-180 seconds per benchmark
- **Warmup**: 3 seconds (Criterion default)

### Polynomial System Configuration
- **Field**: BN254 Fr (256-bit prime field, 4 limbs - optimal for AVX512 IFMA)
- **Degree configuration**:
  - Cyclic/Scaling tests: `DegreeCfg::new(20).with_precompute(10)`
  - Katsura tests: `DegreeCfg::new(15).with_precompute(10)`
  - Max supported degree: 15-20 (depends on test)
  - Precompute degree: 10 (binomial coefficient tables for monomial indexing)
- **Multiplication table**: `(0, 0)` - no precomputed multiplication tables (minimal memory)

### Test Systems Details

#### Cyclic-7
- **Variables**: 7
- **Equations**: 7 (homogeneous cyclic system)
- **Max degree**: ~7
- **Characteristics**: Dense, symmetric

#### Katsura-9
- **Variables**: 9
- **Equations**: 9 (Katsura magnetism problem)
- **Max degree**: ~4-5
- **Characteristics**: Dense quadratic system

#### Scaling Tests
- **Cyclic-5**: 5 variables (baseline)
- **Cyclic-6**: 6 variables (moderate)
- **Cyclic-7**: 7 variables (challenging)

## Interpreting Results

### Criterion Output Location
```
target/criterion/
├── buchberger_cyclic7/
│   └── cyclic7_degrevlex/
│       ├── base/
│       │   ├── estimates.json
│       │   └── sample.json
│       └── report/
│           └── index.html
├── f4_cyclic7/
│   └── cyclic7_degrevlex_avx512/
└── ...
```

### Key Metrics to Compare

1. **Total Time**: Overall Gröbner basis computation time
2. **Throughput**: Iterations per second (if applicable)
3. **Variance**: Consistency of performance

### Expected Performance

Based on ark-algebra AVX512 IFMA documentation:
- **Field multiplication speedup**: ~3.0x (per operation)
- **F4 algorithm benefit**: Batch matrix reduction makes heavy use of field operations
- **Overall speedup**: Depends on proportion of time spent in field arithmetic vs bookkeeping

### Viewing Results

```bash
# Open HTML reports (if gnuplot/plotters available)
firefox target/criterion/buchberger_cyclic7/cyclic7_degrevlex/report/index.html

# Compare raw timing data
cat target/criterion/buchberger_cyclic7/cyclic7_degrevlex/base/estimates.json
cat target/criterion/f4_cyclic7/cyclic7_degrevlex_avx512/base/estimates.json
```

## Troubleshooting

### Benchmark Panics
If benchmarks panic during execution:
```bash
# Run with backtrace
RUST_BACKTRACE=1 cargo bench --bench buchberger_baseline

# Run specific benchmark
cargo bench --bench buchberger_baseline -- cyclic7_degrevlex
```

### Memory Issues
Large polynomial systems may require significant memory:
- Cyclic-7: ~100MB-1GB
- Katsura-9: ~500MB-2GB

If running out of memory, reduce degree limits or use smaller systems.

### AVX512 Not Enabled
Verify AVX512 features are being used:
```bash
# Check if AVX512 instructions are in binary
objdump -d target/release/deps/f4_avx512-* | grep vpmadd52

# Verify CPU supports IFMA
grep avx512_ifma /proc/cpuinfo
```

## Comparison Analysis

### Manual Comparison
Extract timing data and compute speedup:

```bash
# Buchberger Cyclic-7 time
jq '.mean.point_estimate' \
  target/criterion/buchberger_cyclic7/cyclic7_degrevlex/base/estimates.json

# F4 + AVX512 Cyclic-7 time
jq '.mean.point_estimate' \
  target/criterion/f4_cyclic7/cyclic7_degrevlex_avx512/base/estimates.json

# Compute speedup: buchberger_time / f4_avx512_time
```

### Automated Analysis Script (Optional)

Create `compare_benchmarks.sh`:
```bash
#!/bin/bash
BUCH_TIME=$(jq -r '.mean.point_estimate' target/criterion/buchberger_cyclic7/cyclic7_degrevlex/base/estimates.json)
F4_TIME=$(jq -r '.mean.point_estimate' target/criterion/f4_cyclic7/cyclic7_degrevlex_avx512/base/estimates.json)

echo "Buchberger Cyclic-7: ${BUCH_TIME}ns"
echo "F4 + AVX512 Cyclic-7: ${F4_TIME}ns"
echo "Speedup: $(echo "scale=2; $BUCH_TIME / $F4_TIME" | bc)x"
```

## Next Steps

After collecting benchmark data:

1. **Document findings**: Create performance comparison table
2. **Verify correctness**: Both algorithms should produce equivalent Gröbner bases
3. **Profile**: Use `perf` or `cargo flamegraph` to identify hotspots
4. **Optimize**: If speedup is less than expected, investigate AVX512 utilization

## References

- **ark-algebra AVX512 IFMA Guide**: `../ark-algebra/ff/doc/AVX512_GUIDE.md`
- **Criterion.rs Documentation**: https://bheisler.github.io/criterion.rs/
- **F4 Algorithm**: Faugère, J.C. (1999). A new efficient algorithm for computing Gröbner bases (F4)

## Notes

- Benchmarks include polynomial system construction time (unavoidable due to API constraints)
- Clone/copy operations are minimized - systems are reconstructed each iteration
- Results should be reproducible within ~5% variance on the same hardware
- Warm up runs help stabilize CPU frequency/turbo boost
