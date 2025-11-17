# Benchmark Comparison Guide: F4 vs Buchberger & AVX512 Impact

This guide explains how to run and interpret the benchmark suite that tests:
1. **Algorithm comparison**: Buchberger vs F4
2. **AVX512 impact**: With vs without AVX512-IFMA optimizations

## Critical Discovery: AVX512 Is Not Automatically Used

**IMPORTANT**: The AVX512 batch operations in `ark-algebra` require explicit function calls
to `avx512_backend::mont_mul_batch_8()`. The standard field operations used by Buchberger
and F4 do NOT automatically use AVX512 instructions.

Current status:
- ❌ Buchberger: Does not use AVX512 batch ops
- ❌ F4: Does not use AVX512 batch ops
- ✅ Field microbenchmarks: Directly test AVX512 batch operations

To actually benefit from AVX512, you would need to integrate the batch operations into
the matrix reduction code or polynomial arithmetic.

## Benchmark Matrix (2×2)

|  | **No AVX512** | **With AVX512 Feature** |
|--|---------------|-------------------------|
| **Buchberger** | `buchberger_tiny` | `buchberger_tiny_avx512` |
| **F4** | `f4_tiny` | `f4_tiny_avx512` |

Plus:
- **Field ops**: `field_ops` vs `field_ops_avx512` (microbenchmarks)

## Running Benchmarks

### 1. Quick Test: Microbenchmarks Only (~1 minute)

Tests raw field operations to see AVX512's theoretical potential:

```bash
# Without AVX512
cargo bench --bench field_ops

# With AVX512
cargo bench --bench field_ops_avx512 --features avx512-ifma
```

**Expected**: AVX512 batch operations should be ~2-3x faster than sequential operations.

### 2. Fast: Tiny Problems (~10-15 minutes total)

Small Gröbner basis problems (Cyclic-3/4, Katsura-3/4):

```bash
# Buchberger without AVX512
cargo bench --bench buchberger_tiny

# Buchberger with AVX512 (currently same as above)
cargo bench --bench buchberger_tiny_avx512 --features avx512-ifma

# F4 without AVX512
cargo bench --bench f4_tiny

# F4 with AVX512 (currently same as above)
cargo bench --bench f4_tiny_avx512 --features avx512-ifma
```

**Expected**: F4 and Buchberger should have similar performance since AVX512
isn't actually being used. This shows pure algorithm differences.

### 3. Full Suite: Run All 2×2 Combinations

```bash
#!/bin/bash
# Run complete 2x2 matrix

echo "=== 1. Field Operations (No AVX512) ==="
cargo bench --bench field_ops

echo "=== 2. Field Operations (AVX512) ==="
cargo bench --bench field_ops_avx512 --features avx512-ifma

echo "=== 3. Buchberger Tiny (No AVX512) ==="
cargo bench --bench buchberger_tiny

echo "=== 4. Buchberger Tiny (AVX512 feature) ==="
cargo bench --bench buchberger_tiny_avx512 --features avx512-ifma

echo "=== 5. F4 Tiny (No AVX512) ==="
cargo bench --bench f4_tiny

echo "=== 6. F4 Tiny (AVX512 feature) ==="
cargo bench --bench f4_tiny_avx512 --features avx512-ifma

echo "=== Done! Results in target/criterion/ ==="
```

### 4. Large Problems (Original Benchmarks)

For completeness, your original benchmarks with larger systems:

```bash
# Buchberger on Cyclic-7, Katsura-7 (30-120s each)
cargo bench --bench buchberger_baseline

# F4 on Cyclic-7, Katsura-7 (currently running)
cargo bench --bench f4_avx512 --features avx512-ifma
```

## Interpreting Results

### Results Location

```
target/criterion/
├── field_mul/
│   └── bn254_mul_1M/
├── field_mul_avx512/
│   └── bn254_mul_8elem_100K_avx512_batch/
├── buchberger_cyclic3/
├── f4_cyclic3/
└── ...
```

### Extracting Timing Data

```bash
# Get mean time for a benchmark
jq '.mean.point_estimate' target/criterion/<group>/<test>/base/estimates.json

# Example: Buchberger Cyclic-3
jq '.mean.point_estimate / 1000000' \
  target/criterion/buchberger_cyclic3/cyclic3_degrevlex/base/estimates.json
```

Output is in nanoseconds. Divide by 1,000,000 for milliseconds.

### Comparison Template

Create a comparison table:

| Benchmark | Buchberger (ms) | F4 (ms) | Speedup | Winner |
|-----------|-----------------|---------|---------|--------|
| Cyclic-3  | 150 | 200 | 0.75x | Buchberger |
| Cyclic-4  | 800 | 600 | 1.33x | F4 |
| Katsura-3 | 80 | 100 | 0.80x | Buchberger |
| Katsura-4 | 300 | 250 | 1.20x | F4 |

| Field Ops | Sequential (ns) | AVX512 Batch (ns) | Speedup |
|-----------|-----------------|-------------------|---------|
| Mul 8×100K | 500M | 180M | 2.78x |

### Key Questions to Answer

1. **Algorithm Performance**:
   - Is F4 faster than Buchberger on these systems?
   - At what problem size does F4 overtake Buchberger?

2. **AVX512 Theoretical Potential**:
   - How much faster are AVX512 batch ops in microbenchmarks?
   - Is it worth integrating AVX512 into the algorithms?

3. **Current Overhead**:
   - Does F4 have more overhead on tiny problems?
   - Where's the crossover point?

## Expected Findings

### Field Microbenchmarks

✅ **AVX512 batch operations should be 2-3x faster** than sequential operations.

If this doesn't show up, check:
- `objdump -d target/release/deps/field_ops_avx512-* | grep vpmadd52`
- Should see `vpmadd52luq` and `vpmadd52huq` instructions

### Algorithm Benchmarks

❌ **Buchberger and F4 "AVX512" versions likely same speed** as non-AVX512 versions.

Why? Because neither algorithm currently calls the AVX512 batch operations. They just
use standard field arithmetic which compiles to the same code regardless of the feature flag.

### F4 vs Buchberger

**Tiny problems (3-4 vars)**: Buchberger might be faster
- F4 has matrix setup overhead
- Not enough operations to amortize cost

**Larger problems (7+ vars)**: F4 should be faster
- Batch reduction becomes efficient
- Amortizes setup cost

## Next Steps Based on Results

### If Field Ops Show Good AVX512 Speedup (2-3x)

AVX512 is working! To actually benefit:

1. **Option A**: Modify F4 matrix reduction to use `avx512_backend::mont_mul_batch_8()`
2. **Option B**: Create AVX512-aware polynomial multiplication
3. **Option C**: Batch field operations in critical loops

This requires significant integration work.

### If F4 Is Slower Than Buchberger on All Tests

F4 overhead doesn't pay off on these small systems. Options:

1. Test larger systems (Cyclic-8, Cyclic-9)
2. Optimize F4 implementation (reduce allocations, better data structures)
3. Stick with Buchberger for your workload

### If F4 Is Faster on Larger Tests

F4 is the right choice! Next steps:

1. Find the crossover point (what size makes F4 worthwhile?)
2. Consider hybrid: Buchberger for small, F4 for large
3. Profile to find bottlenecks: `cargo flamegraph --bench f4_tiny`

## Verification Checklist

Before drawing conclusions, verify:

- [ ] AVX512 instructions appear in `field_ops_avx512` binary
- [ ] CPU supports AVX512-IFMA: `grep avx512_ifma /proc/cpuinfo`
- [ ] Compiled with correct flags: `target-cpu=native` in `.cargo/config.toml`
- [ ] Benchmark results are consistent across runs (< 5% variance)
- [ ] Sample sizes are adequate (criterion automatic)

## Advanced: Profiling

To see where time is actually spent:

```bash
# Install flamegraph
cargo install flamegraph

# Profile F4 tiny (requires root or CAP_PERFMON)
sudo cargo flamegraph --bench f4_tiny --features avx512-ifma -- --bench cyclic4

# View flamegraph.svg in browser
firefox flamegraph.svg
```

Look for:
- How much time in field operations?
- How much in matrix construction?
- How much in memory allocation?

## Summary

This benchmark suite provides:

1. **2×2 Algorithm Matrix**: Clear comparison of Buchberger vs F4 with/without AVX512 feature
2. **Microbenchmarks**: Show AVX512's theoretical potential
3. **Tiny Problems**: Fast iteration for development
4. **Existing Large Tests**: Full algorithm validation

Key insight: AVX512 feature is currently *not used* by the algorithms, only explicitly
in microbenchmarks. Real AVX512 integration would require code changes.

## Troubleshooting

### Benchmark Won't Compile

```bash
# Missing features?
cargo bench --bench f4_tiny_avx512 --features avx512-ifma --verbose

# Check dependencies
cargo tree
```

### Benchmark Takes Too Long

Reduce sample size or measurement time:
- Edit benchmark file
- Change `group.sample_size(50)` to `group.sample_size(10)`
- Change `Duration::from_secs(10)` to `Duration::from_secs(5)`

### Inconsistent Results

- Close other applications
- Disable CPU frequency scaling: `sudo cpupower frequency-set -g performance`
- Run multiple times and average
- Check thermal throttling: `sensors` (install lm-sensors)

## References

- [Criterion.rs Guide](https://bheisler.github.io/criterion.rs/book/)
- [F4 Algorithm Paper](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf)
- [AVX-512 IFMA](https://en.wikipedia.org/wiki/AVX-512#IFMA)
- [ark-algebra AVX512 Guide](../ark-algebra/ff/doc/AVX512_GUIDE.md) (if exists)
