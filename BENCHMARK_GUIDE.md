# Benchmark Comparison Guide: F4 vs Buchberger

This guide explains how to run and interpret the benchmark suite that compares the Buchberger and F4 algorithms for Gröbner basis computation.

## Available Benchmarks

| Benchmark | Description | Duration |
|-----------|-------------|----------|
| `buchberger_baseline` | Buchberger on Cyclic-7, Katsura-9 | ~30-120s per test |
| `comparison` | Side-by-side F4 vs Buchberger | Varies |
| `buchberger_tiny` | Buchberger on small problems (Cyclic-3/4, Katsura-3/4) | ~5-10 min total |
| `f4_tiny` | F4 on small problems (Cyclic-3/4, Katsura-3/4) | ~5-10 min total |
| `field_ops` | Field operation microbenchmarks | ~1 min |

## Running Benchmarks

### 1. Quick Test: Microbenchmarks (~1 minute)

Tests raw field operations:

```bash
cargo bench --bench field_ops
```

### 2. Fast: Tiny Problems (~10-15 minutes total)

Small Gröbner basis problems (Cyclic-3/4, Katsura-3/4):

```bash
# Buchberger on tiny problems
cargo bench --bench buchberger_tiny

# F4 on tiny problems
cargo bench --bench f4_tiny
```

**Expected**: F4 may have more overhead on very small problems but scales better.

### 3. Full Comparison Suite

```bash
#!/bin/bash
# Run complete benchmark suite

echo "=== 1. Field Operations ==="
cargo bench --bench field_ops

echo "=== 2. Buchberger Tiny ==="
cargo bench --bench buchberger_tiny

echo "=== 3. F4 Tiny ==="
cargo bench --bench f4_tiny

echo "=== 4. Buchberger Baseline ==="
cargo bench --bench buchberger_baseline

echo "=== 5. Comparison ==="
cargo bench --bench comparison

echo "=== Done! Results in target/criterion/ ==="
```

### 4. Large Problems

For larger test systems:

```bash
# Buchberger on Cyclic-7, Katsura-9 (30-120s each)
cargo bench --bench buchberger_baseline

# Direct comparison
cargo bench --bench comparison
```

## Interpreting Results

### Results Location

```
target/criterion/
├── field_mul/
├── buchberger_cyclic3/
├── f4_cyclic3/
├── buchberger_cyclic7/
├── f4_cyclic7/
└── ...
```

### Extracting Timing Data

```bash
# Get mean time for a benchmark
jq '.mean.point_estimate' target/criterion/<group>/<test>/base/estimates.json

# Example: Buchberger Cyclic-3 (convert ns to ms)
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
| Cyclic-7  | 45000 | 15000 | 3.00x | F4 |
| Katsura-3 | 80 | 100 | 0.80x | Buchberger |
| Katsura-4 | 300 | 250 | 1.20x | F4 |

### Key Questions to Answer

1. **Algorithm Performance**:
   - Is F4 faster than Buchberger on these systems?
   - At what problem size does F4 overtake Buchberger?

2. **Current Overhead**:
   - Does F4 have more overhead on tiny problems?
   - Where's the crossover point?

3. **Scaling Behavior**:
   - How does each algorithm scale with system size?
   - Which is more predictable?

## Expected Findings

### F4 vs Buchberger

**Tiny problems (3-4 vars)**: Buchberger might be faster
- F4 has matrix setup overhead
- Not enough operations to amortize cost
- Buchberger's simplicity wins on small inputs

**Medium problems (5-7 vars)**: F4 starts to win
- Batch reduction becomes efficient
- Amortizes setup cost
- Matrix approach pays off

**Larger problems (7+ vars)**: F4 should be significantly faster
- F4 designed for larger systems
- Can handle problems that would timeout with Buchberger
- Better asymptotic complexity

## Next Steps Based on Results

### If F4 Is Slower Than Buchberger on All Tests

F4 overhead doesn't pay off on these small systems. Options:

1. Test larger systems (Cyclic-8, Cyclic-9, Katsura-10)
2. Optimize F4 implementation (reduce allocations, better data structures)
3. Use Buchberger for small problems, F4 for large ones
4. Profile to find bottlenecks

### If F4 Is Faster on Larger Tests

F4 is the right choice! Next steps:

1. Find the crossover point (what size makes F4 worthwhile?)
2. Consider hybrid: Buchberger for small, F4 for large
3. Profile to find bottlenecks: `cargo flamegraph --bench f4_tiny`
4. Optimize the bottlenecks identified

## Verification Checklist

Before drawing conclusions, verify:

- [ ] Compiled in release mode with LTO
- [ ] Benchmark results are consistent across runs (< 5% variance)
- [ ] Sample sizes are adequate (criterion automatic)
- [ ] No thermal throttling during long benchmarks
- [ ] Same polynomial systems used for both algorithms

## Advanced: Profiling

To see where time is actually spent:

```bash
# Install flamegraph
cargo install flamegraph

# Profile F4 tiny (requires root or CAP_PERFMON)
sudo cargo flamegraph --bench f4_tiny -- --bench cyclic4

# View flamegraph.svg in browser
firefox flamegraph.svg
```

Look for:
- How much time in field operations?
- How much in matrix construction?
- How much in memory allocation?
- Where are the bottlenecks?

You can also use `perf` for more detailed analysis:

```bash
# Record performance data
perf record --call-graph dwarf cargo bench --bench f4_tiny -- cyclic4

# Generate report
perf report

# Generate flamegraph-compatible data
perf script | stackcollapse-perf.pl | flamegraph.pl > flamegraph.svg
```

## Optimization Opportunities

Common bottlenecks to investigate:

1. **Field arithmetic**: Can we reduce division operations?
2. **Memory allocation**: Are we allocating/deallocating too often?
3. **Cache efficiency**: Is data access pattern cache-friendly?
4. **Monomial operations**: Can divmask/signature caching help?
5. **Matrix sparsity**: Are we exploiting sparsity effectively?

## Summary

This benchmark suite provides:

1. **Algorithm Comparison**: Clear comparison of Buchberger vs F4
2. **Microbenchmarks**: Isolate field operation performance
3. **Tiny Problems**: Fast iteration for development
4. **Large Tests**: Full algorithm validation and scaling behavior

Key insight: F4 is designed for larger problems where its matrix-based approach
can amortize setup costs. For very small problems, simpler algorithms like
Buchberger may be more efficient.

## Troubleshooting

### Benchmark Won't Compile

```bash
# Verbose output
cargo bench --bench f4_tiny --verbose

# Check dependencies
cargo tree

# Clean and rebuild
cargo clean && cargo build --release
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
- Ensure sufficient memory available

### Out of Memory

Large systems may exhaust memory:
- Monitor with `htop` during benchmarks
- Reduce problem size in benchmark file
- Increase swap space if needed
- Consider testing on machine with more RAM

## References

- [Criterion.rs Guide](https://bheisler.github.io/criterion.rs/book/)
- [F4 Algorithm Paper](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf) - Faugère, J.C. (1999)
- [Gröbner Bases Book](https://www.springer.com/gp/book/9780387979717) - Cox, Little, O'Shea
