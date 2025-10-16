# ark-feanor

A bridge between [arkworks](https://github.com/arkworks-rs)' finite field types and [feanor-math](https://github.com/FeanorTheElf/feanor-math)'s ring system.

[![Crates.io](https://img.shields.io/crates/v/ark-feanor.svg)](https://crates.io/crates/ark-feanor)
[![Documentation](https://docs.rs/ark-feanor/badge.svg)](https://docs.rs/ark-feanor)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)

## Overview

`ark-feanor` enables you to use feanor-math's advanced polynomial and Gröbner basis algorithms with arkworks' cryptographic field implementations. This is particularly useful for:

- Computing Gröbner bases over fields from curves like BLS12-381 and BN254
- Solving polynomial systems in cryptographic contexts
- Performing symbolic computation over finite fields
- Zero-knowledge proof system development

## Features

- **Zero-cost abstraction**: Minimal overhead when wrapping arkworks fields
- **Full ring support**: Complete implementation of feanor-math's ring traits
- **Prime field specialization**: Division and field operations for `PrimeField` types
- **Pre-configured fields**: Ready-to-use wrappers for BN254 and BLS12-381
- **Polynomial support**: Create and manipulate univariate and multivariate polynomials
- **Gröbner basis computation**: Solve polynomial systems using Buchberger's algorithm

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
ark-feanor = "0.1"
```

## Quick Start

### Basic Field Operations

```rust
use ark_feanor::*;

// Use pre-configured BN254 scalar field
let field = &*BN254_FR;
let a = field.int_hom().map(5);
let b = field.int_hom().map(3);
let c = field.add(a, b);
assert_eq!(field.int_hom().map(8), c);

// Or create your own wrapper
use ark_bn254::Fr;
let my_field = ArkFieldWrapper::<Fr>::new();
let x = my_field.from_int(10);
let y = my_field.from_int(20);
let z = my_field.mul_ref(&x, &y);
```

### Polynomial Rings

```rust
use ark_feanor::*;
use feanor_math::rings::poly::dense_poly::DensePolyRing;
use feanor_math::rings::poly::PolyRingStore;

// Create univariate polynomial ring BLS12-381_Fr[x]
let field = &*BLS12_381_FR;
let poly_ring = DensePolyRing::new(field.clone(), "x");

// Create polynomial: 3x² + 2x + 1
let poly = poly_ring.from_terms([
    (field.int_hom().map(1), 0),  // 1
    (field.int_hom().map(2), 1),  // 2x
    (field.int_hom().map(3), 2),  // 3x²
].iter().cloned());

// Evaluate at x = 5
let result = poly_ring.evaluate(&poly, &field.int_hom().map(5), &field.identity());
```

### Gröbner Basis Computation

```rust
use ark_feanor::*;
use feanor_math::algorithms::buchberger::*;
use feanor_math::rings::multivariate::*;

// Create multivariate ring BN254_Fr[x, y]
let field = &*BN254_FR;
let poly_ring = MultivariatePolyRingImpl::new(field, 2);

// Define polynomial system:
// x² + y² - 1 = 0
// xy - 2 = 0
let system = poly_ring.with_wrapped_indeterminates(|[x, y]| {
    vec![
        x.clone().pow(2) + y.clone().pow(2) - 1,
        x.clone() * y.clone() - 2,
    ]
});

// Compute Gröbner basis
let gb = buchberger(&poly_ring, system, DegRevLex, |_| {})?;
println!("Basis has {} polynomials", gb.len());
```

## Supported Fields

### Pre-configured Fields

The library provides convenient access to commonly used cryptographic fields:

- **BN254**:
  - `BN254_FR` / `BN254ScalarField`: Scalar field (Fr)
  - `BN254_FQ` / `BN254BaseField`: Base field (Fq)

- **BLS12-381**:
  - `BLS12_381_FR` / `BLS12_381ScalarField`: Scalar field (Fr)
  - `BLS12_381_FQ` / `BLS12_381BaseField`: Base field (Fq)

### Custom Fields

Any arkworks field implementing `Field` or `PrimeField` can be wrapped:

```rust
use ark_feanor::ArkFieldWrapper;
use ark_bw6_761::Fr;

let field = ArkFieldWrapper::<Fr>::new();
```

## Examples

Check out the `examples/` directory for complete working examples:

- `simple_groebner.rs`: Basic Gröbner basis computation
- `polynomial_system.rs`: Solving polynomial systems

Run examples with:

```bash
cargo run --example simple_groebner
```

## Performance

The wrapper is designed for minimal overhead. Benchmarks show:

- Field arithmetic operations: < 1% overhead
- Polynomial operations: < 5% overhead
- Gröbner basis computation: Comparable to native implementations

Run benchmarks with:

```bash
cargo bench
```

## Architecture

The library is organized into several modules:

- `field_wrapper`: Core `RingBase` trait implementation
- `prime_field`: Specialized traits for prime fields (division, etc.)
- `conversions`: Type conversions between arkworks and feanor-math
- `homomorphisms`: Ring homomorphisms and embeddings
- `common_fields`: Pre-configured field types

## Limitations

- Extension fields (Fq2, Fq6, Fq12) are not yet optimally supported
- Characteristic extraction only works for `PrimeField` types
- Some advanced feanor-math features may not be available

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

### Development

```bash
# Run tests
cargo test

# Run with all features
cargo test --all-features

# Check documentation
cargo doc --open

# Run clippy
cargo clippy -- -D warnings
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Acknowledgments

- [arkworks](https://github.com/arkworks-rs) for the excellent finite field implementations
- [feanor-math](https://github.com/FeanorTheElf/feanor-math) for the comprehensive computer algebra system
- The Rust cryptography community

## Related Projects

- [arkworks-algebra](https://github.com/arkworks-rs/algebra): Algebraic structures for cryptography
- [feanor-math](https://github.com/FeanorTheElf/feanor-math): Computer algebra system in Rust
- [ark-poly](https://github.com/arkworks-rs/poly): Polynomial arithmetic over finite fields
