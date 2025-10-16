//! # ark-feanor
//! 
//! A bridge between arkworks' finite field types and feanor-math's ring system.
//! 
//! This library enables using feanor-math's advanced polynomial and Gröbner basis
//! algorithms with arkworks' cryptographic field implementations, particularly for
//! curves like BLS12-381 and BN254.
//! 
//! ## Features
//! 
//! - **Zero-cost abstraction**: Minimal overhead when wrapping arkworks fields
//! - **Full ring support**: Implements all necessary feanor-math traits
//! - **Prime field specialization**: Additional functionality for prime fields including division
//! - **Common fields**: Pre-configured types for BN254 and BLS12-381
//! - **Polynomial support**: Create and manipulate polynomials over cryptographic fields
//! - **Gröbner basis**: Compute Gröbner bases for polynomial systems
//! 
//! ## Quick Start
//! 
//! ```ignore
//! use ark_feanor::*;
//! use ark_bn254::Fr;
//! 
//! // Method 1: Use pre-configured field
//! let field = &*BN254_FR;
//! let a = field.int_hom().map(5);
//! let b = field.int_hom().map(3);
//! let c = field.add(a, b);
//! 
//! // Method 2: Create your own wrapper
//! let my_field = ArkFieldWrapper::<Fr>::new();
//! let x = my_field.from_int(10);
//! let y = my_field.from_int(20);
//! let z = my_field.mul_ref(&x, &y);
//! ```
//! 
//! ## Polynomial Example
//! 
//! ```ignore
//! use ark_feanor::*;
//! use feanor_math::rings::multivariate::*;
//! 
//! // Create polynomial ring BLS12-381_Fr[x, y]
//! let field = &*BLS12_381_FR;
//! let poly_ring = MultivariatePolyRingImpl::new(field, 2);
//! 
//! // Create polynomials: x² + y² - 1 and xy - 2
//! let polys = poly_ring.with_wrapped_indeterminates(|[x, y]| {
//!     vec![
//!         x.clone().pow(2) + y.clone().pow(2) - 1,
//!         x * y - 2,
//!     ]
//! });
//! ```
//! 
//! ## Gröbner Basis Computation
//!
//! ### Basic Usage
//!
//! ```ignore
//! use ark_feanor::*;
//!
//! let field = &*BN254_FR;
//! let poly_ring = MultivariatePolyRingImpl::new(field, 2);
//!
//! // Define your polynomial system
//! let system = vec![/* your polynomials */];
//!
//! // Compute Gröbner basis with degree reverse lexicographic ordering
//! let gb = buchberger_simple(&poly_ring, system, DegRevLex);
//! ```
//!
//! ### With Degree Limiting (Recommended for Large Systems)
//!
//! ```ignore
//! use ark_feanor::*;
//!
//! let field = &*BN254_FR;
//! let poly_ring = MultivariatePolyRingImpl::new(field, 64);
//!
//! // Configure with degree limit to avoid explosion
//! let config = BuchbergerConfig::new()
//!     .with_max_degree(10);  // Abort if degrees exceed 10
//!
//! let system = vec![/* your polynomials */];
//!
//! // Compute with configuration
//! match buchberger_configured(&poly_ring, system, Lex, config) {
//!     Ok(gb) => println!("Computed GB with {} elements", gb.len()),
//!     Err(GBAborted::DegreeExceeded { max_degree, actual_degree }) => {
//!         println!("Aborted: degree {} exceeded limit {}", actual_degree, max_degree);
//!     }
//!     Err(e) => println!("Error: {}", e),
//! }
//! ```
//!
//! ### Variable Elimination with BlockLex
//!
//! ```ignore
//! use ark_feanor::*;
//!
//! let field = &*BN254_FR;
//! let poly_ring = MultivariatePolyRingImpl::new(field, 65);
//!
//! // Eliminate first 64 variables, keep last variable
//! let order = BlockLex::new(64);
//! let config = BuchbergerConfig::new().with_max_degree(5);
//!
//! let system = vec![/* your polynomials */];
//!
//! // Compute elimination GB
//! let gb = buchberger_configured(&poly_ring, system, order, config)?;
//!
//! // Polynomials involving only variable 64 form the elimination ideal
//! ```
//! 
//! ## Type Aliases
//! 
//! For convenience, this library provides type aliases for common fields:
//! 
//! - `BN254ScalarField` / `BN254_FR`: The BN254 scalar field (Fr)
//! - `BN254BaseField` / `BN254_FQ`: The BN254 base field (Fq)
//! - `BLS12_381ScalarField` / `BLS12_381_FR`: The BLS12-381 scalar field (Fr)
//! - `BLS12_381BaseField` / `BLS12_381_FQ`: The BLS12-381 base field (Fq)
//! 
//! ## Performance
//! 
//! The wrapper is designed to have minimal overhead. Field operations through
//! `ArkFieldWrapper` should be within 5% of direct arkworks operations.
//! 
//! ## Limitations
//! 
//! - Extension fields (Fq2, Fq6, Fq12) are not yet optimally supported
//! - Characteristic extraction only works for `PrimeField` types
//! - Some advanced feanor-math features may not be available for all field types

#![doc(html_root_url = "https://docs.rs/ark-feanor/0.1.0")]
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

// Core modules
pub mod field_wrapper;
pub mod prime_field;
pub mod homomorphisms;
pub mod conversions;
pub mod common_fields;

// Re-export main types
pub use field_wrapper::ArkFieldWrapper;
pub use common_fields::*;
pub use conversions::{
    extract_characteristic,
    biguint_to_field,
    bigint_to_field,
    i64_to_field,
    field_to_biguint,
    field_to_small_int,
};
pub use prime_field::FieldProperties;

// Re-export useful feanor-math types for convenience
pub use feanor_math::{
    ring::*,
    divisibility::*,
    pid::*,
    rings::{
        multivariate::{
            multivariate_impl::MultivariatePolyRingImpl,
            MultivariatePolyRingStore,
            DegRevLex,
            Lex,
            BlockLex,
        },
        poly::{
            PolyRingStore,
            dense_poly::DensePolyRing,
        },
    },
    algorithms::buchberger::{
        buchberger_simple,
        buchberger_configured,
        buchberger_with_sugar,
        BuchbergerConfig,
        GBAborted,
    },
};

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Library name
pub const NAME: &str = env!("CARGO_PKG_NAME");

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    
    #[test]
    fn test_library_basics() {
        // Ensure the library exports work
        let field = ArkFieldWrapper::<Fr>::new();
        let a = field.from_int(5);
        let b = field.from_int(10);
        let c = field.add_ref(&a, &b);
        let expected = field.from_int(15);
        assert!(field.eq_el(&c, &expected));
    }
    
    #[test]
    fn test_version_info() {
        assert_eq!(NAME, "ark-feanor");
        assert!(!VERSION.is_empty());
    }
}
