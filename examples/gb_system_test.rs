#![feature(allocator_api)]

//! Standalone GB test using serialized polynomial systems
//!
//! This tool tests Gröbner basis algorithms (Buchberger vs F4) on real zyga systems
//! without needing to compile circuits or run full setup.
//!
//! Usage:
//!   cargo run --release --features elimination --bin gb_system_test -- <system.json>
//!
//! Example:
//!   cargo run --release --features elimination --bin gb_system_test -- uint4_system.json
//!
//! Options:
//!   --skip-buchberger    Skip Buchberger (for large systems)
//!   --order <Lex|DegRevLex>  Change monomial ordering (default: Lex)

use std::env;
use std::fs;
use std::time::Instant;
use std::alloc::Global;
use serde::{Serialize, Deserialize};
use ark_bn254::Fr;

// Import everything from ark-feanor (includes re-exports from feanor-math)
use ark_feanor::*;
use ark_feanor::f4::f4_simple;

/// Re-export the serializable types from gb_diagnostic
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableSystem {
    pub ring_config: RingConfig,
    pub variables: VariableInfo,
    pub system: R1CSSystem,
    pub gb_config: GBConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RingConfig {
    pub n_vars: usize,
    pub max_degree: usize,
    pub mult_table: (usize, usize),
    pub field: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariableInfo {
    pub names: Vec<String>,
    pub roles: Vec<String>,
    pub public_names: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R1CSSystem {
    pub n_constraints: usize,
    pub a_matrix: Vec<Vec<f64>>,
    pub b_matrix: Vec<Vec<f64>>,
    pub c_matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GBConfig {
    pub max_degree: usize,
    pub s_pair_budget: usize,
    pub monomial_order: String,
    pub parallel: bool,
}

// BN254_FR is imported from ark_feanor (already a static reference)

fn main() {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();

    if args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        println!("Usage: {} [system.json] [--skip-buchberger] [--order <Lex|DegRevLex>]", args[0]);
        println!("\nArguments:");
        println!("  system.json          Path to serialized system (default: examples/uint4_system.json)");
        println!("                       Use examples/amm_system.json for large AMM system");
        println!("\nOptions:");
        println!("  --skip-buchberger    Skip Buchberger test (for large systems)");
        println!("  --order <Lex|DegRevLex>  Monomial ordering (default: Lex)");
        println!("\nExamples:");
        println!("  {}                   # Use default uint4 system", args[0]);
        println!("  {} examples/amm_system.json --skip-buchberger", args[0]);
        println!("  {} examples/uint4_system.json --order DegRevLex", args[0]);
        std::process::exit(0);
    }

    // Determine system path: use first non-flag argument, or default to uint4
    let system_path = args.iter()
        .skip(1)
        .find(|arg| !arg.starts_with("--"))
        .map(|s| s.as_str())
        .unwrap_or("examples/uint4_system.json");

    let skip_buchberger = args.contains(&"--skip-buchberger".to_string());
    let use_degrevlex = args.iter().any(|a| a == "--order")
        && args.iter().position(|a| a == "--order")
            .and_then(|i| args.get(i + 1))
            .map(|s| s == "DegRevLex")
            .unwrap_or(false);

    println!("=== GB System Test ===\n");
    println!("System file: {}", system_path);
    println!("Skip Buchberger: {}", skip_buchberger);
    println!("Monomial order: {}\n", if use_degrevlex { "DegRevLex" } else { "Lex" });

    // Step 1: Load serialized system
    println!("[1/5] Loading serialized system...");
    let start = Instant::now();
    let json = fs::read_to_string(system_path)
        .expect("Failed to read system file");
    let system: SerializableSystem = serde_json::from_str(&json)
        .expect("Failed to parse JSON");
    println!("  Loaded in {:?}", start.elapsed());
    println!("  Variables: {}", system.ring_config.n_vars);
    println!("  Constraints: {}", system.system.n_constraints);
    println!("  Deferred vars: {}", system.variables.roles.iter().filter(|r| *r == "Deferred").count());
    println!();

    // Step 2: Create polynomial ring
    println!("[2/5] Creating polynomial ring...");
    let start = Instant::now();
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new_with_mult_table(
        field,
        system.ring_config.n_vars,
        system.ring_config.max_degree as u16,
        (system.ring_config.mult_table.0 as u16, system.ring_config.mult_table.1 as u16),
        Global,
    );
    println!("  Created in {:?}", start.elapsed());
    println!();

    // Step 3: Build polynomial system from R1CS matrices
    println!("[3/5] Building polynomial system from R1CS...");
    let start = Instant::now();

    // Helper to build linear polynomial from coefficient vector
    let build_linear = |coeffs: &[f64]| {
        let mut terms: Vec<(_, _)> = Vec::new();
        for (var_idx, &c) in coeffs.iter().enumerate() {
            let k = c as i64;
            if k == 0 { continue; }  // Skip zero coefficients

            // Create exponent vector for this variable
            let mut exp = vec![0usize; system.ring_config.n_vars];
            exp[var_idx] = 1;

            let mon = poly_ring.create_monomial(exp.into_iter());
            let coeff = i64_to_field::<Fr>(k);
            terms.push((coeff, mon));
        }
        poly_ring.from_terms(terms.into_iter())
    };

    // Build system: each constraint becomes A(X)*B(X) - C(X)
    let mut poly_system = Vec::with_capacity(system.system.n_constraints);
    for i in 0..system.system.n_constraints {
        let a_lin = build_linear(&system.system.a_matrix[i]);
        let b_lin = build_linear(&system.system.b_matrix[i]);
        let c_lin = build_linear(&system.system.c_matrix[i]);

        let ab = poly_ring.mul(a_lin, b_lin);       // A(X) * B(X)
        let neg_c = poly_ring.negate(c_lin);        // -C(X)
        let f = poly_ring.add(ab, neg_c);           // A*B - C

        poly_system.push(f);
    }

    println!("  Built {} polynomials in {:?}", poly_system.len(), start.elapsed());
    println!();

    // Step 4: Run Buchberger (if not skipped)
    let buchberger_time = if !skip_buchberger {
        println!("[4/5] Running Buchberger algorithm...");

        // Need to rebuild system since polynomials don't implement Clone
        let mut poly_system_buch = Vec::with_capacity(system.system.n_constraints);
        for i in 0..system.system.n_constraints {
            let a_lin = build_linear(&system.system.a_matrix[i]);
            let b_lin = build_linear(&system.system.b_matrix[i]);
            let c_lin = build_linear(&system.system.c_matrix[i]);
            let ab = poly_ring.mul(a_lin, b_lin);
            let neg_c = poly_ring.negate(c_lin);
            let f = poly_ring.add(ab, neg_c);
            poly_system_buch.push(f);
        }

        let start = Instant::now();
        let buchberger_result = if use_degrevlex {
            buchberger_simple(&poly_ring, poly_system_buch, DegRevLex)
        } else {
            buchberger_simple(&poly_ring, poly_system_buch, Lex)
        };
        let elapsed = start.elapsed();

        println!("  Completed in {:?}", elapsed);
        println!("  Basis size: {} polynomials", buchberger_result.len());
        println!();
        Some(elapsed)
    } else {
        println!("[4/5] Buchberger skipped");
        println!();
        None
    };

    // Step 5: Run F4 algorithm
    println!("[5/5] Running F4 algorithm...");
    let start = Instant::now();

    let f4_result = if use_degrevlex {
        f4_simple(&poly_ring, poly_system, DegRevLex)
    } else {
        f4_simple(&poly_ring, poly_system, Lex)
    };

    let f4_time = start.elapsed();
    println!("  Completed in {:?}", f4_time);
    println!("  Basis size: {} polynomials", f4_result.len());
    println!();

    // Summary
    println!("=== Summary ===");
    println!("System: {} vars, {} constraints", system.ring_config.n_vars, system.system.n_constraints);
    if let Some(buch_time) = buchberger_time {
        println!("Buchberger: {:?}", buch_time);
        println!("F4:         {:?}", f4_time);
        let speedup = buch_time.as_secs_f64() / f4_time.as_secs_f64();
        println!("Speedup:    {:.2}x", speedup);
    } else {
        println!("F4:         {:?}", f4_time);
    }
    println!("\n=== Test Complete ===");
}
