//! Generate uint16 polynomial system for 16-bit addition verification
//!
//! This generates a constraint system for verifying 16-bit unsigned integer addition:
//! a + b = sum (mod 2^16) with carry out
//!
//! Each bit position generates constraints for:
//! - Boolean constraints: x * (1 - x) = 0 for all bits
//! - Sum bit: sum_i = a_i XOR b_i XOR carry_{i-1}
//! - Carry bit: carry_i = MAJ(a_i, b_i, carry_{i-1})
//!
//! Usage: cargo run --example gen_uint16 > examples/uint16_system.json

use serde::{Serialize, Deserialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize)]
struct RingConfig {
    n_vars: usize,
    max_degree: usize,
    mult_table: Vec<usize>,
    field: String,
}

#[derive(Serialize, Deserialize)]
struct Variables {
    names: Vec<String>,
    roles: Vec<String>,
    public_names: Vec<String>,
}

#[derive(Serialize, Deserialize)]
struct System {
    n_constraints: usize,
    polynomials: Vec<String>,
}

#[derive(Serialize, Deserialize)]
struct GBConfig {
    max_degree: usize,
    s_pair_budget: usize,
    monomial_order: String,
    parallel: bool,
}

#[derive(Serialize, Deserialize)]
struct PolySystem {
    ring_config: RingConfig,
    variables: Variables,
    system: System,
    gb_config: GBConfig,
}

fn main() {
    let bits = 16;

    // Build variable list
    let mut var_names = vec!["1".to_string()];
    let mut var_roles = vec!["Deferred".to_string()];
    let mut var_map: HashMap<String, usize> = HashMap::new();
    var_map.insert("1".to_string(), 0);

    // Input bits a[0..15]
    for i in 0..bits {
        let name = format!("a[{}]", i);
        var_map.insert(name.clone(), var_names.len());
        var_names.push(name);
        var_roles.push("Private".to_string());
    }

    // Input bits b[0..15]
    for i in 0..bits {
        let name = format!("b[{}]", i);
        var_map.insert(name.clone(), var_names.len());
        var_names.push(name);
        var_roles.push("Deferred".to_string());  // Public outputs
    }

    // Sum bits sum[0..15]
    for i in 0..bits {
        let name = format!("sum[{}]", i);
        var_map.insert(name.clone(), var_names.len());
        var_names.push(name);
        var_roles.push("Private".to_string());
    }

    // Carry bits carry[0..15] (carry[0] is carry-in, usually 0)
    for i in 0..bits {
        let name = format!("carry[{}]", i);
        var_map.insert(name.clone(), var_names.len());
        var_names.push(name);
        var_roles.push("Private".to_string());
    }

    // Intermediate symbols for XOR computation
    // XOR(a,b) = a + b - 2*a*b
    // For 3-input XOR (a XOR b XOR c), we need intermediate
    for i in 0..bits {
        let name = format!("xor_ab[{}]", i);  // a[i] XOR b[i]
        var_map.insert(name.clone(), var_names.len());
        var_names.push(name);
        var_roles.push("Private".to_string());
    }

    // Carry out
    let name = "carry_out".to_string();
    var_map.insert(name.clone(), var_names.len());
    var_names.push(name);
    var_roles.push("Private".to_string());

    let n_vars = var_names.len();

    // Generate polynomial constraints
    let mut polynomials = Vec::new();

    // 1. Boolean constraints for input bits a[i]: a[i] * (1 - a[i]) = 0
    for i in 0..bits {
        // a[i]^2 - a[i] = 0
        polynomials.push(format!("a[{}]^2 - a[{}]", i, i));
    }

    // 2. Boolean constraints for input bits b[i]: b[i] * (1 - b[i]) = 0
    for i in 0..bits {
        polynomials.push(format!("b[{}]^2 - b[{}]", i, i));
    }

    // 3. Boolean constraints for sum bits
    for i in 0..bits {
        polynomials.push(format!("sum[{}]^2 - sum[{}]", i, i));
    }

    // 4. Boolean constraints for carry bits
    for i in 0..bits {
        polynomials.push(format!("carry[{}]^2 - carry[{}]", i, i));
    }

    // 5. Carry[0] = 0 (no carry in)
    polynomials.push("carry[0]".to_string());

    // 6. XOR intermediate: xor_ab[i] = a[i] + b[i] - 2*a[i]*b[i]
    for i in 0..bits {
        polynomials.push(format!("xor_ab[{}] - a[{}] - b[{}] + 2*a[{}]*b[{}]", i, i, i, i, i));
    }

    // 7. Boolean constraints for xor_ab
    for i in 0..bits {
        polynomials.push(format!("xor_ab[{}]^2 - xor_ab[{}]", i, i));
    }

    // 8. Sum bit definition: sum[i] = xor_ab[i] XOR carry[i]
    //    = xor_ab[i] + carry[i] - 2*xor_ab[i]*carry[i]
    for i in 0..bits {
        polynomials.push(format!(
            "sum[{}] - xor_ab[{}] - carry[{}] + 2*xor_ab[{}]*carry[{}]",
            i, i, i, i, i
        ));
    }

    // 9. Next carry definition: carry[i+1] = MAJ(a[i], b[i], carry[i])
    //    = a[i]*b[i] + a[i]*carry[i] + b[i]*carry[i] - 2*a[i]*b[i]*carry[i]
    //    Simplified (for boolean): = a[i]*b[i] + (a[i] + b[i] - 2*a[i]*b[i])*carry[i]
    //                             = a[i]*b[i] + xor_ab[i]*carry[i]
    for i in 0..(bits-1) {
        polynomials.push(format!(
            "carry[{}] - a[{}]*b[{}] - xor_ab[{}]*carry[{}]",
            i+1, i, i, i, i
        ));
    }

    // 10. Carry out definition
    polynomials.push(format!(
        "carry_out - a[{}]*b[{}] - xor_ab[{}]*carry[{}]",
        bits-1, bits-1, bits-1, bits-1
    ));

    // Public names (the b[] outputs for verification)
    let public_names: Vec<String> = (0..bits)
        .map(|i| format!("b[{}]", i))
        .chain(std::iter::once("1".to_string()))
        .collect();

    let system = PolySystem {
        ring_config: RingConfig {
            n_vars,
            max_degree: 3, // We have degree 3 terms from a*b*carry
            mult_table: vec![1, 1],
            field: "BN254_Fr".to_string(),
        },
        variables: Variables {
            names: var_names,
            roles: var_roles,
            public_names,
        },
        system: System {
            n_constraints: polynomials.len(),
            polynomials,
        },
        gb_config: GBConfig {
            max_degree: 32,
            s_pair_budget: 5000000,
            monomial_order: "DegRevLex".to_string(),
            parallel: false,
        },
    };

    // Output JSON
    println!("{}", serde_json::to_string_pretty(&system).unwrap());
}
