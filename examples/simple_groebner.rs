//! Basic example of computing a Gröbner basis over a cryptographic field

use ark_feanor::*;
use feanor_math::algorithms::buchberger::*;
use feanor_math::rings::multivariate::*;

fn main() {
    println!("=== Gröbner Basis Example over BLS12-381 Fr ===\n");
    
    // Use the BLS12-381 scalar field
    let field = &*BLS12_381_FR;
    
    // Create a multivariate polynomial ring Fr[x, y]
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);
    
    println!("Working in polynomial ring: BLS12-381_Fr[x, y]");
    println!("Field characteristic has {} bits\n", {
        use prime_field::FieldProperties;
        let wrapper = ArkFieldWrapper::<ark_bls12_381::Fr>::new();
        wrapper.characteristic().bits()
    });
    
    // Define a polynomial system:
    // - x^2 + y^2 - 1 = 0  (circle equation)
    // - x*y - 2 = 0         (hyperbola equation)
    let system = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        vec![
            x.clone().pow(2) + y.clone().pow(2) - 1,
            x.clone() * y.clone() - 2,
        ]
    });
    
    println!("Input polynomial system:");
    for (i, poly) in system.iter().enumerate() {
        print!("  f_{} = ", i + 1);
        print_polynomial(&poly_ring, poly);
    }
    println!();
    
    // Compute Gröbner basis using degree reverse lexicographic ordering
    println!("Computing Gröbner basis with DegRevLex ordering...");
    
    match buchberger(&poly_ring, system.clone(), DegRevLex, |_| {}) {
        Ok(gb) => {
            println!("Gröbner basis computed successfully!");
            println!("Basis has {} polynomials:\n", gb.len());
            
            for (i, poly) in gb.iter().enumerate() {
                print!("  g_{} = ", i + 1);
                print_polynomial(&poly_ring, poly);
            }
            
            // Check if the system has solutions
            if gb.len() == 1 && poly_ring.is_one(&gb[0]) {
                println!("\nThe system has NO solutions (basis is {1})");
            } else {
                println!("\nThe system may have solutions");
                println!("Further analysis would be needed to find them");
            }
        }
        Err(e) => {
            println!("Failed to compute Gröbner basis: {:?}", e);
        }
    }
    
    println!("\n=== Example with a simpler system ===\n");
    
    // Another example: linear system
    // x + y - 3 = 0
    // x - y - 1 = 0
    // Solution: x = 2, y = 1
    let linear_system = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        vec![
            x.clone() + y.clone() - 3,
            x.clone() - y.clone() - 1,
        ]
    });
    
    println!("Linear system:");
    for (i, poly) in linear_system.iter().enumerate() {
        print!("  f_{} = ", i + 1);
        print_polynomial(&poly_ring, poly);
    }
    println!();
    
    match buchberger(&poly_ring, linear_system, Lex, |_| {}) {
        Ok(gb) => {
            println!("Gröbner basis (Lex ordering):");
            for (i, poly) in gb.iter().enumerate() {
                print!("  g_{} = ", i + 1);
                print_polynomial(&poly_ring, poly);
            }
            
            // Try to extract solution if possible
            println!("\nAnalyzing the basis for solutions...");
            analyze_linear_solution(&poly_ring, &gb);
        }
        Err(e) => {
            println!("Failed to compute Gröbner basis: {:?}", e);
        }
    }
}

fn print_polynomial<R: MultivariatePolyRingStore>(
    ring: &R,
    poly: &R::Element
) where
    R::Type: MultivariatePolyRing,
{
    let terms: Vec<_> = ring.terms(poly).collect();
    
    if terms.is_empty() {
        println!("0");
        return;
    }
    
    for (i, (coeff, monomial)) in terms.iter().enumerate() {
        // Print coefficient
        let coeff_str = format!("{:?}", coeff);
        
        // Print monomial
        let vars = ["x", "y", "z", "w"];  // Variable names
        let mut monomial_parts = Vec::new();
        
        for (var_idx, &power) in monomial.exponents().iter().enumerate() {
            if power > 0 {
                if power == 1 {
                    monomial_parts.push(vars[var_idx].to_string());
                } else {
                    monomial_parts.push(format!("{}^{}", vars[var_idx], power));
                }
            }
        }
        
        // Format the term
        if monomial_parts.is_empty() {
            // Constant term
            print!("{}", coeff_str);
        } else if coeff_str == "1" {
            print!("{}", monomial_parts.join("*"));
        } else if coeff_str == "-1" {
            print!("-{}", monomial_parts.join("*"));
        } else {
            print!("{}*{}", coeff_str, monomial_parts.join("*"));
        }
        
        if i < terms.len() - 1 {
            print!(" + ");
        }
    }
    println!();
}

fn analyze_linear_solution<R: MultivariatePolyRingStore>(
    ring: &R,
    gb: &[R::Element]
) where
    R::Type: MultivariatePolyRing,
{
    // For a linear system, the Gröbner basis might directly give us the solution
    // Look for polynomials of the form "x - c" or "y - c"
    
    for poly in gb {
        let terms: Vec<_> = ring.terms(poly).collect();
        
        // Check if it's a linear polynomial in one variable
        if terms.len() == 2 {
            // Might be of the form "var - constant"
            // This is a simple heuristic; real implementation would be more sophisticated
            println!("  Found potential solution polynomial with {} terms", terms.len());
        }
    }
    
    println!("  (Full solution extraction would require more sophisticated analysis)");
}
