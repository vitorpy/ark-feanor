//! Example of solving polynomial systems over cryptographic fields

use ark_feanor::*;
use feanor_math::algorithms::buchberger::*;
use feanor_math::rings::multivariate::*;
use std::time::Instant;

fn main() {
    println!("=== Polynomial System Solver ===\n");
    
    example_linear_system();
    println!("\n" + "=".repeat(50) + "\n");
    example_circle_parabola();
    println!("\n" + "=".repeat(50) + "\n");
    example_three_variables();
}

/// Solve a simple linear system
fn example_linear_system() {
    println!("Example 1: Linear System over BN254_Fr");
    println!("----------------------------------------");
    
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);
    
    // System:
    // 2x + 3y = 7
    // x - y = 1
    // Solution: x = 2, y = 1
    let system = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        vec![
            2 * x.clone() + 3 * y.clone() - 7,
            x.clone() - y.clone() - 1,
        ]
    });
    
    println!("System of equations:");
    println!("  2x + 3y - 7 = 0");
    println!("  x - y - 1 = 0");
    println!();
    
    let start = Instant::now();
    match buchberger(&poly_ring, system, Lex, |_| {}) {
        Ok(gb) => {
            let elapsed = start.elapsed();
            println!("Gröbner basis computed in {:?}", elapsed);
            println!("Number of polynomials in basis: {}", gb.len());
            
            // In a linear system with Lex ordering, we often get
            // polynomials like "y - 1" and "x - 2"
            for poly in &gb {
                let terms: Vec<_> = poly_ring.terms(poly).collect();
                if terms.len() <= 2 {
                    println!("  Potential solution polynomial found");
                }
            }
        }
        Err(e) => println!("Error: {:?}", e),
    }
}

/// Solve the intersection of a circle and parabola
fn example_circle_parabola() {
    println!("Example 2: Circle-Parabola Intersection over BLS12-381_Fr");
    println!("----------------------------------------------------------");
    
    let field = &*BLS12_381_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 2);
    
    // System:
    // x² + y² = 4  (circle with radius 2)
    // y = x²       (parabola)
    let system = poly_ring.with_wrapped_indeterminates(|[x, y]| {
        vec![
            x.clone().pow(2) + y.clone().pow(2) - 4,
            y.clone() - x.clone().pow(2),
        ]
    });
    
    println!("System of equations:");
    println!("  x² + y² - 4 = 0  (circle)");
    println!("  y - x² = 0        (parabola)");
    println!();
    
    // Try with different monomial orderings
    for (name, ordering) in &[("Lex", Lex), ("DegRevLex", DegRevLex)] {
        println!("Using {} ordering:", name);
        let start = Instant::now();
        
        match buchberger(&poly_ring, system.clone(), *ordering, |_| {}) {
            Ok(gb) => {
                let elapsed = start.elapsed();
                println!("  Computed in {:?}", elapsed);
                println!("  Basis size: {}", gb.len());
                
                // Check if system is consistent
                if gb.len() == 1 && poly_ring.is_one(&gb[0]) {
                    println!("  System is INCONSISTENT (no solutions)");
                } else {
                    println!("  System may have solutions");
                    
                    // Look for univariate polynomials in the basis
                    for poly in &gb {
                        if is_univariate(&poly_ring, poly) {
                            println!("    Found univariate polynomial (can solve directly)");
                        }
                    }
                }
            }
            Err(e) => println!("  Error: {:?}", e),
        }
        println!();
    }
}

/// Solve a system with three variables
fn example_three_variables() {
    println!("Example 3: Three-Variable System over BN254_Fr");
    println!("-----------------------------------------------");
    
    let field = &*BN254_FR;
    let poly_ring = MultivariatePolyRingImpl::new(field, 3);
    
    // System of symmetric polynomials
    // x + y + z = 6
    // xy + yz + xz = 11
    // xyz = 6
    // (Solution: permutations of {1, 2, 3})
    let system = poly_ring.with_wrapped_indeterminates(|[x, y, z]| {
        vec![
            x.clone() + y.clone() + z.clone() - 6,
            x.clone() * y.clone() + y.clone() * z.clone() + x.clone() * z.clone() - 11,
            x.clone() * y.clone() * z.clone() - 6,
        ]
    });
    
    println!("System of equations (symmetric polynomials):");
    println!("  x + y + z = 6");
    println!("  xy + yz + xz = 11");
    println!("  xyz = 6");
    println!();
    println!("Expected solutions: permutations of {{1, 2, 3}}");
    println!();
    
    let start = Instant::now();
    match buchberger(&poly_ring, system, DegRevLex, |step| {
        if step % 10 == 0 {
            print!(".");
            use std::io::{self, Write};
            io::stdout().flush().unwrap();
        }
    }) {
        Ok(gb) => {
            println!();
            let elapsed = start.elapsed();
            println!("Gröbner basis computed in {:?}", elapsed);
            println!("Basis contains {} polynomials", gb.len());
            
            // Analyze the structure of the basis
            let mut min_terms = usize::MAX;
            let mut max_degree = 0;
            
            for poly in &gb {
                let terms: Vec<_> = poly_ring.terms(poly).collect();
                min_terms = min_terms.min(terms.len());
                
                for (_, monomial) in &terms {
                    max_degree = max_degree.max(monomial.total_degree());
                }
            }
            
            println!("Basis analysis:");
            println!("  Minimum terms in a polynomial: {}", min_terms);
            println!("  Maximum total degree: {}", max_degree);
            
            if min_terms == 1 {
                println!("  Note: Found constant polynomial(s)");
            }
        }
        Err(e) => println!("Error computing basis: {:?}", e),
    }
}

/// Check if a polynomial is univariate (involves only one variable)
fn is_univariate<R: MultivariatePolyRingStore>(
    ring: &R,
    poly: &R::Element
) -> bool
where
    R::Type: MultivariatePolyRing,
{
    let mut active_vars = vec![false; ring.variable_count()];
    
    for (_, monomial) in ring.terms(poly) {
        for (i, &exp) in monomial.exponents().iter().enumerate() {
            if exp > 0 {
                active_vars[i] = true;
            }
        }
    }
    
    active_vars.iter().filter(|&&v| v).count() == 1
}
