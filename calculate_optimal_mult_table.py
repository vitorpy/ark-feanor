#!/usr/bin/env python3
"""
Calculate optimal multiplication table configurations for specific memory budgets.
For n=1142 variables with 50GB and 100GB memory constraints.
"""

from math import comb

def num_monomials(n, d):
    """Calculate number of monomials for n variables at degree d."""
    if d == 0:
        return 1
    return comb(n + d - 1, d)

def mult_table_memory(n, d1, d2):
    """
    Calculate total memory for multiplication table with max degrees (d1, d2).
    Table structure: [lhs_deg][rhs_deg][lhs_index][rhs_index] -> u64
    """
    total_bytes = 0
    for lhs_deg in range(d1 + 1):
        for rhs_deg in range(d2 + 1):
            count_lhs = num_monomials(n, lhs_deg)
            count_rhs = num_monomials(n, rhs_deg)
            total_bytes += count_lhs * count_rhs * 8  # u64 = 8 bytes
    return total_bytes

def format_bytes(bytes_val):
    """Format bytes in human-readable form."""
    if bytes_val < 1024:
        return f"{bytes_val} B"
    elif bytes_val < 1024**2:
        return f"{bytes_val / 1024:.2f} KB"
    elif bytes_val < 1024**3:
        return f"{bytes_val / (1024**2):.2f} MB"
    elif bytes_val < 1024**4:
        return f"{bytes_val / (1024**3):.2f} GB"
    else:
        return f"{bytes_val / (1024**4):.2f} TB"

def find_optimal_config(n, max_bytes, max_degree=7):
    """
    Find the optimal (d1, d2) configuration that fits within max_bytes.
    Returns the configuration that maximizes d1 + d2 while staying under budget.
    """
    best_config = (0, 0)
    best_memory = 0
    best_sum = 0

    print(f"\nSearching for optimal config with n={n}, budget={format_bytes(max_bytes)}, max_degree={max_degree}")
    print("=" * 80)

    # Try all combinations up to max_degree
    candidates = []
    for d1 in range(max_degree + 1):
        for d2 in range(d1, max_degree + 1):  # d2 >= d1 is typical
            memory = mult_table_memory(n, d1, d2)
            if memory <= max_bytes:
                candidates.append((d1, d2, memory, d1 + d2))
                if d1 + d2 > best_sum or (d1 + d2 == best_sum and memory > best_memory):
                    best_config = (d1, d2)
                    best_memory = memory
                    best_sum = d1 + d2

    # Sort candidates by sum of degrees (descending), then by memory (descending)
    candidates.sort(key=lambda x: (x[3], x[2]), reverse=True)

    # Show top 10 candidates
    print(f"\nTop configurations that fit within budget:")
    print(f"{'Config':<12} {'Memory':<15} {'Sum':<8} {'Utilization':<12}")
    print("-" * 80)
    for i, (d1, d2, mem, sum_deg) in enumerate(candidates[:10]):
        utilization = (mem / max_bytes) * 100
        print(f"({d1:2}, {d2:2})      {format_bytes(mem):<15} {sum_deg:<8} {utilization:6.2f}%")

    return best_config, best_memory

def main():
    n = 1142  # number of variables
    max_degree = 7  # absolute limit due to u64 overflow at degree 8

    print("=" * 80)
    print(f"MULTIPLICATION TABLE MEMORY CALCULATOR")
    print(f"Variables: {n}")
    print(f"Max degree: {max_degree} (u64 limit)")
    print("=" * 80)

    # Calculate some reference configurations
    print("\nReference configurations:")
    print(f"{'Config':<12} {'Memory':<15}")
    print("-" * 40)
    for d1, d2 in [(0, 0), (1, 1), (2, 2), (2, 3), (2, 4), (2, 5)]:
        mem = mult_table_memory(n, d1, d2)
        print(f"({d1}, {d2})        {format_bytes(mem):<15}")

    # 50GB budget
    budget_50gb = 50 * (1024**3)
    config_50, mem_50 = find_optimal_config(n, budget_50gb, max_degree)

    print("\n" + "=" * 80)
    print(f"OPTIMAL CONFIG FOR 50GB BUDGET:")
    print(f"  Config: ({config_50[0]}, {config_50[1]})")
    print(f"  Memory: {format_bytes(mem_50)}")
    print(f"  Utilization: {(mem_50 / budget_50gb) * 100:.2f}%")
    print("=" * 80)

    # 100GB budget
    budget_100gb = 100 * (1024**3)
    config_100, mem_100 = find_optimal_config(n, budget_100gb, max_degree)

    print("\n" + "=" * 80)
    print(f"OPTIMAL CONFIG FOR 100GB BUDGET:")
    print(f"  Config: ({config_100[0]}, {config_100[1]})")
    print(f"  Memory: {format_bytes(mem_100)}")
    print(f"  Utilization: {(mem_100 / budget_100gb) * 100:.2f}%")
    print("=" * 80)

    # Show monomial counts for context
    print(f"\nMonomial counts for n={n}:")
    print(f"{'Degree':<8} {'Monomials':<20}")
    print("-" * 40)
    for d in range(8):
        count = num_monomials(n, d)
        print(f"{d:<8} {count:>20,}")

    print("\n" + "=" * 80)
    print("RECOMMENDATIONS:")
    print("=" * 80)
    print(f"For 50GB:  Use config ({config_50[0]}, {config_50[1]}) - {format_bytes(mem_50)}")
    print(f"For 100GB: Use config ({config_100[0]}, {config_100[1]}) - {format_bytes(mem_100)}")
    print("\nCode example:")
    print(f"""
let poly_ring = MultivariatePolyRingImpl::new_with_mult_table(
    field,
    {n},          // variables
    {max_degree},            // max degree
    ({config_50[0]}, {config_50[1]}),      // mult table config for 50GB
    Global
);
""")

if __name__ == "__main__":
    main()
