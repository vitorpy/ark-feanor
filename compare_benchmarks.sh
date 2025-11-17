#!/bin/bash
# Quick benchmark comparison script
# Compares F4 vs Buchberger and AVX512 vs baseline

set -e

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   Benchmark Comparison Runner${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Function to extract mean time in milliseconds
get_time_ms() {
    local bench_group=$1
    local bench_name=$2
    local json_file="target/criterion/${bench_group}/${bench_name}/base/estimates.json"

    if [ -f "$json_file" ]; then
        # Extract mean point estimate in nanoseconds, convert to milliseconds
        jq -r '.mean.point_estimate / 1000000' "$json_file" 2>/dev/null || echo "N/A"
    else
        echo "N/A"
    fi
}

# Function to calculate speedup
calc_speedup() {
    local baseline=$1
    local optimized=$2

    if [ "$baseline" != "N/A" ] && [ "$optimized" != "N/A" ]; then
        echo "scale=2; $baseline / $optimized" | bc
    else
        echo "N/A"
    fi
}

# Function to print time value (handles N/A)
print_time() {
    local label=$1
    local value=$2

    if [ "$value" = "N/A" ]; then
        printf "%-30s %10s\n" "$label" "N/A"
    else
        printf "%-30s %10.2f ms\n" "$label" "$value"
    fi
}

echo -e "${YELLOW}Reading existing benchmark results...${NC}"
echo ""

# Uncomment to re-run benchmarks:
# echo -e "${GREEN}1/6${NC} Running field_ops (baseline)..."
# cargo bench --bench field_ops --quiet > /dev/null 2>&1 || true
#
# echo -e "${GREEN}2/6${NC} Running field_ops_avx512..."
# cargo bench --bench field_ops_avx512 --features avx512-ifma --quiet > /dev/null 2>&1 || true
#
# echo -e "${GREEN}3/6${NC} Running buchberger_tiny..."
# cargo bench --bench buchberger_tiny --quiet > /dev/null 2>&1 || true
#
# echo -e "${GREEN}4/6${NC} Running buchberger_tiny_avx512..."
# cargo bench --bench buchberger_tiny_avx512 --features avx512-ifma --quiet > /dev/null 2>&1 || true
#
# echo -e "${GREEN}5/6${NC} Running f4_tiny..."
# cargo bench --bench f4_tiny --quiet > /dev/null 2>&1 || true
#
# echo -e "${GREEN}6/6${NC} Running f4_tiny_avx512..."
# cargo bench --bench f4_tiny_avx512 --features avx512-ifma --quiet > /dev/null 2>&1 || true


echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   RESULTS${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Field operations comparison
echo -e "${YELLOW}Field Operations (8-element batch, 100K iterations):${NC}"
echo "----------------------------------------"

seq_mul=$(get_time_ms "field_mul_sequential" "bn254_mul_8elem_100K_sequential")
batch_mul=$(get_time_ms "field_mul_avx512" "bn254_mul_8elem_100K_avx512_batch")
speedup=$(calc_speedup "$seq_mul" "$batch_mul")

print_time "Sequential (8 muls):" "$seq_mul"
print_time "AVX512 batch (8 muls):" "$batch_mul"
if [ "$speedup" != "N/A" ]; then
    printf "%-30s ${GREEN}%10.2fx${NC}\n" "Speedup:" "$speedup"
fi
echo ""

# Cyclic-3 comparison
echo -e "${YELLOW}Cyclic-3 (3 variables):${NC}"
echo "----------------------------------------"

buch3=$(get_time_ms "buchberger_cyclic3" "cyclic3_degrevlex")
f4_3=$(get_time_ms "f4_cyclic3" "cyclic3_degrevlex")
speedup=$(calc_speedup "$buch3" "$f4_3")

print_time "Buchberger:" "$buch3"
print_time "F4:" "$f4_3"
if [ "$speedup" != "N/A" ]; then
    if (( $(echo "$speedup > 1.0" | bc -l) )); then
        printf "%-30s ${GREEN}%10.2fx${NC} (F4 faster)\n" "Speedup:" "$speedup"
    else
        printf "%-30s ${RED}%10.2fx${NC} (Buchberger faster)\n" "Speedup:" "$speedup"
    fi
fi
echo ""

# Cyclic-4 comparison
echo -e "${YELLOW}Cyclic-4 (4 variables):${NC}"
echo "----------------------------------------"

buch4=$(get_time_ms "buchberger_cyclic4" "cyclic4_degrevlex")
f4_4=$(get_time_ms "f4_cyclic4" "cyclic4_degrevlex")
speedup=$(calc_speedup "$buch4" "$f4_4")

print_time "Buchberger:" "$buch4"
print_time "F4:" "$f4_4"
if [ "$speedup" != "N/A" ]; then
    if (( $(echo "$speedup > 1.0" | bc -l) )); then
        printf "%-30s ${GREEN}%10.2fx${NC} (F4 faster)\n" "Speedup:" "$speedup"
    else
        printf "%-30s ${RED}%10.2fx${NC} (Buchberger faster)\n" "Speedup:" "$speedup"
    fi
fi
echo ""

# Katsura-3 comparison
echo -e "${YELLOW}Katsura-3 (3 variables):${NC}"
echo "----------------------------------------"

buch_k3=$(get_time_ms "buchberger_katsura3" "katsura3_degrevlex")
f4_k3=$(get_time_ms "f4_katsura3" "katsura3_degrevlex")
speedup=$(calc_speedup "$buch_k3" "$f4_k3")

print_time "Buchberger:" "$buch_k3"
print_time "F4:" "$f4_k3"
if [ "$speedup" != "N/A" ]; then
    if (( $(echo "$speedup > 1.0" | bc -l) )); then
        printf "%-30s ${GREEN}%10.2fx${NC} (F4 faster)\n" "Speedup:" "$speedup"
    else
        printf "%-30s ${RED}%10.2fx${NC} (Buchberger faster)\n" "Speedup:" "$speedup"
    fi
fi
echo ""

# Katsura-4 comparison
echo -e "${YELLOW}Katsura-4 (4 variables):${NC}"
echo "----------------------------------------"

buch_k4=$(get_time_ms "buchberger_katsura4" "katsura4_degrevlex")
f4_k4=$(get_time_ms "f4_katsura4" "katsura4_degrevlex")
speedup=$(calc_speedup "$buch_k4" "$f4_k4")

print_time "Buchberger:" "$buch_k4"
print_time "F4:" "$f4_k4"
if [ "$speedup" != "N/A" ]; then
    if (( $(echo "$speedup > 1.0" | bc -l) )); then
        printf "%-30s ${GREEN}%10.2fx${NC} (F4 faster)\n" "Speedup:" "$speedup"
    else
        printf "%-30s ${RED}%10.2fx${NC} (Buchberger faster)\n" "Speedup:" "$speedup"
    fi
fi
echo ""

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   SUMMARY${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Determine overall winner
f4_wins=0
buch_wins=0

for system in "cyclic3" "cyclic4" "katsura3" "katsura4"; do
    if [ "$system" = "cyclic3" ]; then
        b="$buch3"
        f="$f4_3"
    elif [ "$system" = "cyclic4" ]; then
        b="$buch4"
        f="$f4_4"
    elif [ "$system" = "katsura3" ]; then
        b="$buch_k3"
        f="$f4_k3"
    else
        b="$buch_k4"
        f="$f4_k4"
    fi

    if [ "$b" != "N/A" ] && [ "$f" != "N/A" ]; then
        if (( $(echo "$f < $b" | bc -l) )); then
            ((f4_wins++))
        else
            ((buch_wins++))
        fi
    fi
done

echo "Algorithm Performance:"
printf "  F4 wins:        %d/4 tests\n" "$f4_wins"
printf "  Buchberger wins: %d/4 tests\n" "$buch_wins"
echo ""

if [ "$speedup" != "N/A" ] && (( $(echo "$speedup > 1.5" | bc -l) )); then
    echo -e "${GREEN}✓ AVX512 batch operations show good speedup (${speedup}x)${NC}"
    echo "  Consider integrating into F4 matrix reduction"
else
    echo -e "${YELLOW}⚠ AVX512 benefit unclear from current benchmarks${NC}"
fi
echo ""

echo "Full results: target/criterion/"
echo ""
