#!/usr/bin/env bash
set -u

binary=${1:-cpp/build_titin_gate_20260802/motor_validation}
output_dir=${2:-cpp/diagnostics/self_loaded_balance_20260804/suite}
seed=${3:-2080206}
thermal_scale=${4:-1}
heads_per_half=${5:-35}
balanced_duration_ms=${6:-100}

mkdir -p "$output_dir"

failures=0
run_case() {
    local mode=$1
    local duration_ms=$2
    local preeq_ms=$3
    local log="$output_dir/${mode}_seed${seed}_${duration_ms}ms.log"
    local csv="$output_dir/${mode}_seed${seed}_${duration_ms}ms.csv"
    OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE "$binary" \
        "$mode" "$duration_ms" "$seed" "$heads_per_half" 0 0.5 -1 "$preeq_ms" \
        0.16 0.01 100 3 "$thermal_scale" "$csv" >"$log" 2>&1
    local result=$?
    tail -n 9 "$log"
    if (( result != 0 )); then
        failures=$((failures + 1))
    fi
}

run_case free 100 0
run_case linked 300 0
run_case balanced "$balanced_duration_ms" 30

if (( failures != 0 )); then
    printf 'SUITE FAIL: %d/3 scenarios failed\n' "$failures"
    exit 1
fi
printf 'SUITE PASS: 3/3 scenarios passed\n'
