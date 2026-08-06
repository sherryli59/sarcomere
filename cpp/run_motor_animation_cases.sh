#!/usr/bin/env bash
set -u

binary=${1:-cpp/build_titin_gate_20260802/motor_validation}
output_dir=${2:-cpp/diagnostics/self_loaded_balance_20260804/animations}
seed=${3:-2080206}

mkdir -p "$output_dir"

run_case() {
    local label=$1
    local mode=$2
    local duration_ms=$3
    local heads=$4
    local preeq_ms=$5
    local save_every=$6
    local prefix="$output_dir/$label"
    OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE "$binary" \
        "$mode" "$duration_ms" "$seed" "$heads" 0 0.5 -1 "$preeq_ms" \
        0.16 0.01 100 3 1 "$prefix.csv" "$prefix.h5" "$save_every" \
        >"$prefix.log" 2>&1
    local result=$?
    printf '%s exit=%d\n' "$label" "$result"
    tail -n 3 "$prefix.log"
}

run_case free_100ms free 100 35 0 1000
run_case one_sided_300ms linked 300 35 0 3000
run_case balanced_stall_100ms balanced 100 70 30 1000
# This trajectory is intentionally expected to fail the all-actin A-M gate;
# it exposes the ACTN2-invisible pair-centroid mode.
run_case balanced_centroid_failure_300ms balanced 300 70 30 3000
