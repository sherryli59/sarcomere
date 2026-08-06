#!/usr/bin/env bash
set -u

trajectory_dir=${1:-cpp/diagnostics/self_loaded_balance_20260804/animations}
renderer=analysis/animate_minimal_mechanics.py
export MPLCONFIGDIR=/tmp/sarcomere_matplotlib
export HDF5_USE_FILE_LOCKING=FALSE

render() {
    local stem=$1
    local mode=$2
    local save_every=$3
    python "$renderer" "$trajectory_dir/$stem.h5" \
        "$trajectory_dir/$stem.mp4" --mode "$mode" \
        --save-every "$save_every" --stride 4 --fps 12 \
        >"$trajectory_dir/${stem}_render.log" 2>&1
}

render free_100ms single 1000 &
pid_free=$!
render one_sided_300ms catch 3000 &
pid_linked=$!
render balanced_stall_100ms catch 1000 &
pid_stall=$!
render balanced_centroid_failure_300ms catch 3000 &
pid_failure=$!

failures=0
for pid in "$pid_free" "$pid_linked" "$pid_stall" "$pid_failure"; do
    if ! wait "$pid"; then
        failures=$((failures + 1))
    fi
done

if (( failures != 0 )); then
    printf 'render failures: %d\n' "$failures"
    exit 1
fi
ls -lh "$trajectory_dir"/*.mp4
