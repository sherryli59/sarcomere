#!/bin/bash
set -euo pipefail

cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static

# Preferred balanced start: A-A overlap is already 0.10 um.  The intended
# titin endpoint separation is about 0.631 um, inside the 0.64 um slack length,
# so newly captured titin molecules initially exert no force.
exec cpp/build_titin_gate_20260802/motor_validation \
    balanced 1000 20260804 35 0 0.5 -1 0 0.10 0.010 100 3 1 \
    cpp/diagnostics/self_loaded_balance_20260804/distance_capture_perfect_start/perfect_1000ms.csv \
    cpp/diagnostics/self_loaded_balance_20260804/distance_capture_perfect_start/perfect_1000ms.h5 \
    2000 100 0.64 0.80 5000 0.05 3 0.95 1
