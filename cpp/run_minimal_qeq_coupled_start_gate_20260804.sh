#!/bin/bash
set -euo pipefail

cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static

exec cpp/build_titin_gate_20260802/motor_validation \
    balanced 1000 20260804 35 0 0.5 -1 0 0.16 0.010 100 3 1 \
    cpp/diagnostics/self_loaded_balance_20260804/qeq_coupled_start/minimal_qeq.csv \
    cpp/diagnostics/self_loaded_balance_20260804/qeq_coupled_start/minimal_qeq.h5 \
    2000 100 0.64 0.80 5000 0.05 3 0.95 1
