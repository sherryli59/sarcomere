#!/bin/bash
set -euo pipefail

cd /home/groups/rotskoff/sherry/sarcomere/sarcomere

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static

output_dir=cpp/diagnostics/productive_register_qeq_1s_20260805
mkdir -p "${output_dir}"

# One-second extension of the q-equilibrium minimal balanced gate. All model
# parameters and the seed match the 100 ms run; only duration changes.
exec cpp/build_split_detachment_20260805/motor_validation \
    balanced 1000 20260805 35 0 0.5 -1 0 0.16 0.010 100 3 1 \
    "${output_dir}/minimal_productive_register_qeq_1s.csv" \
    "${output_dir}/minimal_productive_register_qeq_1s.h5" \
    2000 100 0.64 0.80 5000 0.05 3 0.95 1 100 120 0
