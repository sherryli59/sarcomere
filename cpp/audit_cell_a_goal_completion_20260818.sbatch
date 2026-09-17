#!/bin/bash
# Final requirement-by-requirement Cell A goal-completion audit.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=audit_cellA_goal
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

base=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818
snapshot=${base}/source_snapshot_nobarrier
paired=${base}/full5s
corrected=${base}/full5s_nobarrier
compat=/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818/myomesin_retained_nobarrier
binary=/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818/bin/sarcomere_kmm0_nobarrier_20260818
sha=5ce40655c5349dc5de079713ea29b32984ef5b4a803264d9e77201a63427f51e
tag=A800_M320_cellA_control_seed20260805_5000ms

python "${snapshot}/audit_cell_a_goal_completion.py" \
    --replication "${corrected}/analysis/cell_a_replication.json" \
    --paired-audit "${paired}/analysis/cellA_sparse_myomesin_audit_v2.json" \
    --corrected-health "${corrected}/analysis/${tag}.live_final.json" \
    --corrected-params "${corrected}/params/${tag}.params" \
    --expected-binary "${binary}" \
    --expected-binary-sha256 "${sha}" \
    --treatment-reference "${compat}/reference_trajectory.h5" \
    --treatment-candidate "${compat}/candidate_trajectory.h5" \
    --output "${base}/cell_a_goal_completion.json"
