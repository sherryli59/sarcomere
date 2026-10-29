#!/bin/bash
#SBATCH --job-name=actn2-standalone
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

# Submit from the repository root, specifying a new results directory and log:
# sbatch --output="$SCRATCH/actn2-standalone-%j.log" \
#   analysis/run_standalone_actn2_contact.sbatch "$SCRATCH/actn2-new-results"
# Aggregates the full experiment and the multi-seed refinement checks.
set -euo pipefail
result_dir=${1:?Specify a new results directory, preferably under SCRATCH}
test ! -e "$result_dir" || { echo "Refusing to overwrite $result_dir" >&2; exit 1; }
source_dir="${SLURM_SUBMIT_DIR:?}/analysis"
python_bin="${GROUP_HOME:?}/sherry/mamba_root/envs/sarcomere/bin/python"
test -f "$source_dir/standalone_actn2_contact.py"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR="${SCRATCH:?}/actn2-standalone-mpl-${SLURM_JOB_ID:?}"
mkdir -p "$result_dir"
"$python_bin" -u "$source_dir/standalone_actn2_contact.py" \
  --output "$result_dir/primary" --replicates 64 --seed 20261028
"$python_bin" -u "$source_dir/validate_standalone_actn2_contact.py" \
  --output "$result_dir/refinement" --replicates 256
