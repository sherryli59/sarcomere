#!/bin/bash
# End-to-end historical Cell A compatibility probe. Both binaries run the
# exact historical command at one OpenMP thread; only output cadence/duration
# are shortened. The comparison requires exact equality of molecular state.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --job-name=cellA_compat
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818/slurm/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818/slurm/%x_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

jemalloc=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${jemalloc}" || { echo "missing jemalloc" >&2; exit 1; }
export LD_PRELOAD="${jemalloc}"

repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818
probe=${CELL_A_COMPAT_PROBE:-initial}
probe_steps=${CELL_A_COMPAT_STEPS:-2000}
probe_save_every=${CELL_A_COMPAT_SAVE_EVERY:-500}
outdir=${root}/${probe}
old_exe=${CELL_A_COMPAT_OLD_EXE:-/scratch/users/xsli/sarcomere/prod_brake_A800_M320_5s_20260807/bin/sarcomere_dt1e5_20260807}
new_exe=${CELL_A_COMPAT_NEW_EXE:-/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/bin/sarcomere_sparse_myomesin_module_20260818}
reference_log=${CELL_A_COMPAT_REFERENCE_LOG:-/scratch/users/xsli/sarcomere/prod_brake_A800_M320_5s_20260807/logs/sarc_brake_A800_M320_5s_38011216.out}
old_out=${outdir}/reference_trajectory.h5
new_out=${outdir}/candidate_trajectory.h5

mkdir -p "${root}/slurm" "${outdir}"
for path in "${old_exe}" "${new_exe}" "${reference_log}"; do
    test -r "${path}" || { echo "missing input: ${path}" >&2; exit 1; }
done
for path in "${old_out}" "${new_out}"; do
    test ! -e "${path}" || { echo "refusing to overwrite ${path}" >&2; exit 1; }
done

# Recover the authoritative historical argument vector directly from its log.
# Output path, duration, and save cadence are the only substitutions.
args=()
while IFS= read -r token; do
    args+=("${token}")
done < <(python - "${reference_log}" <<'PY'
import shlex
import sys

line = next(line for line in open(sys.argv[1]) if line.startswith("command="))
tokens = shlex.split(line[len("command="):])
start = next(i for i, token in enumerate(tokens) if token.startswith("--n_actins="))
args = []
for token in tokens[start:]:
    if token.startswith(("--filename=", "--nsteps=", "--save_every=")):
        continue
    args.append(token)
for value in args:
    print(value)
PY
)

echo "old_sha256=$(sha256sum "${old_exe}" | awk '{print $1}')"
echo "new_sha256=$(sha256sum "${new_exe}" | awk '{print $1}')"
echo "probe=${probe} args=${#args[@]} OMP_NUM_THREADS=${OMP_NUM_THREADS} nsteps=${probe_steps} save_every=${probe_save_every}"

srun --kill-on-bad-exit=1 "${old_exe}" "${args[@]}" \
    --nsteps="${probe_steps}" --save_every="${probe_save_every}" --filename="${old_out}"
srun --kill-on-bad-exit=1 "${new_exe}" "${args[@]}" \
    --nsteps="${probe_steps}" --save_every="${probe_save_every}" --filename="${new_out}"

python "${repo}/analysis/compare_cell_a_backward_compat.py" \
    "${old_out}" "${new_out}" | tee "${outdir}/comparison.txt"
