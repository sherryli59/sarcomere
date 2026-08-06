#!/bin/bash
#SBATCH -p owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --job-name=prodnew_safe
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prodnew/logs/%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prodnew/logs/%j.err

set -euo pipefail

# Submit either trajectory with:
#   sbatch --export=ALL,PROD_VARIANT=new cpp/prodnew_checkpointed.sbatch
#   sbatch --export=ALL,PROD_VARIANT=reg cpp/prodnew_checkpointed.sbatch
#
# A Slurm requeue reruns this script. Existing output is accepted only when it
# has a committed atomic resume snapshot; otherwise the job exits without
# touching the file.

VARIANT="${PROD_VARIANT:-new}"
ROOT=/scratch/users/xsli/sarcomere/prodnew
EXE=/scratch/users/xsli/sarcomere/prod_L7p72_20260723/bin/sarcomere_thr2
BASE_FILE=/scratch/users/xsli/prod_verbatim.txt
TARGET_STEPS=400000

case "${VARIANT}" in
    new)
        OUTPUT="${ROOT}/data/prodnew.h5"
        REGISTER_RATIO=0.0
        ;;
    reg)
        OUTPUT="${ROOT}/data/prodnew_reg.h5"
        REGISTER_RATIO=0.3
        ;;
    *)
        echo "Unknown PROD_VARIANT=${VARIANT}; expected new or reg." >&2
        exit 2
        ;;
esac

test -x "${EXE}" || { echo "Missing executable: ${EXE}" >&2; exit 1; }
test -r "${BASE_FILE}" || { echo "Missing parameter file: ${BASE_FILE}" >&2; exit 1; }
mkdir -p "${ROOT}/data" "${ROOT}/logs"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export LD_PRELOAD=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "Another writer owns ${OUTPUT}; exiting." >&2; exit 0; }

resume_flag=()
if [ -e "${OUTPUT}" ]; then
    checkpoint_step=$(
        /home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin/python \
            - "${OUTPUT}" <<'PY'
import sys

import h5py

with h5py.File(sys.argv[1], "r", locking=False) as handle:
    steps = []
    for path in ("/resume/current_step", "/resume_prev/current_step"):
        if path in handle:
            values = handle[path][...].reshape(-1)
            if values.size:
                steps.append(int(values[-1]))
    if not steps or "/actin/center" not in handle:
        raise SystemExit(3)
    print(max(steps))
PY
    ) || {
        echo "Existing output has no valid atomic resume snapshot; refusing to overwrite: ${OUTPUT}" >&2
        exit 3
    }
    if [ "${checkpoint_step}" -ge "${TARGET_STEPS}" ]; then
        echo "${OUTPUT} is already complete at step ${checkpoint_step}."
        exit 0
    fi
    echo "Resuming ${VARIANT} from atomic checkpoint step ${checkpoint_step} toward ${TARGET_STEPS}."
    resume_flag=(--resume=true)
else
    echo "Starting fresh ${VARIANT} trajectory toward ${TARGET_STEPS}."
fi

mapfile -t base_args <"${BASE_FILE}"

cmd=(
    srun --kill-on-bad-exit=1 "${EXE}"
    "${base_args[@]}"
    --initial_structure=myosin_x_noise
    --n_actins=1280 --n_myosins=640
    --Lx=7.72 --Ly=1.12 --Lz=1.12
    --am_overlap_target_ratio=0.40
    --am_bidirectional_crossbridge
    --aa_axial_register_ratio="${REGISTER_RATIO}"
    --dt=1e-6 --nsteps="${TARGET_STEPS}" --save_every=1000
    --pre_phase_steps=5000 --interaction_ramp_steps=20000
    --seed=2072311
    --filename="${OUTPUT}"
    "${resume_flag[@]}"
)

printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
exec "${cmd[@]}"
