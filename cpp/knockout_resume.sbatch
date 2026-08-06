#!/bin/bash
#SBATCH -p rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --time=1-00:00:00
#SBATCH --job-name=ko_resume
#SBATCH --open-mode=append
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_%a_koresume.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_%a_koresume.out
#SBATCH --array=0-9
# Resume the catch-bond knockout runs (they died at the 12h walltime with NO
# babysitter, ~30-53% done). rotskoff is non-preemptable + 1-day walltime, so each
# arm finishes in one allocation (~1500 more frames ~13h). Continues from the HDF5
# /resume checkpoint. Same arm table as catch_knockout.sbatch.

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp

# arm table (tag, lifetime_coeff, base_lifetime) -- MUST match catch_knockout.sbatch
TAG_LIST=(   catch     const0p001  const0p0025  const0p007  const0p02 )
COEFF_LIST=( 0.01      0.0         0.0          0.0         0.0       )
BASE_LIST=(  0.0       0.001       0.0025       0.007       0.02      )
NREPS=2

arm_idx=$(( SLURM_ARRAY_TASK_ID / NREPS ))
rep_idx=$(( SLURM_ARRAY_TASK_ID % NREPS ))
TAG=${TAG_LIST[$arm_idx]}
LIFETIME_COEFF=${COEFF_LIST[$arm_idx]}
BASE_LIFETIME=${BASE_LIST[$arm_idx]}
rep_label=$(printf "r%02d" "${rep_idx}")

# box + physics matched to catch_knockout.sbatch
Lx=12; Ly=2; Lz=2
n_actins=1500; n_myosins=1000
MYOSIN_RADIUS=0.058; AM_CUTOFF=0.08; AM_OPTIMAL=0.06
MAX_MYOSIN_BONDS=6; K_MM=0; K_BUNDLE_MAX=0.1; PERIODIC_MASK=111
SAVE_EVERY=100
nsteps=300000
SEED=$(( 80000 + rep_idx + 1 ))

OUTPUT_DIR="/scratch/users/xsli/sarcomere/data_intervention"
LOG_DIR="/scratch/users/xsli/sarcomere/logs"
FILENAME="A${n_actins}_M${n_myosins}_${TAG}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"

if [ ! -f "${OUTPUT_FILENAME}" ]; then echo "[${FILENAME}] missing; abort"; exit 1; fi

# completeness check: skip if already at target (frames*save_every >= nsteps)
CURRENT_STEP=$(python - "${OUTPUT_FILENAME}" "${SAVE_EVERY}" <<'PY'
import sys, h5py
fn, se = sys.argv[1], int(sys.argv[2]); step=-1
try:
    with h5py.File(fn,"r") as f:
        step = int(f["actin/center"].shape[0]) * se
except Exception: pass
print(step)
PY
)
CURRENT_STEP=$(echo "${CURRENT_STEP}" | tr -d '[:space:]')
if [[ "${CURRENT_STEP}" =~ ^[0-9]+$ ]] && [ "${CURRENT_STEP}" -ge "${nsteps}" ]; then
    echo "[${FILENAME}] already complete (${CURRENT_STEP}); exit"; exit 0
fi

echo "Resuming ${FILENAME} from step ${CURRENT_STEP}/${nsteps} (coeff=${LIFETIME_COEFF}, base=${BASE_LIFETIME})..."
srun build/sarcomere \
    --myosin_radius=${MYOSIN_RADIUS} --am_cutoff=${AM_CUTOFF} --am_optimal=${AM_OPTIMAL} \
    --n_actins=${n_actins} --n_myosins=${n_myosins} \
    --filename="${OUTPUT_FILENAME}" --periodic_mask=${PERIODIC_MASK} \
    --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
    --nsteps=${nsteps} --lifetime_coeff=${LIFETIME_COEFF} --base_lifetime=${BASE_LIFETIME} \
    --save_every=${SAVE_EVERY} --k_mm=${K_MM} --k_bundle_max=${K_BUNDLE_MAX} \
    --max_myosin_bonds=${MAX_MYOSIN_BONDS} --bundle_ramp_steps=0 --resume=true \
    --seed=${SEED} \
    >> "${LOG_DIR}/${FILENAME}.sim.out" 2>&1
