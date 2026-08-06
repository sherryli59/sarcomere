#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mincpus=16
#SBATCH --exclusive
#SBATCH --job-name=sarco_resume
#SBATCH --open-mode=append
#SBATCH --time=48:00:00
# NOTE: --partition, --error, and --output are now provided dynamically by the babysitter!

# Arguments passed by the babysitter
RESUME_LIST=$1
LOG_DIR=$2
TARGET_TOTAL_STEPS=$3

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

mkdir -p "${LOG_DIR}"

# 1. Select the specific file for this array task from the text file
OUTPUT_FILENAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${RESUME_LIST}")
FILENAME=$(basename "${OUTPUT_FILENAME}" .h5)

echo "[${FILENAME}] Selected for resume task ${SLURM_ARRAY_TASK_ID}: ${OUTPUT_FILENAME}"

# 2. Dynamically parse the parameters from the filename
if [[ "$FILENAME" =~ ^A([0-9]+)_M([0-9]+)_lt([^_]+)_r([0-9]+)_L([0-9]+)x([0-9]+)x([0-9]+)$ ]]; then
    n_actins="${BASH_REMATCH[1]}"
    n_myosins="${BASH_REMATCH[2]}"
    lifetime_label="${BASH_REMATCH[3]}"
    rep_label="${BASH_REMATCH[4]}"
    Lx="${BASH_REMATCH[5]}"
    Ly="${BASH_REMATCH[6]}"
    Lz="${BASH_REMATCH[7]}"
else
    echo "[${FILENAME}] Does not match expected pattern. Skipping."
    exit 1
fi

# 3. Convert lifetime_label back into lifetime_coeff
if [ "${lifetime_label}" = "0" ]; then
    lifetime_coeff=0
else
    log10_val=${lifetime_label#1e}
    lifetime_coeff=$(python - << EOF
import math
x = float("${log10_val}")
print(f"{10**x:.6g}")
EOF
)
fi

SAVE_EVERY=100

# 4. Extract current step using Python
CURRENT_STEP_RAW=$(python - "${OUTPUT_FILENAME}" "${SAVE_EVERY}" << 'PY'
import sys
fname = sys.argv[1]
save_every = int(sys.argv[2])
step = -1
try:
    import h5py
except Exception:
    sys.exit(0)
try:
    with h5py.File(fname, "r") as f:
        for cand in ("/actin/centers", "actin/centers", "/actin/center", "actin/center", "/myosin/center"):
            if cand in f:
                d = f[cand]
                if hasattr(d, 'shape') and len(d.shape) >= 1:
                    step = int(d.shape[0]) * save_every
                    break
        if step < 0:
            for p in ("/resume/current_step", "/state/current_step"):
                if p in f:
                    data = f[p][...]
                    if data.size > 0:
                        step = int(data.reshape(-1)[-1])
except Exception:
    pass
print(step)
PY
)

CURRENT_STEP=$(echo "$CURRENT_STEP_RAW" | tr -d '[:space:]')

if ! [[ "${CURRENT_STEP}" =~ ^-?[0-9]+$ ]] || [ "${CURRENT_STEP}" -lt 0 ]; then
    echo "[${FILENAME}] Could not parse current_step. Skipping."
    exit 0
fi

remaining=$(( TARGET_TOTAL_STEPS - CURRENT_STEP ))
if [ "${remaining}" -le 0 ]; then
    echo "[${FILENAME}] Complete. Exiting."
    exit 0
fi

run_steps=${remaining}
SEED=$(( n_actins + 1000 * (10#$rep_label) )) 

MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
MAX_MYOSIN_BONDS=6
PERIODIC_MASK=111

echo "[${FILENAME}] Resuming from step ${CURRENT_STEP} for ${run_steps} steps."

# 5. Run the engine 
srun build/sarcomere \
    --myosin_radius=${MYOSIN_RADIUS} \
    --am_cutoff=${AM_CUTOFF} \
    --am_optimal=${AM_OPTIMAL} \
    --n_actins=${n_actins} \
    --n_myosins=${n_myosins} \
    --filename="${OUTPUT_FILENAME}" \
    --periodic_mask=${PERIODIC_MASK} \
    --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
    --nsteps=${TARGET_TOTAL_STEPS} \ 
    --lifetime_coeff=${lifetime_coeff} \
    --save_every=${SAVE_EVERY} \
    --base_lifetime=${BASE_LIFETIME} \
    --k_mm=${K_MM} \
    --k_bundle_max=${K_BUNDLE_MAX} \
    --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
    --bundle_ramp_steps=0 \
    --seed=${SEED} \
    --resume=true \
    >> "${LOG_DIR}/${FILENAME}.sim.out" 2>&1