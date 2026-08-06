#!/bin/bash

# Configuration
DATA_DIR=${1:-"/scratch/users/xsli/sarcomere/data0416"}
TARGET_TOTAL_STEPS=400000
SAVE_EVERY=100
SMOKE_TEST_STEPS=500  # Number of steps to run to verify it doesn't crash
SMOKE_TEST_DIR="./smoke_test_run"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

echo "Scanning $DATA_DIR for stalled/incomplete .h5 files..."

FOUND_FILE=""
CURRENT_STEP=0

# 1. Find the first file that needs resuming
for FILE in "${DATA_DIR}"/*.h5; do
    [ -e "$FILE" ] || continue

    CURRENT_STEP_RAW=$(python - "${FILE}" "${SAVE_EVERY}" << 'PY'
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
    STEP=$(echo "$CURRENT_STEP_RAW" | tr -d '[:space:]')

    if [[ "${STEP}" =~ ^-?[0-9]+$ ]] && [ "${STEP}" -ge 0 ] && [ "${STEP}" -lt "${TARGET_TOTAL_STEPS}" ]; then
        FOUND_FILE="$FILE"
        CURRENT_STEP="$STEP"
        break
    fi
done

if [ -z "$FOUND_FILE" ]; then
    echo "No incomplete files found in ${DATA_DIR}. Everything is finished!"
    exit 0
fi

echo "Found incomplete file: $FOUND_FILE (Current step: $CURRENT_STEP)"

# 2. Setup the isolated test directory and copy the file
mkdir -p "${SMOKE_TEST_DIR}"
BASENAME=$(basename "${FOUND_FILE}")
TEST_FILE="${SMOKE_TEST_DIR}/${BASENAME}"
FILENAME="${BASENAME%.h5}"

echo "Copying to isolated test environment: ${TEST_FILE}"
cp "${FOUND_FILE}" "${TEST_FILE}"

# 3. Dynamically parse the parameters from the filename
if [[ "$FILENAME" =~ ^A([0-9]+)_M([0-9]+)_lt([^_]+)_r([0-9]+)_L([0-9]+)x([0-9]+)x([0-9]+)$ ]]; then
    n_actins="${BASH_REMATCH[1]}"
    n_myosins="${BASH_REMATCH[2]}"
    lifetime_label="${BASH_REMATCH[3]}"
    rep_label="${BASH_REMATCH[4]}"
    Lx="${BASH_REMATCH[5]}"
    Ly="${BASH_REMATCH[6]}"
    Lz="${BASH_REMATCH[7]}"
else
    echo "[!] Filename ${FILENAME} does not match expected pattern. Aborting."
    exit 1
fi

# 4. Convert lifetime_label back into lifetime_coeff
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

# 5. Define run constants
SEED=$(( n_actins + 1000 * (10#$rep_label) )) 
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
MAX_MYOSIN_BONDS=6
PERIODIC_MASK=111

echo "======================================================"
echo "Starting local resume test on ${BASENAME}"
echo "Running for ${SMOKE_TEST_STEPS} steps to verify execution..."
echo "======================================================"

# 6. Run the engine directly
build/sarcomere \
    --myosin_radius=${MYOSIN_RADIUS} \
    --am_cutoff=${AM_CUTOFF} \
    --am_optimal=${AM_OPTIMAL} \
    --n_actins=${n_actins} \
    --n_myosins=${n_myosins} \
    --filename="${TEST_FILE}" \
    --periodic_mask=${PERIODIC_MASK} \
    --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
    --nsteps=${SMOKE_TEST_STEPS} \
    --lifetime_coeff=${lifetime_coeff} \
    --save_every=${SAVE_EVERY} \
    --base_lifetime=${BASE_LIFETIME} \
    --k_mm=${K_MM} \
    --k_bundle_max=${K_BUNDLE_MAX} \
    --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
    --bundle_ramp_steps=0 \
    --seed=${SEED} \
    --resume=true 

echo "======================================================"
echo "Smoke test complete. Check ${TEST_FILE} for appended data."