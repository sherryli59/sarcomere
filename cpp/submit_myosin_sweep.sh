#!/bin/bash

# =====================================================================
# Fixed-density MYOSIN-axis dispatcher (collapse / partial-correlation test)
# Run from cpp/ with:  bash submit_myosin_sweep.sh
#
# Goal: hold actin density fixed and sweep myosin number to DECOUPLE motor
# occupancy from density. Combined with the existing 6-density nM=1000 sweep
# (data0430), this lets analyze_particle_sims.py --collapse-test run the
# partial correlations  partial(tau_peak, density | occupancy)  vs
# partial(tau_peak, occupancy | density)  to identify the hump's hidden variable.
#
# Consistency with data0430 (required so tau_peak is comparable):
#   nsteps=400000, save_every=100  -> ~4001 frames (data0430 has 4007).
#   This also matches the babysitter default TARGET_TOTAL_STEPS=400000, so its
#   completion test (current_step < TARGET) is correct for these runs.
# =====================================================================

NREPS=3
n_actins=1800                                   # FIXED density (strong hump, headroom for occupancy)
n_myosins_list=(500 750 1000 1500 2000)         # the occupancy axis
lifetime_log_list=(-5.5 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0)

Lx=12
Ly=2
Lz=2
nsteps=400000                                   # matches data0430 (NOT 100000)

MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=100
BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
MAX_MYOSIN_BONDS=6
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
INITIAL_STRUCTURE=myosin_x_noise

DATE="0624_myo1800"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
LOG_DIR="${SCRATCH_DIR}/logs${DATE}"
OUTPUT_DIR="${SCRATCH_DIR}/data${DATE}"
FRAME_BASE_DIR="${SCRATCH_DIR}/frames${DATE}"
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}" "${FRAME_BASE_DIR}"

# Always submit from cpp/ so the worker's relative `build/sarcomere` path resolves.
cd "$(dirname "$0")"

echo "Starting MYOSIN-axis sweep dispatcher (nA=${n_actins} fixed)..."
echo "Output: ${OUTPUT_DIR}"

myo_idx=0
for n_myosins in "${n_myosins_list[@]}"; do

    for lifetime_idx in {0..12}; do

        # Calculate lifetime coefficient (idx 0 == no catch bond)
        if [ "${lifetime_idx}" -eq 0 ]; then
            lifetime_coeff=0
            lifetime_label="0"
        else
            log10_val=${lifetime_log_list[$((lifetime_idx - 1))]}
            lifetime_coeff=$(python - << EOF
import math
print(f"{10**float('${log10_val}'):.6g}")
EOF
)
            lifetime_label="1e${log10_val}"
        fi

        for rep_idx in $(seq 0 $((NREPS - 1))); do

            rep_label=$(printf "r%02d" "${rep_idx}")
            FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
            OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
            FRAME_DIR="${FRAME_BASE_DIR}/${FILENAME}"

            # Deterministic, unique seed, DISJOINT from the original sweep's range.
            SEED=$(( 2000000 + 100000 * myo_idx + 1000 * lifetime_idx + rep_idx + 1 ))

            if [ -f "${OUTPUT_FILENAME}" ]; then
                echo "[Skipping] File already exists: ${FILENAME}"
            else
                mkdir -p "${FRAME_DIR}"

                PARAM_SCRIPT="${FRAME_DIR}/run_parameters.sh"
                cat << EOF > "${PARAM_SCRIPT}"
#!/bin/bash
n_actins=${n_actins}
n_myosins=${n_myosins}
lifetime_coeff=${lifetime_coeff}
Lx=${Lx}
Ly=${Ly}
Lz=${Lz}
nsteps=${nsteps}
rep_idx=${rep_idx}
seed=${SEED}
myosin_radius=${MYOSIN_RADIUS}
am_cutoff=${AM_CUTOFF}
am_optimal=${AM_OPTIMAL}
save_every=${SAVE_EVERY}
base_lifetime=${BASE_LIFETIME}
k_mm=${K_MM}
k_bundle_max=${K_BUNDLE_MAX}
max_myosin_bonds=${MAX_MYOSIN_BONDS}
bundle_ramp_steps=${BUNDLE_RAMP_STEPS}
periodic_mask=${PERIODIC_MASK}
initial_structure=${INITIAL_STRUCTURE}
output_filename=${OUTPUT_FILENAME}
EOF
                chmod +x "${PARAM_SCRIPT}"

                echo "[Submitting] ${FILENAME}  (seed ${SEED})"

                sbatch \
                    --job-name="${FILENAME}" \
                    --output="${LOG_DIR}/${FILENAME}_%j.out" \
                    --error="${LOG_DIR}/${FILENAME}_%j.err" \
                    sarcomere_worker.sbatch \
                    "${n_actins}" \
                    "${n_myosins}" \
                    "${lifetime_coeff}" \
                    "${Lx}" \
                    "${Ly}" \
                    "${Lz}" \
                    "${nsteps}" \
                    "${SEED}" \
                    "${OUTPUT_FILENAME}" \
                    "${FILENAME}" \
                    "${MYOSIN_RADIUS}" \
                    "${AM_CUTOFF}" \
                    "${AM_OPTIMAL}" \
                    "${SAVE_EVERY}" \
                    "${BASE_LIFETIME}" \
                    "${K_MM}" \
                    "${K_BUNDLE_MAX}" \
                    "${MAX_MYOSIN_BONDS}" \
                    "${BUNDLE_RAMP_STEPS}" \
                    "${PERIODIC_MASK}" \
                    "${INITIAL_STRUCTURE}"
            fi

        done # end rep loop
    done # end lifetime loop

    myo_idx=$((myo_idx + 1))
done # end myosin loop

echo "Dispatch complete!  (expected $(( ${#n_myosins_list[@]} * 13 * NREPS )) jobs)"
