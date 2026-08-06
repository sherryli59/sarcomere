#!/bin/bash

# ==========================================
# Dispatcher Script
# Run with: bash submit_sweep.sh
# ==========================================

NREPS=1
n_actins_list=(2100 1800 1500 1200 900 600)
lifetime_log_list=(-5.5 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0)
n_myosins=1000

Lx=12
Ly=2
Lz=2
nsteps=100000

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

DATE="0528"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
LOG_DIR="${SCRATCH_DIR}/logs${DATE}"
OUTPUT_DIR="${SCRATCH_DIR}/data${DATE}"
FRAME_BASE_DIR="${SCRATCH_DIR}/frames${DATE}"
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}" "${FRAME_BASE_DIR}"

echo "Starting sweep dispatcher..."

actin_idx=0
for n_actins in "${n_actins_list[@]}"; do
    
    for lifetime_idx in {0..12}; do
        
        # Calculate lifetime coefficient
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
            
            # Deterministic, unique seed per (actin, lifetime, rep)
            SEED=$(( 10000 * actin_idx + 100 * lifetime_idx + rep_idx + 1 ))

            # Check if the output .h5 file already exists
            if [ -f "${OUTPUT_FILENAME}" ]; then
                echo "[Skipping] File already exists: ${FILENAME}"
            else
                mkdir -p "${FRAME_DIR}"
                
                # Write parameters script
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

                echo "[Submitting] Job for: ${FILENAME}"
                
                # Submit the worker job to SLURM, passing variables as arguments
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
    
    actin_idx=$((actin_idx + 1))
done # end actin loop

echo "Dispatch complete!"
