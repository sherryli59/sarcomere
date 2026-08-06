#!/bin/bash

#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mincpus=16
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_2d_sweep_additional/%j_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_2d_sweep_additional/%j_%a.out
#SBATCH --job-name=sarcomere_sweep
#SBATCH --open-mode=append
#SBATCH --time=48:00:00
#SBATCH --array=0-311

# 8 actins * 13 lifetimes * 3 reps = 312 tasks
NREPS=3
N_LIFETIMES=13

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

# Fixed simulation parameters:
Lx=12
Ly=2
Lz=2
nsteps=400000

# Bash arrays must be space-separated
n_actins_list=(450 750 1050 1350 1650 1800 1950 2250)
n_myosins_list=(1000)
lifetime_log_list=(-5.5 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0)

task_id=${SLURM_ARRAY_TASK_ID}

# Unpack 1D task_id into 3D indices
rep_idx=$(( task_id % NREPS ))
remainder=$(( task_id / NREPS ))
lifetime_idx=$(( remainder % N_LIFETIMES ))
actin_idx=$(( remainder / N_LIFETIMES ))

# map lifetime_idx to lifetime_coeff
if [ "${lifetime_idx}" -eq 0 ]; then
    lifetime_coeff=0
    lifetime_label="0"
else
    log10_val=${lifetime_log_list[$((lifetime_idx - 1))]}
    lifetime_coeff=$(python - << EOF
import math
x = float("${log10_val}")
print(f"{10**x:.6g}")
EOF
)
    lifetime_label="1e${log10_val}"
fi

# Extract the number of actins dynamically based on actin_idx
n_actins=${n_actins_list[$actin_idx]}
n_myosins=${n_myosins_list[0]}

# Static Directory Naming
SWEEP_NAME="2d_sweep_additional"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
LOG_DIR="${SCRATCH_DIR}/logs_${SWEEP_NAME}"
OUTPUT_DIR="${SCRATCH_DIR}/data_${SWEEP_NAME}"
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}"

# include replicate in names
rep_label=$(printf "r%02d" "${rep_idx}")

FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"

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

# deterministic, unique seed per (actin, lifetime, rep)
SEED=$(( 10000 * actin_idx + 100 * lifetime_idx + rep_idx + 1 ))

# Check if the output .h5 file already exists
if [ -f "${OUTPUT_FILENAME}" ]; then
    echo "Output file ${OUTPUT_FILENAME} already exists. Skipping simulation."
else
    echo "Running simulation for ${FILENAME}..."
    
    # Save the param script directly into the data folder alongside the .h5
    PARAM_SCRIPT="${OUTPUT_DIR}/${FILENAME}_params.sh"
    cat << EOF > "${PARAM_SCRIPT}"
#!/bin/bash

# Simulation Parameters:
n_actins=${n_actins}
n_myosins=${n_myosins}
lifetime_coeff=${lifetime_coeff}
Lx=${Lx}
Ly=${Ly}
Lz=${Lz}
nsteps=${nsteps}
rep_idx=${rep_idx}
seed=${SEED}

# Physics Parameters:
MYOSIN_RADIUS=${MYOSIN_RADIUS}
AM_CUTOFF=${AM_CUTOFF}
AM_OPTIMAL=${AM_OPTIMAL}
SAVE_EVERY=${SAVE_EVERY}
BASE_LIFETIME=${BASE_LIFETIME}
K_MM=${K_MM}
K_BUNDLE_MAX=${K_BUNDLE_MAX}
MAX_MYOSIN_BONDS=${MAX_MYOSIN_BONDS}
BUNDLE_RAMP_STEPS=${BUNDLE_RAMP_STEPS}
PERIODIC_MASK=${PERIODIC_MASK}

OUTPUT_FILENAME="${OUTPUT_FILENAME}"

srun build/sarcomere \\
    --myosin_radius=${MYOSIN_RADIUS} \\
    --am_cutoff=${AM_CUTOFF} \\
    --am_optimal=${AM_OPTIMAL} \\
    --n_actins=${n_actins} \\
    --n_myosins=${n_myosins} \\
    --filename="${OUTPUT_FILENAME}" \\
    --periodic_mask=${PERIODIC_MASK} \\
    --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \\
    --nsteps=${nsteps} \\
    --lifetime_coeff=${lifetime_coeff} \\
    --initial_structure=myosin_x_noise \\
    --save_every=${SAVE_EVERY} \\
    --base_lifetime=${BASE_LIFETIME} \\
    --k_mm=${K_MM} \\
    --k_bundle_max=${K_BUNDLE_MAX} \\
    --max_myosin_bonds=${MAX_MYOSIN_BONDS} \\
    --bundle_ramp_steps=${BUNDLE_RAMP_STEPS} \\
    --seed=${SEED}
EOF
    chmod +x "${PARAM_SCRIPT}"

    # Execute
    srun build/sarcomere \
        --myosin_radius=${MYOSIN_RADIUS} \
        --am_cutoff=${AM_CUTOFF} \
        --am_optimal=${AM_OPTIMAL} \
        --n_actins=${n_actins} \
        --n_myosins=${n_myosins} \
        --filename="${OUTPUT_FILENAME}" \
        --periodic_mask=${PERIODIC_MASK} \
        --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
        --nsteps=${nsteps} \
        --lifetime_coeff=${lifetime_coeff} \
        --initial_structure=myosin_x_noise \
        --save_every=${SAVE_EVERY} \
        --base_lifetime=${BASE_LIFETIME} \
        --k_mm=${K_MM} \
        --k_bundle_max=${K_BUNDLE_MAX} \
        --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
        --bundle_ramp_steps=${BUNDLE_RAMP_STEPS} \
        --seed=${SEED} \
        > "${LOG_DIR}/${FILENAME}.sim.out"
fi