#!/bin/bash

#SBATCH -p hns 
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=16
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_%a.out
#SBATCH --job-name=sarcomere_sweep
#SBATCH --open-mode=append
#SBATCH --time=16:00:00
#SBATCH --array=0

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

# Fixed simulation parameters:
Lx=12
Ly=2
Lz=2
nsteps=200000

# Parameter arrays:
n_actins_list=(2200)
n_myosins_list=(1000)
# log10 values for lifetime_coeff (excluding 0)
lifetime_log_list=(-3)


task_id=${SLURM_ARRAY_TASK_ID}

if [ "${task_id}" -eq 9 ]; then
    lifetime_coeff=0
    lifetime_label="0"
else
    log10_val=${lifetime_log_list[$task_id]}
    lifetime_coeff=$(python - << EOF
import math
x = float("${log10_val}")
print(f"{10**x:.6g}")
EOF
)
lifetime_label="1e${log10_val}"
fi


n_actins=${n_actins_list[0]}  # Use the first value from the list
n_myosins=${n_myosins_list[0]}  # Use the first value from the list

# Define date string (MMDD) and attach to directory names:
DATE=$(date +%m%d)
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
LOG_DIR="${SCRATCH_DIR}/logs${DATE}"
OUTPUT_DIR="${SCRATCH_DIR}/data${DATE}"
FRAME_BASE_DIR="${SCRATCH_DIR}/frames${DATE}"
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}" "${FRAME_BASE_DIR}"

# Construct file and directory names:
FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
FRAME_DIR="${FRAME_BASE_DIR}/${FILENAME}"
mkdir -p "${FRAME_DIR}"

# Define simulation parameters as variables
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

# Optionally, create a parameter script that records these settings:
PARAM_SCRIPT="${FRAME_DIR}/run_parameters.sh"
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

# Define simulation parameters as variables
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

# Output Files:
OUTPUT_FILENAME="${OUTPUT_FILENAME}"
FRAME_DIR="${FRAME_DIR}"

# Command Used:
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
    --bundle_ramp_steps=${BUNDLE_RAMP_STEPS}
EOF
chmod +x "${PARAM_SCRIPT}"

# Run the main simulation:
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
    --seed=1 \
	--resume \
    > "${LOG_DIR}/${FILENAME}.sim.out"

# # Run the visualization script:
# srun python ../analysis/visualize_traj.py \
#     --filename="${OUTPUT_FILENAME}" \
#     --frame_dir="${FRAME_DIR}" \
#     --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
#     --myosin_radius=0.05 \
#     --sample_every=10 \
#     --start_frame=0 \
#     --myosin_display cb_attached \
#     > "${LOG_DIR}/${FILENAME}.vis.out"
