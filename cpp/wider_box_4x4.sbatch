#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wb4_reg
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_%a.out
#SBATCH --array=0-3
# ===========================================================================
# RESCOPED wider-box run (replaces the INFEASIBLE 6x6 wider_box_registration).
#
# WHY THE 6x6 FAILED: per-step cost scales super-linearly with bond count and
# GROWS as the system assembles. Measured after ~1.3 days: mmb2 ~9 s/step (~40
# days to 2000 frames), mmb6 ~80 s/step (~1 YEAR). 9x particles bought 9-80x the
# per-step cost. Infeasible for reaching assembled steady state.
#
# FIX: 4x4 cross-section = 4x transverse stats vs the 12x2x2 data0430 (enough to
# lift R off the ~0.3 finite-size floor and shrink WSINDy fluctuations) but ~2x
# fewer particles than 6x6 -> tractable. mmb2 ONLY on the wide box (mmb6 is the
# lethal arm; keep the mmb2-vs-6 mechanism contrast at the cheap 2x2 scale).
# 4 replicas (was 2): more reps directly attacks the WSINDy held-out-R^2 variance.
#
# PAIR WITH: a rebuild that strips the hot-loop debug printfs
#   interaction.cpp:462-467 (endstop_force), sarcomere.cpp:703-704 (zero-xlink
#   warning), sarcomere.cpp:1364-1367 (k_off) -- 2.3M prints/run, ~10-30% of step
#   time + multi-GB logs. All pure output (values computed before the print).
# ===========================================================================

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp

# box: same axial length, 2x wider transverse (4x cross-sectional area vs 2x2)
Lx=12; Ly=4; Lz=4
VOL_SCALE=4                      # (4*4)/(2*2)

# density-matched to the A1500 lt1e-2 operating point
n_actins=$((1500 * VOL_SCALE))   # 6000
n_myosins=$((1000 * VOL_SCALE))  # 4000
lifetime_coeff=0.01
lifetime_label="1e-2"

MAX_MYOSIN_BONDS=6               # full multi-binding (as data0430); mmb2-vs-6 contrast stays at 2x2
NREPS=4                          # 4 reps -> better WSINDy held-out statistics
rep_idx=${SLURM_ARRAY_TASK_ID}
rep_label=$(printf "r%02d" "${rep_idx}")

MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=200
BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
nsteps=300000                    # assembly plateau (~200k in data0430) + margin; tune after a step-rate read

DATE="wb4x4"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
OUTPUT_DIR="${SCRATCH_DIR}/data_${DATE}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_mmb${MAX_MYOSIN_BONDS}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"

# name task after its file so a babysitter's "worker alive?" guard can match it.
scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true
SEED=$(( 90000 + rep_idx + 1 ))

if [ -f "${OUTPUT_FILENAME}" ]; then
    echo "Output ${OUTPUT_FILENAME} exists. Skipping."
else
    echo "Running ${FILENAME}  (mmb=${MAX_MYOSIN_BONDS}, seed=${SEED})..."
    srun build/sarcomere \
        --myosin_radius=${MYOSIN_RADIUS} --am_cutoff=${AM_CUTOFF} --am_optimal=${AM_OPTIMAL} \
        --n_actins=${n_actins} --n_myosins=${n_myosins} \
        --filename="${OUTPUT_FILENAME}" --periodic_mask=${PERIODIC_MASK} \
        --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
        --nsteps=${nsteps} --lifetime_coeff=${lifetime_coeff} \
        --initial_structure=myosin_x_noise --save_every=${SAVE_EVERY} \
        --base_lifetime=${BASE_LIFETIME} --k_mm=${K_MM} \
        --k_bundle_max=${K_BUNDLE_MAX} --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
        --bundle_ramp_steps=${BUNDLE_RAMP_STEPS} --seed=${SEED} \
        > "${LOG_DIR}/${FILENAME}.sim.out"
fi
