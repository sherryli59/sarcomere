#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=widebox_reg
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_%a.out
#SBATCH --array=0-3
# -p rotskoff,owners: runs on whichever frees first. owners is PREEMPTABLE, so
# --no-requeue is set (Slurm will NOT auto-requeue on preemption); the dedicated
# babysitter (wider_box_babycron.sbatch) detects the stall and resumes from the
# HDF5 checkpoint instead. Worker names each task after its output file so the
# babysitter can match it.
# 2 max_myosin_bonds values x 2 replicas = 4 tasks. (Add --qos=long + bump
# --time for >2 day runs; extend the arrays below for more replicas.)

# ---------------------------------------------------------------------------
# Wider-cross-section registration / mechanism test.
#
# WHY: data0430 is 12x2x2 -- too thin to resolve transverse coherence (R floors
# ~0.3, finite-size-inflated; can't tell registered crystal from dephased
# smectic). This widens the cross-section to 6x6 (9x area, density-matched) so
# transverse_coherence.py can map the dephasing length xi_perp, and contrasts
# max_myosin_bonds = 6 (full lateral multi-binding) vs 2 (minimal functional,
# one antiparallel pair per myosin) -- both assemble; isolates lateral coupling.
# (NOT vs 1, which disables pair-bridging and would prevent assembly.) The CORRECT
# registration-mechanism knob (k_bundle is mere transverse confinement).
#
# COST CAVEATS (review before submitting):
#   - ~9x the particles of an A1500 run -> ~9x slower/step + bigger system to
#     equilibrate. Budget days; consider a 1-task pilot first.
#   - File size ~9x: with save_every=200 (-> ~2000 frames) each .h5 is ~60GB;
#     4 tasks ~ 240GB on $SCRATCH. Raise save_every to shrink.
# ---------------------------------------------------------------------------

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp

# box: same axial length, 3x wider transverse (9x cross-sectional area)
Lx=12; Ly=6; Lz=6
VOL_SCALE=9                      # (6*6)/(2*2)

# density-matched to a well-ordered 12x2x2 operating point (A1500, lt1e-2)
n_actins=$((1500 * VOL_SCALE))   # 13500
n_myosins=$((1000 * VOL_SCALE))  # 9000
lifetime_coeff=0.01              # lt1e-2: ordered, away from lt1e0 trapping
lifetime_label="1e-2"

# mechanism contrast: full lateral multi-binding (6, as data0430) vs
# minimal-FUNCTIONAL multi-binding (2 = one actin per side = bridges one
# antiparallel pair). NOT mmb=1: that binds a single actin, cannot bridge a pair
# (code splits binding into 2 sides, sarcomere.cpp:850), so the motor can't
# contract and nothing assembles -- which would confound the registration test.
MMB_LIST=(6 2)
NREPS=2

mmb_idx=$(( SLURM_ARRAY_TASK_ID / NREPS ))
rep_idx=$(( SLURM_ARRAY_TASK_ID % NREPS ))
MAX_MYOSIN_BONDS=${MMB_LIST[$mmb_idx]}
rep_label=$(printf "r%02d" "${rep_idx}")

# fixed params matched to data0430
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=200                   # ~2000 frames at nsteps=400000 (file ~60GB)
BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
nsteps=400000                    # may need more for the larger box to settle

DATE="widebox"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
OUTPUT_DIR="${SCRATCH_DIR}/data_${DATE}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_mmb${MAX_MYOSIN_BONDS}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"

# name this array task after its output file so the babysitter's "is the worker
# still alive?" guard (squeue -n FILENAME) can match it.
scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true

# unique seed per (mmb, rep)
SEED=$(( 70000 + 100 * mmb_idx + rep_idx + 1 ))

if [ -f "${OUTPUT_FILENAME}" ]; then
    echo "Output ${OUTPUT_FILENAME} exists. Skipping."
else
    echo "Running ${FILENAME}  (mmb=${MAX_MYOSIN_BONDS}, seed=${SEED})..."
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
