#!/bin/bash
#SBATCH -p rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH --job-name=catch_ko
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_%a_ko.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_%a_ko.out
#SBATCH --array=0-9
# -p rotskoff only (NOT owners): short data0430-scale runs, no checkpoint infra,
# so we avoid preemption/requeue complexity -- they run straight to completion.

# ===========================================================================
# CATCH-BOND KNOCKOUT / INTERVENTION  (mechanism vs description)
#
# The WSINDy review's core point: regression on stationary trajectories cannot
# separate mechanism from description. The one intervention oracle is the particle
# model itself. Here we KNOCK OUT the load-dependence of the actin-actin catch
# bond and ask: does assembly survive?
#
# The lifetime law (sarcomere.cpp:1356) decomposes cleanly:
#     lifetime  = base_lifetime  +  lifetime_coeff * f_load
#                 \___constant___/    \___load-dependent (catch)___/
#     k_off = dt / lifetime,   f_load = |cos(theta)| * min(f_load_i, f_load_j)
#
# data0430 operating point = pure catch: base_lifetime=0, lifetime_coeff=0.01
# (lifetime = 0.01*f_load; nonzero because f_load>0 while bonded).
#
# INTERVENTION = term swap holding *mean* bond strength fixed:
#     lifetime_coeff=0, base_lifetime=tau0  ->  CONSTANT lifetime, ZERO load-dep.
# Same average bond, but no selective stabilization of high-load/well-oriented
# bonds -- the essence of "catch". Bond FORMATION geometry (antiparallel gating,
# sarcomere.cpp:1280) is untouched; only the load-dependent LIFETIME is removed.
#
# WHY A BRACKET, not one matched control: the catch-bond lifetime-at-break
# distribution (from lt1e-2 logs) is mean~0.0025, median~0.0011, p90~0.0074 --
# the assembly-forming bonds live in the HIGH-LOAD TAIL. A single constant
# lifetime can't reproduce that selectivity, so we bracket tau0 from the mean up
# past the tail. If assembly fails across the WHOLE bracket -> catch is NECESSARY
# (mechanism). If it succeeds at the high end -> it was only "enough total
# binding" (description). Arm 0 is the catch reference (positive control).
#
# 5 arms x 2 reps = 10 tasks. All params matched to data0430 A1500 lt1e-2 EXCEPT
# the two lifetime knobs. Seed depends on rep only (same IC across arms = paired).
# ===========================================================================

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp

# ---- arm table: (tag, lifetime_coeff, base_lifetime) -----------------------
TAG_LIST=(   catch     const0p001  const0p0025  const0p007  const0p02 )
COEFF_LIST=( 0.01      0.0         0.0          0.0         0.0       )
BASE_LIST=(  0.0       0.001       0.0025       0.007       0.02      )
NARMS=5
NREPS=2

arm_idx=$(( SLURM_ARRAY_TASK_ID / NREPS ))
rep_idx=$(( SLURM_ARRAY_TASK_ID % NREPS ))
TAG=${TAG_LIST[$arm_idx]}
LIFETIME_COEFF=${COEFF_LIST[$arm_idx]}
BASE_LIFETIME=${BASE_LIST[$arm_idx]}
rep_label=$(printf "r%02d" "${rep_idx}")

# ---- box + physics matched to data0430 A1500 lt1e-2 ------------------------
Lx=12; Ly=2; Lz=2
n_actins=1500
n_myosins=1000
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
MAX_MYOSIN_BONDS=6
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
SAVE_EVERY=100                   # data0430 cadence: dt_per_frame = 1e-3
nsteps=300000                    # 3000 frames -> well past the ~2000-frame plateau

DATE="intervention"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
OUTPUT_DIR="${SCRATCH_DIR}/data_${DATE}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

FILENAME="A${n_actins}_M${n_myosins}_${TAG}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"

# same IC/noise seed across arms for a given rep (paired A/B; only bond law differs)
SEED=$(( 80000 + rep_idx + 1 ))

if [ -f "${OUTPUT_FILENAME}" ]; then
    echo "Output ${OUTPUT_FILENAME} exists. Skipping."
else
    echo "Running ${FILENAME}  (coeff=${LIFETIME_COEFF}, base=${BASE_LIFETIME}, seed=${SEED})..."
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
        --lifetime_coeff=${LIFETIME_COEFF} \
        --base_lifetime=${BASE_LIFETIME} \
        --initial_structure=myosin_x_noise \
        --save_every=${SAVE_EVERY} \
        --k_mm=${K_MM} \
        --k_bundle_max=${K_BUNDLE_MAX} \
        --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
        --bundle_ramp_steps=${BUNDLE_RAMP_STEPS} \
        --seed=${SEED} \
        > "${LOG_DIR}/${FILENAME}.sim.out"
fi
