#!/bin/bash
#SBATCH -p owners,normal,hns
#SBATCH -C "CPU_GEN:GEN|CPU_GEN:MLN"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
# Do NOT let Slurm auto-requeue this job on preemption (owners is PreemptMode=REQUEUE).
# The worker runs with --resume=false and calls new_file(), so a requeue would restart from
# step 0 and lose all progress. Instead, a preempted worker is cancelled and the babysitter
# resumes it from the last checkpoint via sarcomere_resume_new.sbatch (--resume=true).
#SBATCH --no-requeue

set -euo pipefail

# Retrieve simulation arguments passed by the dispatcher.
if [ "$#" -ne 21 ]; then
    echo "Usage: sbatch sarcomere_worker.sbatch n_actins n_myosins lifetime_coeff Lx Ly Lz nsteps seed output_filename filename myosin_radius am_cutoff am_optimal save_every base_lifetime k_mm k_bundle_max max_myosin_bonds bundle_ramp_steps periodic_mask initial_structure" >&2
    exit 2
fi

n_actins=$1
n_myosins=$2
lifetime_coeff=$3
Lx=$4
Ly=$5
Lz=$6
nsteps=$7
SEED=$8
OUTPUT_FILENAME=$9
FILENAME=${10}
MYOSIN_RADIUS=${11}
AM_CUTOFF=${12}
AM_OPTIMAL=${13}
SAVE_EVERY=${14}
BASE_LIFETIME=${15}
K_MM=${16}
K_BUNDLE_MAX=${17}
MAX_MYOSIN_BONDS=${18}
BUNDLE_RAMP_STEPS=${19}
PERIODIC_MASK=${20}
INITIAL_STRUCTURE=${21}

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

echo "Starting simulation for: ${FILENAME}"
echo "Using Seed: ${SEED}"
echo "Parameters: n_actins=${n_actins}, n_myosins=${n_myosins}, lifetime_coeff=${lifetime_coeff}, box=${Lx}x${Ly}x${Lz}, nsteps=${nsteps}, save_every=${SAVE_EVERY}, base_lifetime=${BASE_LIFETIME}, k_mm=${K_MM}, k_bundle_max=${K_BUNDLE_MAX}, max_myosin_bonds=${MAX_MYOSIN_BONDS}, bundle_ramp_steps=${BUNDLE_RAMP_STEPS}, periodic_mask=${PERIODIC_MASK}, initial_structure=${INITIAL_STRUCTURE}"

# Run the simulation
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
    --initial_structure="${INITIAL_STRUCTURE}" \
    --save_every=${SAVE_EVERY} \
    --base_lifetime=${BASE_LIFETIME} \
    --k_mm=${K_MM} \
    --k_bundle_max=${K_BUNDLE_MAX} \
    --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
    --bundle_ramp_steps=${BUNDLE_RAMP_STEPS} \
    --seed=${SEED}

echo "Simulation complete for: ${FILENAME}"
