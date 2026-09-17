#!/bin/bash
#SBATCH --job-name=thread_equiv
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --array=0-71%18
#SBATCH --output=/scratch/users/xsli/sarcomere/thread_equiv/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/thread_equiv/logs/%x_%A_%a.err
#
# Multi-thread equivalence rig.
#
# Bitwise comparison is impossible above one thread: the force reductions run in
# whatever order the team happens to take, so two identical runs decorrelate.
# That rules out the instrument used for the per-pair gates, and it is why a
# barrier restructuring could not be validated so far.
#
# The construction here is a PAIRED NULL.  Thermal noise is drawn serially from
# one GSL stream before the parallel region (langevin.cpp sample_step), so with
# a fixed seed the only source of run-to-run difference is reduction order.
# Therefore two runs of the SAME code (ctrl_a, ctrl_b) and a run of an
# EQUIVALENT variant (treat) are exchangeable: d(ctrl_a, treat) must be drawn
# from the same distribution as d(ctrl_a, ctrl_b).  A real behavioural change
# shows up as a divergence curve outside that null.
#
# Arms pert_* are a deliberate perturbation ladder: actin_steric_k scaled by
# 1+eps.  They are not controls, they CALIBRATE the rig -- the smallest eps the
# rig separates from the null is its detection threshold.  Without them, "no
# difference found" would be uninformative.
set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export SARC_EQUIV_ROOT=${SARC_EQUIV_ROOT:-/scratch/users/xsli/sarcomere/thread_equiv}
mkdir -p "$SARC_EQUIV_ROOT"/{data,logs}

seeds=(20260801 20260802 20260803 20260804 20260805 20260806
       20260807 20260808 20260809 20260810 20260811 20260812)
arms=(ctrl_a ctrl_b treat pert1em8 pert1em6 pert1em4)

n_seeds=${#seeds[@]}
task=${SLURM_ARRAY_TASK_ID}
arm=${arms[$((task / n_seeds))]}
seed=${seeds[$((task % n_seeds))]}

# The treat arm is what a candidate change is wired into.  Left empty it is a
# third control, which is the right default: it proves the rig reports "no
# difference" when there is none.
extra=()
case "$arm" in
    ctrl_a|ctrl_b) ;;
    treat)    read -r -a extra <<< "${SARC_EQUIV_TREAT_ARGS:-}" ;;
    pert1em8) extra=(--actin_steric_k=2000.00002) ;;    # 1e-8 relative
    pert1em6) extra=(--actin_steric_k=2000.002) ;;      # 1e-6 relative
    pert1em4) extra=(--actin_steric_k=2000.2) ;;        # 1e-4 relative
    *) echo "unknown arm $arm" >&2; exit 2 ;;
esac

nsteps=${SARC_EQUIV_NSTEPS:-2000}
save_every=${SARC_EQUIV_SAVE_EVERY:-25}

bash /scratch/users/xsli/sarcomere/thread_equiv_run.sh \
    "${arm}_seed${seed}" "$seed" "$OMP_NUM_THREADS" "$nsteps" "$save_every" \
    ${extra[@]+"${extra[@]}"}
