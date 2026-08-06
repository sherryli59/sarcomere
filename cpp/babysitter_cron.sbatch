#!/bin/bash
#SBATCH -p normal,hns,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=babycron_myo1800
#SBATCH -o /scratch/users/xsli/sarcomere/logs0624_myo1800/babycron_%j.out

# Self-chaining watchdog (Sherlock-safe: no sleep, exits immediately after one scan).
# Each run: scan the sweep dir once via babysitter.sh, then schedule the next run with
# --begin=now+INTERVAL so the wait happens in the Slurm queue, not in a held allocation.
# The chain stops itself once no sweep/resume jobs remain.
#
# TARGET_TOTAL_STEPS MUST equal the dispatcher's nsteps (400000).

set -uo pipefail
cd "${SLURM_SUBMIT_DIR}"                 # cpp/  (so babysitter's relative paths resolve)

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere            # h5py for babysitter's stall-detection reader
unset PYTHONPATH

export TARGET_TOTAL_STEPS=400000
DATA_DIR="/scratch/users/xsli/sarcomere/data0624_myo1800"
INTERVAL_MIN="${INTERVAL_MIN:-30}"

echo "babycron scan $(date)"
bash babysitter.sh "${DATA_DIR}" || echo "  (babysitter scan returned nonzero; continuing)"

# Count my jobs still in the system, EXCLUDING this watchdog chain (job name babycron_myo1800).
# That leaves the sweep jobs (A1800_M*) and any resume-array jobs babysitter just launched.
OTHERS=$(squeue --me -h -o "%j" | grep -vc '^babycron_myo1800' || true)
if [ "${OTHERS}" -gt 0 ]; then
    echo "${OTHERS} sweep/resume jobs remain -> next scan in ${INTERVAL_MIN} min."
    sbatch --begin="now+${INTERVAL_MIN}minutes" babysitter_cron.sbatch
else
    echo "No sweep/resume jobs remain -> babysitter chain stopping."
fi
