#!/bin/bash
#SBATCH -p normal,hns,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=wb_babycron
#SBATCH -o /scratch/users/xsli/sarcomere/logs/wb_babycron_%j.out

# Self-chaining watchdog for the wider-box run (Sherlock-safe: no sleep; one scan
# then reschedule via --begin=now+INTERVAL). Scoped to wider-box only so it is
# independent of the sarco_resume sweeps. Order matters: scan (which submits any
# needed resumes) BEFORE counting, so a just-preempted file's resume is already
# in the queue when we decide whether to stop.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR}"                 # cpp/
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere            # h5py for the babysitter's readers
unset PYTHONPATH

INTERVAL_MIN="${INTERVAL_MIN:-30}"

echo "wb-babycron scan $(date)"
bash wider_box_babysitter.sh || echo "  (babysitter returned nonzero; continuing)"

# Count wider-box jobs: pending/just-started main tasks (name widebox_reg) OR
# running mains + resumes (renamed to their A13500_M9000_... filename).
REMAIN=$(squeue --me -h -o "%j" | grep -cE 'widebox_reg|^A13500_M9000_' || true)
if [ "${REMAIN}" -gt 0 ]; then
    echo "${REMAIN} wider-box job(s) remain -> next scan in ${INTERVAL_MIN} min."
    sbatch --begin="now+${INTERVAL_MIN}minutes" wb_babycron.sbatch
else
    echo "no wider-box jobs remain -> babysitter chain stopping."
fi
