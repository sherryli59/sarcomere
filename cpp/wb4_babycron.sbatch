#!/bin/bash
#SBATCH -p normal,hns,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=wb4_babycron
#SBATCH -o /scratch/users/xsli/sarcomere/logs/wb4_babycron_%j.out

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

echo "wb4-babycron scan $(date)"
bash wb4_babysitter.sh || echo "  (babysitter returned nonzero; continuing)"

# STOP CONDITION FIX (race bug 2026-07-04): the chain previously stopped when
# "no jobs in queue", but the babysitter SKIPS files modified <30min ago (they might
# still be alive). When the mains die at walltime, that same 30-min window has zero
# queued jobs AND no resume yet submitted -> the old logic declared the run "done" at
# 36%. FIX: stop ONLY when every target file is actually COMPLETE (frames>=TARGET).
# The chain then persists across allocation gaps and the babysitter resumes stalled
# files on the NEXT scan (once they're >30min stale).
TARGET_FRAMES=1500      # 300000 steps / save_every 200
INCOMPLETE=$(python - <<'PY'
import glob, h5py
n=0
for f in glob.glob("/scratch/users/xsli/sarcomere/data_wb4x4/A6000_M4000_*_L12x4x4.h5"):
    try:
        with h5py.File(f,"r") as h:
            if h["actin/center"].shape[0] < 1500: n+=1
    except Exception:
        n+=1     # unreadable (locked/mid-write) -> treat as incomplete, keep chaining
print(n)
PY
)
INCOMPLETE=$(echo "${INCOMPLETE}" | tr -d '[:space:]')
if [ "${INCOMPLETE:-1}" -gt 0 ]; then
    echo "${INCOMPLETE} wider-box file(s) still incomplete -> next scan in ${INTERVAL_MIN} min."
    sbatch --begin="now+${INTERVAL_MIN}minutes" wb4_babycron.sbatch
else
    echo "all wider-box files complete -> babysitter chain stopping."
fi
