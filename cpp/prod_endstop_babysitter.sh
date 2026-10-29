#!/bin/bash
# Babysitter for prod_endstop_20260930.  Resubmits any array task whose
# trajectory has not reached TARGET steps and which is not already queued or
# running.  Resubmission is deliberate, so it exports PROD_ALLOW_RESUME=1 to
# let the production sbatch continue a partial file instead of refusing.
#
# Called from prod_endstop_babycron.sbatch, which runs as a JOB and re-submits
# itself -- /etc/agents/AGENTS.md forbids unattended processes on login nodes.
set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
ROOT=/scratch/users/xsli/sarcomere/prod_endstop_20260930
DATA="${ROOT}/data"
TARGET=500000
JOBNAME=prodES
export HDF5_USE_FILE_LOCKING=FALSE

# Must mirror submit_prod_endstop.sbatch exactly.
sigmas=(0.40 0.40 0.40 0.40 0.40 0.40 0.65 0.65 0.65 0.65 0.65 0.65)
kappas=(3    3    3    0    0    0    3    3    3    0    0    0)
seeds=(20260805 20260815 20260825 20260805 20260815 20260825 20260805 20260815 20260825 20260805 20260815 20260825)
tag_for() {  # task -> output tag
    local t=$1
    echo "A1600_M320_Lx6_intermediate_sig${sigmas[$t]}_kaa${kappas[$t]}_seed${seeds[$t]}_5000ms_endstop"
}

# Tasks already queued or running, so we never double-submit.
live=$(squeue --me -h -n "${JOBNAME}" -o "%K" 2>/dev/null \
       | tr ',' '\n' | sed 's/%.*//' \
       | awk '/^[0-9]+$/{print} /^[0-9]+-[0-9]+$/{split($0,a,"-"); for(i=a[1];i<=a[2];i++) print i}' \
       | sort -un | tr '\n' ' ')
echo "live tasks: [${live:-none}]"

progress() {  # committed steps in a trajectory, or -1
    python - "$1" <<'PY'
import sys
try:
    import h5py, numpy as np
    with h5py.File(sys.argv[1], "r", locking=False) as h:
        cf = int(np.asarray(h["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, h["state/current_step"].shape[0])
        print(int(np.asarray(h["state/current_step"][cf-1]).reshape(-1)[-1]) if cf > 0 else 0)
except Exception:
    print(-1)
PY
}

retry=(); done_n=0
for t in $(seq 0 11); do
    f="${DATA}/$(tag_for "$t").h5"
    if [[ -e "$f" ]]; then st=$(progress "$f"); else st=0; fi
    if [[ "$st" -ge "$TARGET" ]]; then
        done_n=$((done_n+1)); printf "  task %2d  %8d steps  COMPLETE\n" "$t" "$st"; continue
    fi
    if [[ " ${live} " == *" ${t} "* ]]; then
        printf "  task %2d  %8d steps  live\n" "$t" "$st"; continue
    fi
    printf "  task %2d  %8d steps  RETRY\n" "$t" "$st"
    retry+=("$t")
done

echo "${done_n}/12 complete"
if [[ "${#retry[@]}" -gt 0 ]]; then
    ids=$(IFS=,; echo "${retry[*]}")
    echo "resubmitting task IDs ${ids}"
    PROD_ALLOW_RESUME=1 sbatch --export=ALL,PROD_ALLOW_RESUME=1 \
        --array="${ids}" "${CPP_DIR}/submit_prod_endstop.sbatch"
else
    echo "nothing to resubmit"
fi

if [[ "${done_n}" -eq 12 ]]; then
    touch "${ROOT}/.campaign_complete"
    echo "CAMPAIGN COMPLETE"
fi
