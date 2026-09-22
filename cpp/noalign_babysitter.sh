#!/bin/bash
# Single-pass scan for the noalign production array.  Follows the established
# babysitter conventions in this repo: .inflight markers checked for QUEUE
# LIVENESS (never an age timer), and a hard skip when the worker itself is
# still queued.
#
# Adapted from the per-file babysitters: this campaign is an ARRAY, so the unit
# of work is an array index, not a filename.  Index -> (arm, seed) mirrors the
# launcher: arm = idx/6, seed = idx%6.
set -uo pipefail

ROOT=${NOALIGN_ROOT:-/scratch/users/xsli/sarcomere/production_noalign_5s_20260921}
SBATCH_FILE=${NOALIGN_SBATCH:-${ROOT}/source_snapshot/submit_production_noalign_5s_20260921.sbatch}
PARTITIONS=${NOALIGN_PARTITIONS:-rotskoff,owners}
TARGET_FRAMES=${NOALIGN_TARGET_FRAMES:-1001}
JOBNAME=${NOALIGN_JOBNAME:-noalign}
CPUS=${NOALIGN_CPUS:-32}
MEM=${NOALIGN_MEM:-16G}
WALL=${NOALIGN_WALL:-1-23:00:00}
INFLIGHT="${ROOT}/.inflight"
mkdir -p "${INFLIGHT}" "${ROOT}/logs"
LOG="${ROOT}/logs/babysitter.log"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

SEEDS=(20260805 20260815 20260825 20260835 20260845 20260855)
ARMS=(titinon titinoff)

# Array indices currently queued or running under our own job name.  GUARD 1:
# owners is PreemptMode=REQUEUE, so a preempted worker comes BACK under the
# same job; resuming alongside it would double-write the HDF5.
LIVE_IDX=$(squeue --me -h -n "${JOBNAME}" -o "%K" 2>/dev/null \
           | tr ',' '\n' | sed 's/%.*//' \
           | awk -F- '/^[0-9]+$/{print} /^[0-9]+-[0-9]+$/{split($0,a,"-"); for(i=a[1];i<=a[2];i++) print i}')

NEED=()
for idx in $(seq 0 11); do
    arm=${ARMS[$(( idx / 6 ))]}
    seed=${SEEDS[$(( idx % 6 ))]}
    f=$(ls "${ROOT}"/data/*_"${arm}"_seed"${seed}"_*.h5 2>/dev/null | head -1)
    cf=0
    if [ -n "${f}" ]; then
        cf=$(python3 - "$f" <<'PY'
import sys, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
try:
    import h5py
    with h5py.File(sys.argv[1], "r", locking=False) as h:
        print(int(h["meta/committed_frames"][0]) if "meta/committed_frames" in h else 0)
except Exception:
    print(0)
PY
)
    fi
    cf=$(echo "$cf" | tr -d '[:space:]'); [[ "$cf" =~ ^[0-9]+$ ]] || cf=0
    if [ "$cf" -ge "$TARGET_FRAMES" ]; then continue; fi

    # GUARD 1: worker still in the queue (running, or requeued after preemption)
    if echo "${LIVE_IDX}" | grep -qx "${idx}"; then
        say "skip idx ${idx} (${arm}/${seed}, ${cf}/${TARGET_FRAMES}): worker still queued"
        continue
    fi
    # GUARD 2: a resume we launched is still in the queue.  Liveness, not age:
    # under cluster priority a resume can sit PD for hours and an age timer
    # would fire duplicates.
    M="${INFLIGHT}/${idx}"
    if [ -f "$M" ]; then
        jid=$(cat "$M" 2>/dev/null)
        if [ -n "$jid" ] && squeue --me -h -j "$jid" -o "%i" 2>/dev/null | grep -q .; then
            say "skip idx ${idx} (${arm}/${seed}, ${cf}/${TARGET_FRAMES}): resume ${jid} still queued"
            continue
        fi
        rm -f "$M"   # stale: recorded resume no longer exists
    fi
    say "NEEDS RESUME idx ${idx} (${arm}/${seed}): ${cf}/${TARGET_FRAMES} frames"
    NEED+=("$idx")
done

if [ "${#NEED[@]}" -eq 0 ]; then
    say "nothing to resume"
    exit 0
fi
LIST=$(IFS=,; echo "${NEED[*]}")
# --requeue on purpose: owners has GraceTime=0 so a trap cannot fire before the
# kill, but PreemptMode=REQUEUE means Slurm will put the job back itself, and
# the launcher resumes from committed frames (SLURM_RESTART_COUNT bypasses its
# overwrite guard).
NEW=$(sbatch --parsable --partition="${PARTITIONS}" --requeue \
        --time="${WALL}" --mem="${MEM}" --cpus-per-task="${CPUS}" \
        --export=ALL,ASSEMBLY_PILOT_ROOT="${ROOT}",ASSEMBLY_PILOT_RESUME_EXISTING=true \
        --array="${LIST}" \
        --output="${ROOT}/logs/%x_%A_%a.out" --error="${ROOT}/logs/%x_%A_%a.err" \
        "${SBATCH_FILE}" 2>>"$LOG")
if [ -n "${NEW}" ]; then
    say "submitted resume array ${NEW} for indices [${LIST}]"
    for idx in "${NEED[@]}"; do echo "${NEW}_${idx}" > "${INFLIGHT}/${idx}"; done
else
    say "sbatch FAILED for [${LIST}]"
fi
