#!/bin/bash
#SBATCH --partition=hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=watch_A1280_M640
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_L7p72_20260723/logs/watch_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_L7p72_20260723/logs/watch_%j.err

set -u

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
export HDF5_USE_FILE_LOCKING=FALSE

REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_L7p72_20260723}"
WORKER="${REPO}/cpp/production_L7p72_A1280_M640.sbatch"
MONITOR="${REPO}/analysis/monitor_L7p72_production.py"
OUTPUT="${PROD_OUTPUT:-${ROOT}/data/A1280_M640_L7p72_seed2072301.h5}"
CSV="${ROOT}/monitor/health.csv"
TARGET="${PROD_NSTEPS:-1000000}"
LOCK="${ROOT}/watch.lock"
INTERVAL_MINUTES="${PROD_WATCH_INTERVAL_MINUTES:-10}"
mkdir -p "${ROOT}/logs" "${ROOT}/monitor"

exec 8>"${LOCK}"
flock -n 8 || exit 0

worker_line=$(squeue -u "${USER}" -h -n prod_A1280_M640_L7p72 -o '%i|%T|%M' | head -1)
worker_id=""
elapsed_seconds=0
if [ -n "${worker_line}" ]; then
    IFS='|' read -r worker_id worker_state worker_elapsed <<<"${worker_line}"
    elapsed_seconds=$(python - "${worker_elapsed}" <<'PY'
import sys
text=sys.argv[1]
days=0
if "-" in text:
    d,text=text.split("-",1); days=int(d)
parts=[int(x) for x in text.split(":")]
if len(parts)==3: h,m,s=parts
elif len(parts)==2: h=0; m,s=parts
else: h=0; m=0; s=parts[0]
print(days*86400+h*3600+m*60+s)
PY
)
    echo "worker=${worker_id} state=${worker_state} elapsed=${worker_elapsed}"
else
    echo "no active worker found"
fi

monitor_status=0
if [ -s "${OUTPUT}" ]; then
    python "${MONITOR}" "${OUTPUT}" --csv "${CSV}" \
        --job-id "${worker_id}" --elapsed-seconds "${elapsed_seconds}" \
        --target-steps "${TARGET}" || monitor_status=$?
fi

if [ "${monitor_status}" -eq 42 ] || [ "${monitor_status}" -eq 43 ]; then
    reason="nonfinite_state"
    [ "${monitor_status}" -eq 43 ] && reason="accelerating_force_caps"
    echo "SAFETY STOP: ${reason}"
    if [ -n "${worker_id}" ]; then
        scancel "${worker_id}"
    fi
    touch "${ROOT}/monitor/STOPPED_${reason}"
    exit 0
fi

complete=0
if [ -s "${OUTPUT}" ]; then
    complete=$(python - "${OUTPUT}" "${TARGET}" <<'PY'
import h5py, sys
try:
    with h5py.File(sys.argv[1], "r", locking=False) as h:
        values=[]
        for path in ("/resume/current_step", "/resume_prev/current_step", "/state/current_step"):
            if path in h:
                values.append(int(h[path][...].reshape(-1)[-1]))
        print(1 if values and max(values) >= int(sys.argv[2]) else 0)
except (OSError, KeyError, ValueError):
    print(0)
PY
)
fi

if [ "${complete}" -eq 1 ]; then
    echo "Production trajectory reached ${TARGET} steps; watcher stopping."
    exit 0
fi

if [ -z "${worker_id}" ]; then
    echo "Incomplete trajectory has no active worker; resubmitting checkpoint-aware worker."
    sbatch "${WORKER}"
fi

sbatch --begin="now+${INTERVAL_MINUTES}minutes" "$0"
