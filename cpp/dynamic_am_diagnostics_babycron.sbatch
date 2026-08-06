#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=am_diag_cron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_am_diagnostics/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_am_diagnostics/%x_%j.err
#SBATCH --no-requeue

set -euo pipefail

MODE=${DIAG_MODE:?DIAG_MODE must be single, crosslink, or memory}
CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SCRATCH_DIR=/scratch/users/xsli/sarcomere
TARGET=${DIAG_NSTEPS:-200000}
INTERVAL_MIN=${DIAG_WATCHDOG_INTERVAL_MIN:-15}
WATCHDOG=${CPP_DIR}/dynamic_am_diagnostics_babysitter.sh

case "${MODE}" in
    single) DATA_DIR=${AM_SINGLE_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_single_motor} ;;
    crosslink) DATA_DIR=${AM_XLINK_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_minimal_crosslink} ;;
    memory) DATA_DIR=${AM_MEMORY_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_status2_memory} ;;
    *) exit 1 ;;
esac

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
bash "${WATCHDOG}"

COMPLETE=$(HDF5_USE_FILE_LOCKING=FALSE python - "${MODE}" "${DATA_DIR}" "${TARGET}" <<'PY'
import sys
from pathlib import Path
import h5py

mode, directory, target = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
complete = 0
total_tasks = 160 if mode == "memory" else 80
for task_id in range(total_tasks):
    arm, rep = divmod(task_id, 20)
    if mode == "single":
        stem = f"am_single_F{('0','0p5','1','2')[arm]}_r{rep:02d}"
    elif mode == "crosslink":
        stem = f"am_xlink_{('free_noaa','free_persistentaa','fixed_noaa','fixed_persistentaa')[arm]}_r{rep:02d}"
    else:
        geometry = ("free", "fixed")[arm // 4]
        memory = ("0", "1ms", "5ms", "20ms")[arm % 4]
        stem = f"am_xlink_{geometry}_memory{memory}_r{rep:02d}"
    path = directory / f"{stem}.h5"
    try:
        with h5py.File(path, "r", locking=False) as handle:
            resume = "/resume/current_step" if "/resume/current_step" in handle else "/resume_prev/current_step"
            complete += int(int(handle[resume][...].reshape(-1)[-1]) >= target)
    except Exception:
        pass
print(complete)
PY
)

TOTAL=80
if [ "${MODE}" = memory ]; then TOTAL=160; fi
if [ "${COMPLETE}" -ge "${TOTAL}" ]; then
    echo "all ${TOTAL} ${MODE} diagnostic trajectories complete"
    if [ "${MODE}" = single ]; then
        python /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/aggregate_single_motor_diagnostics.py \
            "${DATA_DIR}" --dt "${AM_SINGLE_DT:-1e-5}" --k-strong-off "${AM_SINGLE_K_STRONG_OFF:-80}"
        python /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/check_am_kmc_statistics.py \
            "${DATA_DIR}" --dt "${AM_SINGLE_DT:-1e-5}" \
            --k-on "${AM_SINGLE_K_ON:-1000}" \
            --k-weak-off "${AM_SINGLE_K_WEAK_OFF:-200}" \
            --k-powerstroke "${AM_SINGLE_K_POWERSTROKE:-200}" \
            --k-strong-off "${AM_SINGLE_K_STRONG_OFF:-80}" \
            --output-prefix "${DATA_DIR}/explicit_crossbridge_kmc_check"
    elif [ "${MODE}" = crosslink ]; then
        python /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/aggregate_minimal_crosslink.py \
            "${DATA_DIR}" --dt "${AM_XLINK_DT:-1e-5}"
    else
        python /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/aggregate_minimal_crosslink.py \
            "${DATA_DIR}" --dt "${AM_MEMORY_DT:-1e-5}" \
            --output-prefix "${DATA_DIR}/status2_memory"
    fi
else
    echo "${COMPLETE}/${TOTAL} ${MODE} trajectories complete; scheduling next scan"
    sbatch --partition=owners --account=rotskoff \
        --begin="now+${INTERVAL_MIN}minutes" \
        --export="ALL,DIAG_MODE=${MODE},DIAG_NSTEPS=${TARGET}" \
        "${CPP_DIR}/dynamic_am_diagnostics_babycron.sbatch"
fi
