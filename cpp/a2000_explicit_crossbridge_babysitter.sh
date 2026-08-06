#!/bin/bash

set -euo pipefail
CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${EXPLICIT_CB_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_explicit_crossbridge_v3}"
TARGET="${EXPLICIT_CB_NSTEPS:-300000}"
STALE_MIN="${EXPLICIT_CB_STALE_MINUTES:-30}"
mkdir -p "${DATA_DIR}"
exec 8>/scratch/users/xsli/sarcomere/.a2000_explicit_crossbridge_watchdog.lock
flock -n 8 || exit 0

live=$(squeue --me -h -o '%j' 2>/dev/null || true)
if grep -Fxq a2k_explicit_cb <<<"${live}"; then
    echo "production array still has pending elements; deferring"
    exit 0
fi

mapfile -t retry < <(HDF5_USE_FILE_LOCKING=FALSE python - \
    "${DATA_DIR}" "${TARGET}" "${STALE_MIN}" "${live}" <<'PY'
import sys, time
from pathlib import Path
import h5py

directory, target, stale, live = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), set(sys.argv[4].splitlines())
now = time.time()
for rep in range(3):
    stem = f"A2000_M1000_explicit_cb_lt1e-2_r{rep:02d}_L12x2x2"
    path = directory / f"{stem}.h5"
    if stem in live:
        continue
    if not path.exists() or path.stat().st_size == 0:
        print(rep); continue
    if now - path.stat().st_mtime < stale * 60:
        continue
    try:
        with h5py.File(path, "r", locking=False) as h:
            group = "/resume" if "/resume/current_step" in h else "/resume_prev"
            step = int(h[f"{group}/current_step"][...].reshape(-1)[-1])
    except Exception:
        step = -1
    if step < target:
        print(rep)
PY
)

if [ "${#retry[@]}" -gt 0 ]; then
    ids=$(IFS=,; echo "${retry[*]}")
    echo "retrying production task IDs ${ids}"
    sbatch --array="${ids}" "${CPP_DIR}/a2000_explicit_crossbridge.sbatch"
fi

complete=$(HDF5_USE_FILE_LOCKING=FALSE python - "${DATA_DIR}" "${TARGET}" <<'PY'
import sys
from pathlib import Path
import h5py
d, target = Path(sys.argv[1]), int(sys.argv[2])
n = 0
for rep in range(3):
    p=d/f"A2000_M1000_explicit_cb_lt1e-2_r{rep:02d}_L12x2x2.h5"
    try:
        with h5py.File(p,"r",locking=False) as h:
            g="/resume" if "/resume/current_step" in h else "/resume_prev"
            n += int(int(h[f"{g}/current_step"][...].reshape(-1)[-1]) >= target)
    except Exception:
        pass
print(n)
PY
)
if [ "${complete}" -eq 3 ]; then
    if ! squeue --me -h -n a2k_explicit_analysis | grep -q .; then
        sbatch "${CPP_DIR}/a2000_explicit_crossbridge_analysis.sbatch"
    fi
    echo "all production trajectories complete; analysis submitted"
else
    echo "${complete}/3 complete"
fi
