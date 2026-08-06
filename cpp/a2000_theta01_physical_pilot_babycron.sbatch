#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=a2k_theta01_babycron
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_theta01_physical_pilot/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_theta01_physical_pilot/babycron_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
WORKER="${CPP_DIR}/a2000_theta01_physical_pilot.sbatch"
OUT="${THETA_PILOT_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_theta01_physical_pilot_v2}"
LOCK=/scratch/users/xsli/sarcomere/a2000_theta01_physical_pilot_babycron.lock
NSTEPS="${THETA_PILOT_NSTEPS:-10000}"
SEED="${THETA_PILOT_SEED:-203001}"
INTERVAL="${THETA_PILOT_BABY_INTERVAL:-10}"
mkdir -p "${OUT}" "$(dirname "${LOCK}")"

exec 8>"${LOCK}"
flock -n 8 || exit 0

active_names=$(squeue -u "${USER}" -h -o '%100j' || true)
missing=()
for idx in 0 1; do
    if [ "${idx}" -eq 0 ]; then condition=ungated; else condition=theta0p1; fi
    tag="A2000_M1000_physical_v2_${condition}_seed${SEED}_L12x2x2"
    output="${OUT}/${tag}.h5"

    complete=0
    if [ -s "${output}" ] && python - "${output}" "${NSTEPS}" <<'PY'
import h5py
import sys

path, target = sys.argv[1], int(sys.argv[2])
with h5py.File(path, "r", locking=False) as handle:
    steps = []
    for key in ("/resume/current_step", "/resume_prev/current_step", "/state/current_step"):
        if key in handle:
            steps.append(int(handle[key][...].reshape(-1)[-1]))
raise SystemExit(0 if steps and max(steps) >= target else 1)
PY
    then
        complete=1
    fi

    if [ "${complete}" -eq 1 ]; then
        echo "[${tag}] complete"
    elif grep -Fxq "${tag}" <<<"${active_names}"; then
        echo "[${tag}] worker still queued/running"
    else
        echo "[${tag}] missing or incomplete; scheduling task ${idx}"
        missing+=("${idx}")
    fi
done

if [ "${#missing[@]}" -gt 0 ]; then
    array_spec=$(IFS=,; echo "${missing[*]}")
    sbatch --array="${array_spec}" "${WORKER}"
fi

# Chain one future audit while work remains. Singleton locking and worker-name
# checks prevent duplicate writers even if two audits briefly overlap.
if [ "${#missing[@]}" -gt 0 ] || grep -Eq 'A2000_M1000_physical_v2_(ungated|theta0p1)_seed' <<<"${active_names}"; then
    sbatch --begin="now+${INTERVAL}minutes" "$0"
else
    echo "Both paired pilot trajectories are complete; babysitter stopping."
fi
