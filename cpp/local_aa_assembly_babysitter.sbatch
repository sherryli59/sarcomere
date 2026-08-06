#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=local_aa_babycron
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/babycron_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
WORKER="${CPP_DIR}/local_aa_assembly_pilot.sbatch"
ROOT="${LOCAL_AA_ROOT:-/scratch/users/xsli/sarcomere/data_local_aa_assembly_v1}"
LOCK="${ROOT}/.babysitter.lock"
NSTEPS="${LOCAL_AA_NSTEPS:-200000}"
N_ACTINS="${LOCAL_AA_N_ACTINS:-2000}"
N_MYOSINS="${LOCAL_AA_N_MYOSINS:-1000}"
SEED_BASE="${LOCAL_AA_SEED_BASE:-204000}"
LX="${LOCAL_AA_LX:-12}"; LY="${LOCAL_AA_LY:-2}"; LZ="${LOCAL_AA_LZ:-2}"
INTERVAL="${LOCAL_AA_BABY_INTERVAL:-20}"
mkdir -p "${ROOT}"
exec 8>"${LOCK}"; flock -n 8 || exit 0

active=$(squeue -u "${USER}" -h -o '%100j' || true)
worker_active=0
if grep -Fxq "local_aa_pilot" <<<"${active}"; then worker_active=1; fi
conditions=(aa_off one_nocatch two_nocatch two_status2 two_status1)
missing=()
for task in $(seq 0 14); do
    ci=$((task / 3)); rep=$((task % 3)); seed=$((SEED_BASE + rep))
    tag="A${N_ACTINS}_M${N_MYOSINS}_${conditions[$ci]}_r${rep}_seed${seed}_L${LX}x${LY}x${LZ}"
    output="${ROOT}/${tag}.h5"
    if [ -s "${output}" ] && python - "${output}" "${NSTEPS}" <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], "r", locking=False) as h:
    steps=[]
    for key in ("/resume/current_step", "/resume_prev/current_step", "/state/current_step"):
        if key in h: steps.append(int(h[key][...].reshape(-1)[-1]))
raise SystemExit(0 if steps and max(steps) >= int(sys.argv[2]) else 1)
PY
    then
        echo "[${tag}] complete"
    elif grep -Fxq "${tag}" <<<"${active}"; then
        echo "[${tag}] active"
    else
        missing+=("${task}")
    fi
done

if [ "${#missing[@]}" -gt 0 ] && [ "${worker_active}" -eq 0 ]; then
    spec=$(IFS=,; echo "${missing[*]}")
    sbatch --export=ALL --array="${spec}%6" "${WORKER}"
fi
if [ "${#missing[@]}" -gt 0 ] || [ "${worker_active}" -eq 1 ] || grep -Eq "A${N_ACTINS}_M${N_MYOSINS}_(aa_off|one_nocatch|two_nocatch|two_status2|two_status1)" <<<"${active}"; then
    sbatch --export=ALL --begin="now+${INTERVAL}minutes" "$0"
else
    analysis_marker="${ROOT}/.analysis_submitted"
    if [ ! -e "${analysis_marker}" ]; then
        analysis_job=$(sbatch --parsable --export=ALL "${CPP_DIR}/local_aa_assembly_analysis.sbatch")
        printf '%s\n' "${analysis_job}" >"${analysis_marker}"
        echo "submitted terminal analysis ${analysis_job}"
    fi
fi
