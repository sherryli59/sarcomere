#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --time=04:00:00
#SBATCH --job-name=local_aa_analysis
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/analysis_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/analysis_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
export HDF5_USE_FILE_LOCKING=FALSE

REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
ROOT="${LOCAL_AA_ROOT:-/scratch/users/xsli/sarcomere/data_local_aa_assembly_v1}"
N_ACTINS="${LOCAL_AA_N_ACTINS:-2000}"
N_MYOSINS="${LOCAL_AA_N_MYOSINS:-1000}"
OUT="${LOCAL_AA_ANALYSIS:-${REPO}/analysis/local_aa_assembly_A${N_ACTINS}_M${N_MYOSINS}_v1}"
PRE="${LOCAL_AA_PRE_STEPS:-50000}"
DT="${LOCAL_AA_DT:-1e-6}"
python "${REPO}/analysis/analyze_local_aa_assembly_pilot.py" \
    "${ROOT}" "${OUT}" --dt "${DT}" --pre-steps "${PRE}" --phi-stride 10

# A reduced campaign is a numerical/kinetic gate for the 2000/1000 array.
# Promotion is automatic only after every reduced replicate passes.
if [ "${N_ACTINS}" -lt 2000 ]; then
    python "${REPO}/analysis/gate_local_aa_assembly_pilot.py" \
        "${OUT}/local_aa_assembly_per_run.csv"
    FULL_ROOT=/scratch/users/xsli/sarcomere/data_local_aa_assembly_v1
    FULL_PARAM=/scratch/users/xsli/sarcomere/params_local_aa_assembly_v1
    marker="${FULL_ROOT}/.submitted_from_reduced_gate"
    mkdir -p "${FULL_ROOT}" "${FULL_PARAM}"
    if [ ! -e "${marker}" ]; then
        full_job=$(sbatch --parsable --export=ALL,LOCAL_AA_ROOT="${FULL_ROOT}",LOCAL_AA_PARAM="${FULL_PARAM}",LOCAL_AA_N_ACTINS=2000,LOCAL_AA_N_MYOSINS=1000,LOCAL_AA_LX=12,LOCAL_AA_LY=2,LOCAL_AA_LZ=2,LOCAL_AA_PRE_STEPS=50000,LOCAL_AA_RAMP_STEPS=25000,LOCAL_AA_NSTEPS=200000,LOCAL_AA_SAVE_EVERY=1000 "${REPO}/cpp/local_aa_assembly_pilot.sbatch")
        baby_job=$(sbatch --parsable --export=ALL,LOCAL_AA_ROOT="${FULL_ROOT}",LOCAL_AA_PARAM="${FULL_PARAM}",LOCAL_AA_N_ACTINS=2000,LOCAL_AA_N_MYOSINS=1000,LOCAL_AA_LX=12,LOCAL_AA_LY=2,LOCAL_AA_LZ=2,LOCAL_AA_PRE_STEPS=50000,LOCAL_AA_RAMP_STEPS=25000,LOCAL_AA_NSTEPS=200000,LOCAL_AA_SAVE_EVERY=1000 "${REPO}/cpp/local_aa_assembly_babysitter.sbatch")
        printf '%s %s\n' "${full_job}" "${baby_job}" >"${marker}"
        echo "promoted full-scale jobs: worker=${full_job} babysitter=${baby_job}"
    fi
fi
