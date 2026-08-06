#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --job-name=a2k_explicit_analysis
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%j_analysis.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%j_analysis.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

ROOT=/home/groups/rotskoff/sherry/sarcomere/sarcomere
INPUT="${EXPLICIT_CB_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_explicit_crossbridge_v3}"
OUT="${EXPLICIT_CB_ANALYSIS_DIR:-${ROOT}/analysis/a2000_explicit_crossbridge_v3}"
mkdir -p "${OUT}"
saved_only=()
if [ "${EXPLICIT_CB_ANALYSIS_SAVED_ONLY:-0}" = 1 ]; then
    saved_only=(--saved-only)
fi
structural_only=0
if [ "${EXPLICIT_CB_ANALYSIS_STRUCTURAL_ONLY:-0}" = 1 ]; then
    structural_only=1
fi

for trajectory in "${INPUT}"/A2000_M1000_explicit_cb_*.h5; do
    [ -e "${trajectory}" ] || continue
    stem=$(basename "${trajectory}" .h5)
    python "${ROOT}/analysis/analyze_actn2_force_projection_run.py" \
        "${trajectory}" --output-prefix "${OUT}/${stem}" --dt 1e-5 --box 12 2 2 \
        --frame-stride "${EXPLICIT_CB_ANALYSIS_STRIDE:-1}"
    if [ "${structural_only}" != 1 ]; then
        python "${ROOT}/analysis/summarize_am_transitions.py" \
            "${trajectory}" --dt 1e-5 --output-prefix "${OUT}/${stem}" \
            "${saved_only[@]}"
        python "${ROOT}/analysis/summarize_aa_flicker.py" \
            "${trajectory}" --dt 1e-5 --tau1 1e-3 --output-prefix "${OUT}/${stem}" \
            "${saved_only[@]}"
    fi
done

aggregate_mode=()
if [ "${structural_only}" = 1 ]; then
    aggregate_mode=(--structural-only)
fi
python "${ROOT}/analysis/aggregate_explicit_crossbridge_runs.py" \
    "${INPUT}" "${OUT}" --dt 1e-5 \
    --k-weak-off 200 --k-powerstroke 200 --k-strong-off 80 --force-scale 1 \
    "${saved_only[@]}" "${aggregate_mode[@]}"
touch "${OUT}/.analysis_complete"
