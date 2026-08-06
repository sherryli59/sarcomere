#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --job-name=sarc_3d_render
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/render_logs/%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/render_logs/%j.err

set -euo pipefail
if [ "$#" -lt 5 ]; then
    echo "usage: sbatch $0 TRAJECTORY OUTPUT_DIR LX LY LZ [STRIDE] [FPS]" >&2
    exit 2
fi

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
export HDF5_USE_FILE_LOCKING=FALSE
export PYVISTA_OFF_SCREEN=true

REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
TRAJECTORY=$1; OUTPUT_DIR=$2; LX=$3; LY=$4; LZ=$5
STRIDE=${6:-2}; FPS=${7:-15}
mkdir -p "${OUTPUT_DIR}" /scratch/users/xsli/sarcomere/render_logs
stem=$(basename "${TRAJECTORY}" .h5)

python "${REPO}/analysis/animate_sarcomere_3d.py" "${TRAJECTORY}" \
    "${OUTPUT_DIR}/${stem}_status_isometric.mp4" --box "${LX}" "${LY}" "${LZ}" \
    --stride "${STRIDE}" --fps "${FPS}" --camera isometric \
    --style clean --selection loaded --actin-color status

python "${REPO}/analysis/animate_sarcomere_3d.py" "${TRAJECTORY}" \
    "${OUTPUT_DIR}/${stem}_status_transverse.mp4" --box "${LX}" "${LY}" "${LZ}" \
    --stride "${STRIDE}" --fps "${FPS}" --camera transverse \
    --style clean --selection loaded --parallel-projection --actin-color status

python "${REPO}/analysis/animate_sarcomere_3d.py" "${TRAJECTORY}" \
    "${OUTPUT_DIR}/${stem}_polarity_isometric.mp4" --box "${LX}" "${LY}" "${LZ}" \
    --stride "${STRIDE}" --fps "${FPS}" --camera isometric \
    --style clean --selection loaded --actin-color polarity --show-polarity
