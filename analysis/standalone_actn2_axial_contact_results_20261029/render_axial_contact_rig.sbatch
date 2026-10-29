#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --job-name=actn2_rig_render
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/render_logs/%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/render_logs/%j.err

# Render the weighted vs axial-stop rig with analysis/animate_sarcomere_3d.py (camera zoom wrapper only).
# usage: sbatch render_axial_contact_rig.sbatch RESULTS_DIR
set -euo pipefail
RESULTS=${1:?results directory holding trajectories.npz}
REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
ENV=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere
export PATH="$ENV/bin:$PATH" HDF5_USE_FILE_LOCKING=FALSE PYVISTA_OFF_SCREEN=true
export MPLCONFIGDIR=/tmp/sarcomere_matplotlib_$SLURM_JOB_ID
OUT="$RESULTS/movies"; H5="$RESULTS/h5"
mkdir -p "$OUT" "$H5"
python "$RESULTS/export_rig_h5.py" "$RESULTS/trajectories.npz" "$H5"

BOX=(3.3 0.5 0.5)
COMMON=(--box "${BOX[@]}" --style clean --selection all --actin-color status --show-links
        --show-forces both --actin-width 5 --myosin-width 12 --link-width 5
        --fps 20 --no-dashboard --width 1200 --height 700 --font-size 11)
for case in escape_run1 c1_removed; do
  for law in weighted axial; do
    for camera in transverse axial isometric; do
      pp=""; [ "$camera" != isometric ] && pp=--parallel-projection
      python "$RESULTS/render_rig_with_pipeline.py" "$H5/${case}_${law}.h5" \
        "$OUT/${case}_${law}_${camera}.mp4" --camera "$camera" "${COMMON[@]}" $pp
      python "$RESULTS/render_rig_with_pipeline.py" "$H5/${case}_${law}.h5" \
        "$OUT/${case}_${law}_${camera}_final.png" --camera "$camera" "${COMMON[@]}" $pp
    done
  done
  # Side-by-side comparison per camera: weighted (left) | axial stop (right).
  label() { echo "drawtext=text='$1':x=20:y=h-50:fontsize=32:fontcolor=black:box=1:boxcolor=white@0.85:boxborderw=8"; }
  for camera in transverse axial isometric; do
    name=$camera; [ "$camera" = axial ] && name=endon
    ffmpeg -y -loglevel error -i "$OUT/${case}_weighted_${camera}.mp4" -i "$OUT/${case}_axial_${camera}.mp4" \
      -filter_complex "[0:v]$(label 'weighted energy, alpha=3')[a];[1:v]$(label 'axial stop')[b];[a][b]hstack" \
      -c:v libx264 -crf 18 -pix_fmt yuv420p -movflags +faststart "$OUT/${case}_${name}_side_by_side.mp4"
    ffmpeg -y -loglevel error -sseof -0.1 -i "$OUT/${case}_${name}_side_by_side.mp4" -frames:v 1 "$OUT/.${case}_${name}.png"
  done
  ffmpeg -y -loglevel error -i "$OUT/.${case}_transverse.png" -i "$OUT/.${case}_endon.png" -i "$OUT/.${case}_isometric.png" \
    -filter_complex vstack=inputs=3 -frames:v 1 "$OUT/${case}_side_by_side_final.png"
  rm -f "$OUT"/.${case}_*.png
done
ls -la "$OUT"
