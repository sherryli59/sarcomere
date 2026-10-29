#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --job-name=actn2_rig_oblique
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/render_logs/%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/render_logs/%j.err

# Old vs new contact rule from a side-oblique camera, filaments at full length.
# usage: sbatch render_rig_oblique.sbatch RESULTS_DIR   (needs h5/ from render_axial_contact_rig.sbatch)
set -euo pipefail
RESULTS=${1:?results directory}
ENV=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere
export PATH="$ENV/bin:$PATH" HDF5_USE_FILE_LOCKING=FALSE PYVISTA_OFF_SCREEN=true
export MPLCONFIGDIR=/tmp/sarcomere_matplotlib_$SLURM_JOB_ID
OUT="$RESULTS/movies"; H5="$RESULTS/h5"
mkdir -p "$OUT"
COMMON=(--camera oblique --box 3.3 0.5 0.5 --style clean --selection all --actin-color status
        --show-links --show-forces myosin --actin-width 6 --myosin-width 16 --link-width 6
        --fps 20 --no-dashboard --width 1400 --height 560 --font-size 11)
label() { echo "drawtext=text='$1':x=(w-text_w)/2:y=h-58:fontsize=34:fontcolor=black:box=1:boxcolor=white@0.9:boxborderw=10"; }
for case in c1_removed escape_run1; do
  for law in weighted axial; do
    python "$RESULTS/render_rig_with_pipeline.py" "$H5/${case}_${law}.h5" "$OUT/${case}_${law}_oblique.mp4" "${COMMON[@]}"
  done
  ffmpeg -y -loglevel error -i "$OUT/${case}_weighted_oblique.mp4" -i "$OUT/${case}_axial_oblique.mp4" \
    -filter_complex "[0:v]$(label 'OLD rule (weighted energy)')[a];[1:v]$(label 'NEW rule (axial stop)')[b];[a][b]vstack" \
    -c:v libx264 -crf 18 -pix_fmt yuv420p -movflags +faststart "$OUT/${case}_oblique_old_vs_new.mp4"
  ffmpeg -y -loglevel error -i "$OUT/${case}_oblique_old_vs_new.mp4" \
    -vf "select=eq(n\,0)+eq(n\,20)+eq(n\,60)+eq(n\,200),scale=900:-1,tile=4x1" -frames:v 1 -fps_mode vfr \
    "$OUT/${case}_oblique_strip.png"
done
ls -la "$OUT" | grep oblique
