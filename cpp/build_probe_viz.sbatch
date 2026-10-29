#!/bin/bash
#SBATCH --job-name=probeviz
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/viz_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/viz_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
echo "host=$(hostname) start=$(date '+%F %T')"
V=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/probe-hph
OUT=/scratch/users/xsli/sarcomere/probe_hph_20261005/analysis
mkdir -p "$OUT"
H9=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910/data/span3_tilt60_seed1.h5
H27=/scratch/users/xsli/sarcomere/probe_hph_20261005/data/hph27_tilt60_seed1.h5
for f in "$H9" "$H27" "$V/build-probe-hph.py" "$V/probe-hph-template.html"; do
  [ -r "$f" ] || { echo "ABORT missing $f"; exit 1; }
done
# --azimuth auto is REQUIRED: this binary does not persist am_head_azimuth and
# the builder's absent-attr default would assume ON, mis-modelling the radial
# channel. sarcomere.h:642 makes the real default false.
python3 "$V/build-probe-hph.py" "$H9" "$H27" \
  --labels "hph=9  (225 nm crowns) - does NOT repair" "hph=27 (75 nm crowns) - REPAIRS" \
  --azimuth auto --stride 20 \
  --title "Confined probe, actin started at 60 deg: can the AM heads repair it? (seed 1)" \
  --template "$V/probe-hph-template.html" \
  --out "$OUT/probe-hph-tilt60-traj.html"
echo "--- result"; ls -la "$OUT/probe-hph-tilt60-traj.html"
python3 -c "
s=open('$OUT/probe-hph-tilt60-traj.html',encoding='ascii').read()
print('ASCII:', s.isascii())
print('markers replaced:', '/*__RUNS__*/' not in s and '__TITLE__' not in s)
print('correct rig:', 'single_am_tilted' in s and 'myosin_x_noise' not in s)
"
echo "end=$(date '+%F %T')"
