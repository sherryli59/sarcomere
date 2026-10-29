#!/bin/bash
#SBATCH --job-name=hphsnap
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=32G --time=2:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/hphnorm-snapshots/build-%j.out
set -eu
export OMP_NUM_THREADS=2 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/hphnorm-snapshots
D=/scratch/users/xsli/sarcomere/gaphunt_20261023/data
python3 build-latest-snapshots.py $D/hph27m24_seed20260805.h5 $D/hph54m24_seed20260815.h5 $D/hph81m24_seed20260805.h5 \
  --last 1 --labels "27 heads/half (seed 805), 2 s" "54 heads/half (seed 815), 2 s" "81 heads/half (seed 805), 2 s" \
  --title "Head-count production runs at 2 s, sliding force fixed" \
  --template latest-snapshots-template.html --out hphnorm-snapshots.html
echo SNAP_DONE
