#!/bin/bash
#SBATCH --job-name=tiltviz
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=1:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/tiltbuild-%j.out
set -u
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
grep -c "exit=0" /scratch/users/xsli/sarcomere/gaphunt_20261023/logs/tiltscan-49539937_*.out
python3 tilt_an.py
python3 build-gaphunt-rig.py --groups "Tilt scan: 1 actin + 1 FREE myosin, strong heads" \
   --title "Tilt scan: one actin, one free myosin, strong heads" --out tilt.html
python3 validate-minimal-wrap.py tilt.html | grep -E "^PASS|FAIL"
echo TILTBUILD_DONE
