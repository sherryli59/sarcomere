#!/bin/bash
#SBATCH --job-name=tpviz
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=1:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/tpbuild-%j.out
set -u
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 tilt_an.py tp tiltp_summary.json
python3 build-gaphunt-rig.py --groups 'Tilt scan, PRODUCTIVE polarity, pointed end outside the mid zone' --title 'Tilt scan, productive polarity' --out tiltp.html
python3 validate-minimal-wrap.py tiltp.html | grep -E '^PASS|FAIL'
echo TPBUILD_DONE
