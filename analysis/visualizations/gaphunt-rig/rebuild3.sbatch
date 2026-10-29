#!/bin/bash
#SBATCH --job-name=ghviz3
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=1:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/rebuild3-%j.out
set -eu
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 build-gaphunt-rig.py --groups 'One actin shared by 3 parallel myosins (vs 1)' --title 'One actin, three parallel myosins' --out tri.html
python3 build-gaphunt-rig.py --groups 'Tilt scan, productive polarity, Z-disc stop gated on ACTN2' --title 'Tilt scan, Z-disc stop gated on ACTN2' --out tiltq.html
python3 build-gaphunt-rig.py --groups 'Tilt scan: 1 actin + 1 FREE myosin, strong heads' --title 'Tilt scan: one actin, one free myosin, strong heads' --out tilt.html
for f in tri tiltq tilt; do python3 validate-minimal-wrap.py $f.html | grep -E '^PASS|FAIL'; grep -c 'starts at 36' $f.html || true; done
echo REBUILD3_DONE
