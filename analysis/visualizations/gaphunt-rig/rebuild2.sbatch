#!/bin/bash
#SBATCH --job-name=ghviz2
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=2:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/rebuild2-%j.out
set -eu
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 -m py_compile build-gaphunt-rig.py
python3 tilt_an.py > tilt_an.out
python3 build-gaphunt-rig.py --groups "Tilt scan: 1 actin + 1 FREE myosin, strong heads" \
   --title "Tilt scan: one actin, one free myosin, strong heads" --out tilt.html
i=0
while IFS= read -r G; do
  i=$((i+1)); python3 build-gaphunt-rig.py --groups "$G" --title "$G" --out grp${i}.html
done <<"EOF"
1 actin + 1 myosin (pose = standalone)
Partner ladder: N frozen parallel myosins
Scaffold disorder: 6 frozen myosins, translation off
Idealisation ladder at sigma = 0
Dwell sweep, FREE myosins
ACTN2 attempt v1 (geometry failed)
EOF
for f in tilt grp1 grp2 grp3 grp4 grp5 grp6; do python3 validate-minimal-wrap.py $f.html | grep -E "^PASS|FAIL"; done
echo REBUILD2_DONE
