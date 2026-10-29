#!/bin/bash
#SBATCH --job-name=ghviz
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=2:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/build-%j.out
set -eu
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 -m py_compile build-gaphunt-rig.py
i=0
while IFS= read -r G; do
  i=$((i+1))
  python3 build-gaphunt-rig.py --groups "$G" --title "$G" --out grp${i}.html
done <<"EOF"
1 actin + 1 myosin (pose = standalone)
Partner ladder: N frozen parallel myosins
Scaffold disorder: 6 frozen myosins, translation off
Idealisation ladder at sigma = 0
Dwell sweep, FREE myosins
ACTN2 attempt v1 (geometry failed)
EOF
ls -la *.html
echo BUILD_DONE
