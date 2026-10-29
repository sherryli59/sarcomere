#!/bin/bash
#SBATCH --job-name=tqviz
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=1:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/tqbuild-%j.out
set -u
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
grep -h requires_link /scratch/users/xsli/sarcomere/gaphunt_20261023/logs/tq30_seed20260805.argv
python3 - <<'PY'
import h5py,glob,numpy as np
mx=0; n=0
for p in glob.glob('/scratch/users/xsli/sarcomere/gaphunt_20261023/data/tq[0-9]*_seed*.h5'):
    f=h5py.File(p,'r',locking=False); k=int(f['meta/committed_frames'][0]); n+=1
    mx=max(mx,float(f['actin/aa_barbed_tip_endstop_force'][:k,0,0].max()))
print('runs %d, max Z-disc stop force applied %.3f pN (must be 0)'%(n,mx))
PY
python3 tilt_an.py tq tiltq_summary.json
python3 build-gaphunt-rig.py --groups 'Tilt scan, productive polarity, Z-disc stop gated on ACTN2' --title 'Tilt scan, Z-disc stop gated on ACTN2' --out tiltq.html
python3 validate-minimal-wrap.py tiltq.html | grep -E '^PASS|FAIL'
echo TQBUILD_DONE
