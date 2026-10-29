#!/bin/bash
#SBATCH --job-name=snapviz
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/amalign_20261005/logs/snapviz_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/amalign_20261005/logs/snapviz_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
echo "host=$(hostname) start=$(date '+%F %T')"
V=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/latest-snapshots
OUT=/scratch/users/xsli/sarcomere/amalign_20261005/analysis
mkdir -p "$OUT"
S=20260805
A=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_capsule_normal.h5
B=/scratch/users/xsli/sarcomere/hph_ladder_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_spokes1.h5
C=/scratch/users/xsli/sarcomere/hph_ladder_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_hph27.h5
A9=/scratch/users/xsli/sarcomere/amalign_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_align9.h5
C27=/scratch/users/xsli/sarcomere/amalign_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_align27.h5
for f in "$A" "$B" "$C" "$A9" "$C27" "$V/build-latest-snapshots.py" "$V/latest-snapshots-template.html"; do
  [ -r "$f" ] || { echo "ABORT missing $f"; exit 1; }
done
echo "--- committed time of each run (these DIFFER; the viewer shows it per frame)"
python3 - "$A" "$B" "$C" "$A9" "$C27" <<'PY'
import sys, os
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py, numpy as np
for p in sys.argv[1:]:
    with h5py.File(p,"r",locking=False) as h:
        cf=int(np.asarray(h["meta/committed_frames"]).reshape(-1)[-1])
        cf=min(cf,h["state/current_step"].shape[0])
        st=int(np.asarray(h["state/current_step"][cf-1]).reshape(-1)[-1])
        print("   %-10s frames=%-5d step=%-7d t=%7.1f ms" %
              (os.path.basename(p).split("_5000ms_")[1].replace(".h5",""), cf, st, st*float(h["meta"].attrs["dt"])*1e3))
PY
python3 "$V/build-latest-snapshots.py" "$A" "$B" "$C" "$A9" "$C27" \
  --labels "A  hph9 spokes3  k=0  (baseline)" \
           "B  hph9 spokes1  k=0  (75nm crowns, same heads)" \
           "C  hph27 spokes3 k=0  (75nm crowns, 3x heads)" \
           "A9 hph9 spokes3  k=30 (alignment on)" \
           "C27 hph27 spokes3 k=30 (3x heads + alignment)" \
  --last 3 \
  --title "Latest committed frame per arm, seed 20260805 - PANELS ARE AT DIFFERENT TIMES" \
  --template "$V/latest-snapshots-template.html" \
  --out "$OUT/latest-snapshots.html"
echo "--- result"; ls -la "$OUT/latest-snapshots.html"
python3 -c "
s=open('$OUT/latest-snapshots.html',encoding='ascii').read()
print('ASCII:', s.isascii())
print('markers replaced:', '/*__RUNS__*/' not in s and '__TITLE__' not in s)
print('time warning present:', 'DIFFERENT SIMULATED TIMES' in s)
"
echo "end=$(date '+%F %T')"
