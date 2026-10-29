#!/bin/bash
#SBATCH --job-name=spokesviz
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/hph_ladder_20261005/logs/spokesviz_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/hph_ladder_20261005/logs/spokesviz_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
echo "host=$(hostname) start=$(date '+%F %T')"
V=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/ksoff27-prod
OUT=/scratch/users/xsli/sarcomere/hph_ladder_20261005/analysis
mkdir -p "$OUT"
S=20260835
A=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_capsule_normal.h5
B=/scratch/users/xsli/sarcomere/hph_ladder_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed${S}_5000ms_spokes1.h5
for f in "$A" "$B" "$V/build-ksoff27-prod.py" "$V/ksoff27-prod-template.html"; do
  [ -r "$f" ] || { echo "ABORT missing $f"; exit 1; }
done
# both must be complete, or the two panels are at different times
python3 - "$A" "$B" <<'PY'
import sys, os
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py, numpy as np
for p in sys.argv[1:]:
    with h5py.File(p,"r",locking=False) as h:
        cf=int(np.asarray(h["meta/committed_frames"]).reshape(-1)[-1]); cf=min(cf,h["state/current_step"].shape[0])
        st=int(np.asarray(h["state/current_step"][cf-1]).reshape(-1)[-1])
        a=h["meta"].attrs
        print("   %-16s frames=%-5d t=%7.1f ms  hph=%-3d spokes=%-2s align_k=%.0f"
              % (os.path.basename(p).split("_5000ms_")[1].replace(".h5",""), cf, st*float(a["dt"])*1e3,
                 int(a["am_heads_per_half"]), int(float(a.get("am_spokes_per_crown",3))), float(a["am_alignment_k"])))
        assert cf >= 1001, "run not complete: %d frames" % cf
print("   both complete at 5000 ms -- matched time, matched seed")
PY
python3 "$V/build-ksoff27-prod.py" "$A" "$B" \
  --labels "A  spokes=3  (3 crowns/half, 225 nm)" \
           "B  spokes=1  (9 crowns/half, 75 nm, SAME 18 heads)" \
  --stride 50 \
  --title "am_spokes_per_crown 3 vs 1 at fixed head count, seed 20260835, 5000 ms" \
  --template "$V/ksoff27-prod-template.html" \
  --out "$OUT/spokes-traj.html"
echo "--- result"; ls -la "$OUT/spokes-traj.html"
python3 -c "
s=open('$OUT/spokes-traj.html',encoding='ascii').read()
print('ASCII:', s.isascii())
print('markers replaced:', '/*__RUNS__*/' not in s and '__TITLE__' not in s)
print('production rig named:', 'myosin_x_noise' in s)
"
echo "end=$(date '+%F %T')"
