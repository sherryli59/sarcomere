#!/bin/bash
#SBATCH --job-name=probesmoke
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:25:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/smoke_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/smoke_%j.err
set -euo pipefail
export OMP_NUM_THREADS=1
src=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910
root=/scratch/users/xsli/sarcomere/probe_hph_20261005
mkdir -p "${root}"/{data,params,logs,smoke}
for hph in 9 27; do
  "${src}/run_conf2.sh" "${root}/smoke/smoke_hph${hph}.h5" \
    --Lx=4 --Ly=4 --Lz=4 --am_heads_per_half=${hph} \
    --am_head_transverse_k=100 --am_alignment_k=0 \
    --probe_tilt_deg=60 --probe_confine_k=50 \
    --nsteps=2000 --seed=1 >"${root}/logs/smoke_hph${hph}.log" 2>&1
  echo "ran hph=${hph}"
done
python3 - <<'PY'
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py
R="/scratch/users/xsli/sarcomere/probe_hph_20261005/smoke"
for hph in (9,27):
    p=f"{R}/smoke_hph{hph}.h5"
    with h5py.File(p,"r",locking=False) as h:
        a=h["meta"].attrs
        hs=h["state/am_head_state_current"]
        need=["am_align_drive_minus","am_align_drive_plus",
              "am_align_head_count_minus","am_align_head_count_plus"]
        miss=[k for k in need if k not in h["myosin"]]
        print("hph=%-3d attr_hph=%-3s tilt=%-5s confine_k=%-5s head_axis=%-3d cap/half=%s  missing=%s"
              % (hph, a.get("am_heads_per_half"), a.get("probe_tilt_deg"),
                 a.get("probe_confine_k"), hs.shape[2],
                 a.get("am_max_heads_per_actin_half"), miss or "none"))
        assert int(a["am_heads_per_half"])==hph, "hph did not reach the binary"
        assert hs.shape[2]==2*hph, "head axis %d != 2*%d" % (hs.shape[2],hph)
        assert not miss, miss
        bound=((np.asarray(hs)==1)|(np.asarray(hs)==2)).sum()
        print("      bound head-frames=%d  (must be >0 or the probe never engages)" % bound)
        assert bound>0, "no heads ever bound"
print("SMOKE PASS")
PY
