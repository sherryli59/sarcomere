#!/bin/bash
#SBATCH --job-name=hphsmoke
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/hph_ladder_20261005/logs/smoke_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/hph_ladder_20261005/logs/smoke_%j.err
set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
root=/scratch/users/xsli/sarcomere/hph_ladder_20261005
exe=/scratch/users/xsli/sarcomere/barbed_capsule_normal_20260930/bin/sarcomere_barbed_capsule_normal
mkdir -p "$root"/{data,logs,smoke}
test -x "$exe"
# Only the head-geometry knobs are exercised; everything else is irrelevant to
# the question "does the knob reach the binary and produce 9 crowns".
for spec in "spokes1:9:1" "hph27:27:3"; do
  IFS=: read -r lab hph spk <<< "$spec"
  srun --kill-on-bad-exit=1 "$exe" \
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
    --actin_length=1 --myosin_length=1.5 --actin_radius=0.015 --myosin_radius=0.05 \
    --initial_structure=myosin_x_noise --actin_direction_noise=0.40 \
    --myosin_direction_noise=0.40 --am_head_azimuth=false --am_spatial_heads=true \
    --am_heads_per_half="$hph" --am_spokes_per_crown="$spk" \
    --am_max_heads_per_actin_half=24 --am_max_heads_per_actin=24 \
    --dynamic_am=true --kappa_aa=3 --k_aa=1000 --am_k_strong_off=160 \
    --dt=1e-5 --nsteps=3000 --save_every=500 --seed=20260805 \
    --filename="$root/smoke/smoke_${lab}.h5" >"$root/logs/smoke_${lab}.log" 2>&1
  echo "ran $lab"
done
python3 - <<'PY'
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py
BEAR, Lm = 0.45, 1.5
R="/scratch/users/xsli/sarcomere/hph_ladder_20261005/smoke"
for lab, hph, spk in (("spokes1",9,1), ("hph27",27,3)):
    with h5py.File(f"{R}/smoke_{lab}.h5","r",locking=False) as h:
        a=h["meta"].attrs
        got_h=int(a["am_heads_per_half"]); got_s=int(float(a.get("am_spokes_per_crown",3)))
        axis=h["state/am_head_state_current"].shape[2]
        cc=max(1,(got_h+got_s-1)//got_s)
        print("%-8s attr hph=%-3d spokes=%d  head_axis=%-3d crowns/half=%-2d spacing=%5.1f nm  heads/myo=%d"
              % (lab, got_h, got_s, axis, cc, BEAR*Lm/cc*1000, 2*got_h))
        assert got_h==hph, "hph did not reach the binary"
        assert got_s==spk, "spokes did not reach the binary"
        assert axis==2*hph, "head axis %d != 2*%d" % (axis,hph)
        assert cc==9, "expected 9 crowns/half, got %d" % cc
        assert abs(float(a["box_x"])-6.0)<1e-9, "Lx must be 6"
        bound=((np.asarray(h["state/am_head_state_current"])==1)|
               (np.asarray(h["state/am_head_state_current"])==2)).sum()
        aa=float(np.asarray(h["state/aa_pairs_current_count"])[-1].reshape(-1)[-1])
        print("         bound head-frames=%-8d ACTN2 pairs at last frame=%.0f" % (bound, aa))
        assert bound>0 and aa>0, "nothing engaged"
print("SMOKE PASS: both arms give 9 crowns/half at 75 nm with the intended head budgets")
PY
