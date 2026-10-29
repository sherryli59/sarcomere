#!/bin/bash
#SBATCH --job-name=stallsmoke
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/stallforce_20261005/logs/smoke_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/stallforce_20261005/logs/smoke_%j.err
set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
root=/scratch/users/xsli/sarcomere/stallforce_20261005
exe=/scratch/users/xsli/sarcomere/barbed_capsule_normal_20260930/bin/sarcomere_barbed_capsule_normal
mkdir -p "$root"/{data,logs,smoke}; test -x "$exe"
# Paired short runs, identical except am_stall_force, so the RETENTION confound
# is visible: if bound heads collapse at 1.5 the arm cannot be interpreted.
for sf in 4.5 1.5; do
  srun --kill-on-bad-exit=1 "$exe" \
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
    --actin_length=1 --myosin_length=1.5 --actin_radius=0.015 --myosin_radius=0.05 \
    --initial_structure=myosin_x_noise --actin_direction_noise=0.40 \
    --myosin_direction_noise=0.40 --am_head_azimuth=false --am_spatial_heads=true \
    --am_heads_per_half=27 --am_spokes_per_crown=3 \
    --am_max_heads_per_actin_half=24 --am_max_heads_per_actin=24 \
    --dynamic_am=true --kappa_aa=3 --k_aa=1000 --am_k_strong_off=160 \
    --am_motor_k=1000 --am_power_stroke=0.005 --am_stall_force="$sf" \
    --dt=1e-5 --nsteps=6000 --save_every=500 --seed=20260805 \
    --filename="$root/smoke/sf${sf}.h5" >"$root/logs/smoke_sf${sf}.log" 2>&1
  echo "ran stall=$sf"
done
python3 - <<'PY'
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py
R="/scratch/users/xsli/sarcomere/stallforce_20261005/smoke"
res={}
for sf in ("4.5","1.5"):
    with h5py.File(f"{R}/sf{sf}.h5","r",locking=False) as h:
        a=h["meta"].attrs
        hs=np.asarray(h["state/am_head_state_current"])
        bound=int(((hs==1)|(hs==2)).sum()); strong=int((hs==2).sum())
        links=float(np.asarray(h["state/aa_pairs_current_count"])[-1].reshape(-1)[-1])
        res[sf]=(float(a["am_stall_force"]), bound, strong, links, hs.shape)
        print("stall=%-4s attr=%-4.1f  bound head-frames=%-8d strong=%-8d ACTN2 last=%.0f  head_axis=%d"
              % (sf, a["am_stall_force"], bound, strong, links, hs.shape[2]))
        assert abs(float(a["am_stall_force"])-float(sf))<1e-9, "stall force did not reach the binary"
        assert int(a["am_heads_per_half"])==27 and hs.shape[2]==54
b45,b15 = res["4.5"][1], res["1.5"][1]
s45,s15 = res["4.5"][2], res["1.5"][2]
print()
print("RETENTION CHECK  bound %d -> %d (%.2fx)   strong %d -> %d (%.2fx)"
      % (b45,b15,b15/b45, s45,s15,s15/s45))
assert b15 > 0.5*b45, "bound heads collapsed >2x: arm would be confounded by retention loss"
print("SMOKE PASS: stall force reaches the binary, hph=27 intact, retention not collapsed")
PY
