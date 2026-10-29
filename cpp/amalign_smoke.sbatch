#!/bin/bash
#SBATCH --job-name=alsmoke
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/amalign_20261005/logs/smoke_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/amalign_20261005/logs/smoke_%j.err
set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
root=/scratch/users/xsli/sarcomere/amalign_20261005
exe=/scratch/users/xsli/sarcomere/barbed_capsule_normal_20260930/bin/sarcomere_barbed_capsule_normal
mkdir -p "$root"/{data,logs,smoke}; test -x "$exe"
# At am_alignment_k=0 the tb_am_align channel is EXACTLY 0.000 (measured). So a
# nonzero tb_am_align is proof the term actually activated, not just that an
# attribute was recorded.
for spec in "9:0" "9:30" "27:0" "27:30"; do
  IFS=: read -r hph k <<< "$spec"
  srun --kill-on-bad-exit=1 "$exe" \
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
    --actin_length=1 --myosin_length=1.5 --actin_radius=0.015 --myosin_radius=0.05 \
    --initial_structure=myosin_x_noise --actin_direction_noise=0.40 \
    --myosin_direction_noise=0.40 --am_head_azimuth=false --am_spatial_heads=true \
    --am_heads_per_half="$hph" --am_spokes_per_crown=3 \
    --am_max_heads_per_actin_half=24 --am_max_heads_per_actin=24 \
    --dynamic_am=true --kappa_aa=3 --k_aa=1000 --am_k_strong_off=160 \
    --am_alignment_k="$k" --am_alignment_slack_deg=2 --am_tether_saturation_heads=6 \
    --torque_breakdown=true \
    --dt=1e-5 --nsteps=6000 --save_every=500 --seed=20260805 \
    --filename="$root/smoke/h${hph}_k${k}.h5" >"$root/logs/smoke_h${hph}_k${k}.log" 2>&1
  echo "ran hph=$hph k=$k"
done
python3 - <<'PY'
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
import h5py
R="/scratch/users/xsli/sarcomere/amalign_20261005/smoke"
for hph in (9,27):
    for k in (0,30):
        with h5py.File(f"{R}/h{hph}_k{k}.h5","r",locking=False) as h:
            a=h["meta"].attrs
            tb=np.stack([np.asarray(h[f"actin/tb_am_align_{c}"],float) for c in "xyz"],-1)
            m=np.linalg.norm(tb,axis=-1); nz=m[m>0]
            print("hph=%-3d k=%-3d attr_k=%-5.1f  |tb_am_align| nonzero=%-8d mean(nonzero)=%.4f pN.um  S=%s"
                  % (hph,k,a["am_alignment_k"],nz.size,nz.mean() if nz.size else 0.0,""))
            assert abs(float(a["am_alignment_k"])-k)<1e-9, "k did not reach the binary"
            if k==0: assert nz.size==0, "tb_am_align must be exactly 0 at k=0, got %d nonzero" % nz.size
            else:    assert nz.size>0, "tb_am_align is still zero at k=30: the term did NOT activate"
print("SMOKE PASS: the alignment term activates at k=30 and is exactly inert at k=0, both head counts")
PY
