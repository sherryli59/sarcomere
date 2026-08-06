#!/bin/bash
#SBATCH -p normal,hns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=cg_fields
#SBATCH --open-mode=append
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_cg.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_cg.out

# Coarse-grain a representative subset and measure q_dom(t), period, and
# transverse coherence. Spans:
#   - the catch-bond strength axis (lt0 .. lt1e0) at fixed A=1200  -> does the
#     pattern depend on the directional catch bond?
#   - the density axis (A1200 vs A1500) at fixed lt               -> is the band
#     spacing steric (density-independent) or an emergent wavelength?

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

DATA=/scratch/users/xsli/sarcomere/data0430
OUT=/scratch/users/xsli/sarcomere/cg_fields0430
SCRIPT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/coarse_grain_fields.py
mkdir -p "$OUT"

# catch-bond strength sweep, A=1200, replicate r00
LT_FILES=(
  "$DATA/A1200_M1000_lt0_r00_L12x2x2.h5"
  "$DATA/A1200_M1000_lt1e-3_r00_L12x2x2.h5"
  "$DATA/A1200_M1000_lt1e-2_r00_L12x2x2.h5"
  "$DATA/A1200_M1000_lt1e-1_r00_L12x2x2.h5"
  "$DATA/A1200_M1000_lt1e0_r00_L12x2x2.h5"
)
# density contrast at fixed lt (add A1500 partners)
DENS_FILES=(
  "$DATA/A1500_M1000_lt0_r00_L12x2x2.h5"
  "$DATA/A1500_M1000_lt1e-2_r00_L12x2x2.h5"
)

python "$SCRIPT" "${LT_FILES[@]}" "${DENS_FILES[@]}" \
    --out_dir "$OUT" \
    --stride 5 \
    --nbins 240 \
    --smooth 0.05 \
    --n_snapshots 6 \
    --population cb2

echo "done; outputs in $OUT"
