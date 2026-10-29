#!/bin/bash
#SBATCH --job-name=kaa3seedviz
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/kaa3seedviz_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/kaa3seedviz_%j.err

# Viewer for the kappa_aa=3 capsule-normal production run at seed 20260805 --
# the run whose S(t) course reaches 0.7954 at 4200 ms -- paired against seed
# 20260835 of the SAME arm, so the seed-to-seed difference is visible. Same
# campaign, same config, same binary: only the seed differs.
set -euo pipefail
echo "host=$(hostname)  start=$(date '+%F %T')"
export HDF5_USE_FILE_LOCKING=FALSE
VDIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/ksoff27-prod
OUTDIR=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/analysis
mkdir -p "$OUTDIR"
D=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data
S805=$D/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260805_5000ms_capsule_normal.h5
S835=$D/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260835_5000ms_capsule_normal.h5
for f in "$S805" "$S835" "$VDIR/build-ksoff27-prod.py" "$VDIR/ksoff27-prod-template.html"; do
    [ -r "$f" ] || { echo "ABORT: missing $f"; exit 1; }
done
python3 "$VDIR/build-ksoff27-prod.py" "$S805" "$S835" \
    --labels "kappa_aa=3 seed 20260805 (S=0.795 @ 4200 ms)" "kappa_aa=3 seed 20260835" \
    --stride 50 \
    --title "kappa_aa=3 capsule-normal, k_soff=160: seed 20260805 vs 20260835 (seed scatter)" \
    --template "$VDIR/ksoff27-prod-template.html" \
    --out "$OUTDIR/kaa3-seed-scatter-traj.html"
echo "--- result"; ls -la "$OUTDIR/kaa3-seed-scatter-traj.html"
python3 -c "
s=open('$OUTDIR/kaa3-seed-scatter-traj.html',encoding='ascii').read()
print('ASCII-safe:', s.isascii())
print('markers replaced:', '/*__RUNS__*/' not in s and '__TITLE__' not in s)
print('stale minimal-rig text absent:', 'minimal_catch_balanced' not in s)
"
echo "end=$(date '+%F %T')"
