#!/bin/bash
#SBATCH --job-name=ksoff27viz
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/ksoff27viz_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/ksoff27viz_%j.err

# Interactive trajectory viewer for the k_soff=27 PRODUCTION run, paired
# seed-for-seed against its k_soff=160 baseline. Run as a job, not on a login
# node: the builder reads the full ACTN2 link arrays out of two ~14 GB files,
# which is exactly the intensive I/O /etc/agents/AGENTS.md says must be a job.
set -euo pipefail
echo "host=$(hostname)  start=$(date '+%F %T')"

export HDF5_USE_FILE_LOCKING=FALSE
VDIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/ksoff27-prod
OUTDIR=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/analysis
mkdir -p "$OUTDIR"

K27=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260835_5000ms_ksoff27.h5
K160=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260835_5000ms_capsule_normal.h5

for f in "$K27" "$K160" "$VDIR/build-ksoff27-prod.py" "$VDIR/ksoff27-prod-template.html"; do
    [ -r "$f" ] || { echo "ABORT: missing $f"; exit 1; }
done

python3 "$VDIR/build-ksoff27-prod.py" \
    "$K160" "$K27" \
    --labels "k_soff=160 baseline (seed 835)" "k_soff=27 (seed 835)" \
    --stride "${STRIDE:-50}" \
    --title "k_soff = 27 vs 160, A1600/M320 Lx=6, kappa_aa=3, seed 20260835" \
    --template "$VDIR/ksoff27-prod-template.html" \
    --out "$OUTDIR/ksoff27-prod-traj-stride${STRIDE:-50}.html"

echo "--- result"
ls -la "$OUTDIR/ksoff27-prod-traj-stride${STRIDE:-50}.html"
python3 -c "
s=open('$OUTDIR/ksoff27-prod-traj-stride${STRIDE:-50}.html',encoding='ascii').read()
print('ASCII-safe:', s.isascii())
print('no unreplaced markers:', '/*__RUNS__*/' not in s and '__TITLE__' not in s)
"
echo "end=$(date '+%F %T')"
