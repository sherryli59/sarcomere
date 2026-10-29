#!/bin/bash
#SBATCH --partition=normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --job-name=am_partner_dist
#SBATCH --output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/analysis/am_partner_dist_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/analysis/am_partner_dist_%j.err
set -euo pipefail

repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
baseline=/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data
ksoff=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/data
output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/analysis
mkdir -p "${output}"

python "${repo}/analysis/analyze_am_partner_distribution.py" \
  --pair "20260805=${baseline}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260805_5000ms_capsule_normal.h5,${ksoff}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260805_5000ms_ksoff27.h5" \
  --pair "20260815=${baseline}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260815_5000ms_capsule_normal.h5,${ksoff}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260815_5000ms_ksoff27.h5" \
  --pair "20260825=${baseline}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260825_5000ms_capsule_normal.h5,${ksoff}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260825_5000ms_ksoff27.h5" \
  --pair "20260835=${baseline}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260835_5000ms_capsule_normal.h5,${ksoff}/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed20260835_5000ms_ksoff27.h5" \
  --burn-in 0.5 --tail 0.5 \
  --output-json "${output}/am_partner_distribution.json" \
  --output-csv "${output}/am_partner_distribution.csv"
