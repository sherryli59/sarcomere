#!/bin/bash
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:40:00
#SBATCH --job-name=besviz
#SBATCH --exclude=sh03-07n53,sh03-07n55,sh03-07n57,sh03-07n58,sh03-10n30,sh03-10n33,sh03-10n35
#SBATCH --output=/scratch/users/xsli/sarcomere/barbed_endstop_kappa_5s_20260930/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/barbed_endstop_kappa_5s_20260930/logs/%x_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}

# Reading three 764 MB trajectories is real I/O, so this runs as a job rather
# than on a login node (/etc/agents/AGENTS.md).
B=/scratch/users/xsli/sarcomere/barbed_endstop_kappa_5s_20260930/data
N=/scratch/users/xsli/sarcomere/novetrun_20260928/data
OUT=/scratch/users/xsli/sarcomere/barbed_endstop_kappa_5s_20260930/barbed-endstop-traj.html
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/minimal-wrap-animation

python build-minimal-wrap-animation.py \
  "$N/novetrun_vetoon_seed20260835.h5" \
  "$B/barbed_endstop_kappa0_seed20260835_5s.h5" \
  "$B/barbed_endstop_kappa3_seed20260835_5s.h5" \
  --stride 20 --box 4 1 1 \
  --out "$OUT" \
  --template minimal-wrap-animation-template.html \
  --title "Barbed-tip end-stop, seed 20260835: no end-stop vs end-stop kappa_aa=0 vs end-stop kappa_aa=3"
ls -la "$OUT"
