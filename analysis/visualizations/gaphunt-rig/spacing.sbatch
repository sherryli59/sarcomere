#!/bin/bash
#SBATCH --job-name=spacing
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=1:00:00
#SBATCH --output=/scratch/users/xsli/sarcomere/gaphunt_20261023/logs/%x-%j.out
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin/python spacing_dist.py
