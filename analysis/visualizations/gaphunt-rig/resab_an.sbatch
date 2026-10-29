#!/bin/bash
#SBATCH --job-name=resaban
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G --time=3:00:00
#SBATCH --output=/scratch/users/xsli/sarcomere/gaphunt_20261023/logs/%x-%j.out
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin/python resab_an.py /scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/resab_an.json
