#!/bin/bash
#SBATCH --job-name=chain
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G --time=2:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/chain-%j.out
export OMP_NUM_THREADS=2 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 chain_an.py
echo CHAIN_DONE
