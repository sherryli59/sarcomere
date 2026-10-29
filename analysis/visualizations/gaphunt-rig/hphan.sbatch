#!/bin/bash
#SBATCH --job-name=hphan
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G --time=2:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/hphan-%j.out
export OMP_NUM_THREADS=2 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
python3 hph_an.py 27:20260805,54:20260815,81:20260805 400 hph_test.json
echo TEST_DONE
python3 hph_an.py 27:20260805,54:20260815,81:20260805 20 hph_an3.json
echo HPHAN_DONE
