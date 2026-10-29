#!/bin/bash
#SBATCH --job-name=trisviz
#SBATCH --partition=hns,normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=1:00:00
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/trisbuild-%j.out
set -u
export OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig
grep -h -o -- '--aa_zdisc_tip_contact_requires_link=[a-z]*\|--n_myosins=[0-9]*\|--n_fixed_myosins=[0-9]*\|--pre_phase_steps=[0-9]*' /scratch/users/xsli/sarcomere/gaphunt_20261023/logs/tris30_seed20260805.argv
python3 tris_an.py
python3 build-gaphunt-rig.py --groups 'One actin shared by 3 TILTED myosins (sigma 0.40)' --title 'One actin, three tilted myosins' --out tris.html
python3 validate-minimal-wrap.py tris.html | grep -E '^PASS|FAIL'
echo TRISBUILD_DONE
