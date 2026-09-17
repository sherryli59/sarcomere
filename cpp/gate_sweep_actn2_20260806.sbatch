#!/bin/bash
# Minimal-rig ACTN2 factorial gate: 4 actins + 2 myosins, ~35 s per run.
# {links 1,3} x {perp_all 0,1} x {rupture 0,0.12} x {max_bonds 2,3} x 6 seeds = 96
# ref_links is tied to links so the pair-edge stiffness is held fixed.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --job-name=aa_gate
#SBATCH --no-requeue
#SBATCH --array=0-95
#SBATCH --output=/scratch/users/xsli/sarcomere/gate_sweep_20260806/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/gate_sweep_20260806/logs/%A_%a.err

set -uo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

BIN=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/motor_validation
ROOT=/scratch/users/xsli/sarcomere/gate_sweep_20260806
mkdir -p "${ROOT}"/{data,logs}

ID=${SLURM_ARRAY_TASK_ID}
SEEDS=(20260805 20260806 20260807 20260808 20260809 20260810)
LINKS_OPT=(1 3); PERP_OPT=(0 1); RUP_OPT=(0 0.12); BONDS_OPT=(2 3)

SEED=${SEEDS[$(( ID % 6 ))]}
BONDS=${BONDS_OPT[$(( (ID / 6) % 2 ))]}
RUP=${RUP_OPT[$(( (ID / 12) % 2 ))]}
PERP=${PERP_OPT[$(( (ID / 24) % 2 ))]}
LINKS=${LINKS_OPT[$(( (ID / 48) % 2 ))]}

TAG="L${LINKS}_P${PERP}_R${RUP}_B${BONDS}_s${SEED}"
MS=60
SAVE=200

MV_AA_LINKS=${LINKS} MV_AA_REF_LINKS=${LINKS} \
MV_AA_PERP_ALL=${PERP} MV_AA_RUPTURE=${RUP} MV_AA_MAX_BONDS=${BONDS} \
"${BIN}" balanced ${MS} ${SEED} 35 0 0.5 -1 0 0.16 0.010 100 3 1 \
    "${ROOT}/data/${TAG}.csv" "${ROOT}/data/${TAG}.h5" ${SAVE} \
    100 0.64 0.80 5000 0.05 3 0.95 1 100 120 0
echo "TAG=${TAG} rc=$?"
