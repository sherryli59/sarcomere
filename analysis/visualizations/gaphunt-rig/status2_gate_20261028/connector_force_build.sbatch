#!/bin/bash
#SBATCH --job-name=connector-force-check
#SBATCH --partition=normal
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
set -euo pipefail
export CONDA_PREFIX="$GROUP_HOME/sherry/mamba_root/envs/sarcomere"
export PATH="$CONDA_PREFIX/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
task_build="$SCRATCH/sarcomere/connector_force_20261028"
mkdir -p "$task_build"
cmake -S "$SCRATCH/sarcomere/zdisclink_20261020/cpp" -B "$task_build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON > "$task_build/configure.log" 2>&1
cmake --build "$task_build" -j "$SLURM_CPUS_PER_TASK" > "$task_build/build.log" 2>&1
cd "$task_build"
ctest --output-on-failure -j 1 > ctest.log 2>&1
sha256sum sarcomere > binary.sha256
cat ctest.log
