#!/bin/bash
#SBATCH --job-name=connector-perf
#SBATCH --partition=normal
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
set -euo pipefail
export CONDA_PREFIX="$GROUP_HOME/sherry/mamba_root/envs/sarcomere"
export PATH="$CONDA_PREFIX/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
task_build="$SCRATCH/sarcomere/connector_perf_20261028"
task_scripts="$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/status2_gate_20261028"
cmake -S "$SCRATCH/sarcomere/zdisclink_20261020/cpp" -B "$task_build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON > "$task_build/configure.log" 2>&1
cmake --build "$task_build" -j "$SLURM_CPUS_PER_TASK" > "$task_build/build.log" 2>&1
cd "$task_build"
ctest --output-on-failure -j 1 > ctest.log 2>&1
cp "$task_scripts/reference_geometry_before_perf.h" "$task_build/reference_geometry.h"
g++ -std=c++17 -O3 -UNDEBUG -fopenmp -I"$SCRATCH/sarcomere/zdisclink_20261020/cpp/include" -I"$task_build" "$task_scripts/deletion_benchmark.cpp" -o deletion_benchmark
./deletion_benchmark "$SCRATCH/sarcomere/audit_20261028" > deletion_benchmark.json
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
python -u "$task_scripts/benchmark_steps.py" "$task_build"
export OMP_NUM_THREADS=1
python -u "$task_scripts/connector_metrics.py" --output "$task_build/connector_metrics.json"
sha256sum sarcomere > binary.sha256
cat deletion_benchmark.json
tail -4 ctest.log
