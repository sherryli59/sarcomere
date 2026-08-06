#!/bin/bash
#SBATCH --job-name=sarc_rng_ab
#SBATCH --partition=owners
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/users/xsli/sarcomere/perf_counter_dynamic_20260723/ab_%j.out

set -euo pipefail

OLD_EXE=/scratch/users/xsli/sarcomere/perf_openmp_exact_20260723/bin/sarcomere_exact_optimized_v2
NEW_EXE=/scratch/users/xsli/sarcomere/perf_counter_dynamic_20260723/bin/sarcomere_counter_dynamic_fused
JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
RESULT_ROOT=/scratch/users/xsli/sarcomere/perf_counter_dynamic_20260723/ab_${SLURM_JOB_ID}

mkdir -p "${RESULT_ROOT}/logs"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export LD_PRELOAD="${JEMALLOC}"

csv="${RESULT_ROOT}/timings.csv"
printf 'build,repeat,threads,wall_s,user_s,sys_s,max_rss_kb,exit_code\n' >"${csv}"

run_one() {
    local build=$1 exe=$2 repeat=$3 threads=$4
    local output="${RESULT_ROOT}/${build}_r${repeat}_t${threads}.h5"
    local timing="${RESULT_ROOT}/${build}_r${repeat}_t${threads}.time"
    local log="${RESULT_ROOT}/logs/${build}_r${repeat}_t${threads}.log"
    local status=0

    /usr/bin/time -f '%e,%U,%S,%M' -o "${timing}" \
        env OMP_NUM_THREADS="${threads}" OMP_SCHEDULE=static \
        "${exe}" \
        --n_actins=1280 --n_myosins=640 \
        --Lx=7.72 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
        --initial_structure=myosin_x_noise \
        --pre_phase_steps=100 --interaction_ramp_steps=2000 \
        --physical_mobility=true --temperature_kelvin=310 \
        --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000 \
        --k_aa=3000 --kappa_aa=50 \
        --dt=1e-6 --nsteps=200 --save_every=1000 \
        --seed=1731 --filename="${output}" >"${log}" 2>&1 || status=$?

    printf '%s,%s,%s,%s,%s\n' \
        "${build}" "${repeat}" "${threads}" \
        "$(tr '\n' ' ' <"${timing}")" "${status}" >>"${csv}"
}

# Interleave builds and thread counts so both see nearly the same node state.
for repeat in 1 2; do
    run_one old "${OLD_EXE}" "${repeat}" 1
    run_one new "${NEW_EXE}" "${repeat}" 1
    run_one old "${OLD_EXE}" "${repeat}" 8
    run_one new "${NEW_EXE}" "${repeat}" 8
done

python - "${csv}" "${RESULT_ROOT}/summary.csv" <<'PY'
import csv
import statistics
import sys
from collections import defaultdict

source, target = sys.argv[1:3]
groups = defaultdict(list)
with open(source, newline="") as stream:
    for row in csv.DictReader(stream):
        if int(row["exit_code"]) == 0:
            groups[(row["build"], int(row["threads"]))].append(float(row["wall_s"]))

with open(target, "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("build", "threads", "n", "median_wall_s", "speedup"))
    medians = {key: statistics.median(values) for key, values in groups.items()}
    for key in sorted(medians):
        build, threads = key
        writer.writerow((build, threads, len(groups[key]), medians[key],
                         medians[(build, 1)] / medians[key]))
PY

cat "${RESULT_ROOT}/summary.csv"
