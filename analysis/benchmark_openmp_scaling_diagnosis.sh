#!/bin/bash
set -euo pipefail

REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
EXE="${SARC_BENCH_EXE:-${REPO}/cpp/build/sarcomere}"
JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
ROOT="${SARC_DIAG_ROOT:-/scratch/users/xsli/sarcomere/perf_openmp_diagnosis_20260723}"
mkdir -p "${ROOT}/logs"
exec 9>"${ROOT}/benchmark.lock"
flock -n 9 || exit 0

export HDF5_USE_FILE_LOCKING=FALSE
export OMP_DYNAMIC=FALSE
scratch_dir=$(mktemp -d /tmp/sarc_omp_diag.XXXXXX)
trap 'rm -rf "${scratch_dir}"' EXIT
csv="${ROOT}/timings.csv"
printf 'test,repeat,threads,schedule,n_actins,n_myosins,pre_steps,nsteps,wall_s,user_s,sys_s,max_rss_kb,exit_code\n' >"${csv}"

run_one() {
    local test_name=$1 repeat=$2 threads=$3 schedule=$4
    local n_actins=$5 n_myosins=$6 lx=$7 ly=$8 lz=$9
    local pre_steps=${10} nsteps=${11}
    local label="${test_name}_r${repeat}_t${threads}_${schedule//[,]/_}"
    local output="${scratch_dir}/${label}.h5"
    local timing="${scratch_dir}/${label}.time"
    local log="${ROOT}/logs/${label}.log"
    local status=0
    /usr/bin/time -f '%e,%U,%S,%M' -o "${timing}" \
        env OMP_NUM_THREADS="${threads}" OMP_SCHEDULE="${schedule}" \
        OPENBLAS_NUM_THREADS=1 LD_PRELOAD="${JEMALLOC}" \
        "${EXE}" \
        --n_actins="${n_actins}" --n_myosins="${n_myosins}" \
        --Lx="${lx}" --Ly="${ly}" --Lz="${lz}" --dimension=3 --pbc_mask=111 \
        --initial_structure=myosin_x_noise \
        --pre_phase_steps="${pre_steps}" --interaction_ramp_steps=2000 \
        --physical_mobility=true --temperature_kelvin=310 \
        --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000 \
        --k_aa=3000 --kappa_aa=50 \
        --dt=1e-6 --nsteps="${nsteps}" --save_every=1000 \
        --seed=1731 --filename="${output}" >"${log}" 2>&1 || status=$?
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${test_name}" "${repeat}" "${threads}" "${schedule}" \
        "${n_actins}" "${n_myosins}" "${pre_steps}" "${nsteps}" \
        "$(tr '\n' ' ' <"${timing}")" "${status}" >>"${csv}"
}

for repeat in 1 2 3; do
    for schedule in static dynamic,1 dynamic,8 guided,1; do
        run_one schedule "${repeat}" 8 "${schedule}" 320 160 7.72 1 1 200 400
    done
done

for repeat in 1 2; do
    for threads in 1 2 4 8; do
        run_one production_scale "${repeat}" "${threads}" static \
            1280 640 7.72 2 2 100 200
    done
done

python - "${csv}" "${ROOT}/summary.csv" <<'PY'
import csv
import statistics
import sys
from collections import defaultdict

source, target = sys.argv[1:3]
groups = defaultdict(list)
with open(source, newline="") as stream:
    for row in csv.DictReader(stream):
        if int(row["exit_code"]) == 0:
            key = (row["test"], int(row["threads"]), row["schedule"])
            groups[key].append(float(row["wall_s"]))

with open(target, "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("test", "threads", "schedule", "n", "median_wall_s", "min_wall_s"))
    for key in sorted(groups):
        values = groups[key]
        writer.writerow((*key, len(values), statistics.median(values), min(values)))
PY
