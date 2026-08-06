#!/bin/bash
set -euo pipefail

REPO=/home/groups/rotskoff/sherry/sarcomere/sarcomere
EXE="${SARC_BENCH_EXE:-${REPO}/cpp/build/sarcomere}"
JEMALLOC="${SARC_JEMALLOC:-/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so}"
ROOT="${SARC_BENCH_ROOT:-/scratch/users/xsli/sarcomere/perf_allocator_blas_20260723}"
REPEATS="${SARC_BENCH_REPEATS:-3}"
NSTEPS="${SARC_BENCH_NSTEPS:-400}"
mkdir -p "${ROOT}/logs" "${ROOT}/perf"

test -x "${EXE}"
test -r "${JEMALLOC}"
exec 9>"${ROOT}/benchmark.lock"
flock -n 9 || {
    echo "Another allocator/BLAS benchmark owns ${ROOT}/benchmark.lock"
    exit 0
}
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static

scratch_dir=$(mktemp -d /tmp/sarc_allocator_blas.XXXXXX)
trap 'rm -rf "${scratch_dir}"' EXIT

csv="${ROOT}/timings.csv"
printf 'repeat,threads,condition,wall_s,user_s,sys_s,max_rss_kb,exit_code\n' >"${csv}"

run_one() {
    local repeat=$1
    local threads=$2
    local condition=$3
    local output="${scratch_dir}/${condition}_t${threads}_r${repeat}.h5"
    local timing="${scratch_dir}/${condition}_t${threads}_r${repeat}.time"
    local log="${ROOT}/logs/${condition}_t${threads}_r${repeat}.log"
    local -a runtime_env=(
        env
        -u OPENBLAS_NUM_THREADS
        -u GOTO_NUM_THREADS
        -u MKL_NUM_THREADS
        -u LD_PRELOAD
        "OMP_NUM_THREADS=${threads}"
    )
    case "${condition}" in
        baseline)
            ;;
        blas1)
            runtime_env+=("OPENBLAS_NUM_THREADS=1")
            ;;
        jemalloc)
            runtime_env+=("LD_PRELOAD=${JEMALLOC}")
            ;;
        both)
            runtime_env+=("OPENBLAS_NUM_THREADS=1" "LD_PRELOAD=${JEMALLOC}")
            ;;
        *)
            return 2
            ;;
    esac

    local status=0
    /usr/bin/time -f '%e,%U,%S,%M' -o "${timing}" \
        "${runtime_env[@]}" "${EXE}" \
        --n_actins=320 --n_myosins=160 \
        --Lx=7.72 --Ly=1 --Lz=1 --dimension=3 --pbc_mask=111 \
        --initial_structure=myosin_x_noise \
        --pre_phase_steps=200 --interaction_ramp_steps=2000 \
        --physical_mobility=true --temperature_kelvin=310 \
        --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000 \
        --k_aa=3000 --kappa_aa=50 \
        --dt=1e-6 --nsteps="${NSTEPS}" --save_every=1000 \
        --seed=1731 --filename="${output}" >"${log}" 2>&1 || status=$?
    printf '%s,%s,%s,%s,%s\n' \
        "${repeat}" "${threads}" "${condition}" "$(tr '\n' ' ' <"${timing}")" "${status}" \
        >>"${csv}"
    if [ "${status}" -ne 0 ]; then
        return "${status}"
    fi
}

conditions=(baseline blas1 jemalloc both)
threads=(1 2 4 8)
for repeat in $(seq 1 "${REPEATS}"); do
    for thread_count in "${threads[@]}"; do
        offset=$(( (repeat + thread_count) % 4 ))
        for index in 0 1 2 3; do
            condition="${conditions[$(( (index + offset) % 4 ))]}"
            run_one "${repeat}" "${thread_count}" "${condition}"
        done
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
            groups[(row["condition"], int(row["threads"]))].append(float(row["wall_s"]))

baseline_one = statistics.median(groups[("baseline", 1)])
with open(target, "w", newline="") as stream:
    fields = (
        "condition", "threads", "n", "median_wall_s", "min_wall_s",
        "speedup_vs_baseline_1t", "parallel_efficiency",
    )
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for condition in ("baseline", "blas1", "jemalloc", "both"):
        for threads in (1, 2, 4, 8):
            values = groups[(condition, threads)]
            median = statistics.median(values)
            speedup = baseline_one / median
            writer.writerow({
                "condition": condition,
                "threads": threads,
                "n": len(values),
                "median_wall_s": median,
                "min_wall_s": min(values),
                "speedup_vs_baseline_1t": speedup,
                "parallel_efficiency": speedup / threads,
            })
PY

profile_one() {
    local condition=$1
    local data="${ROOT}/perf/${condition}_omp8.data"
    local report="${ROOT}/perf/${condition}_omp8.txt"
    local output="${scratch_dir}/perf_${condition}.h5"
    local -a runtime_env=(
        env -u OPENBLAS_NUM_THREADS -u GOTO_NUM_THREADS
        -u MKL_NUM_THREADS -u LD_PRELOAD OMP_NUM_THREADS=8
    )
    if [ "${condition}" = both ]; then
        runtime_env+=("OPENBLAS_NUM_THREADS=1" "LD_PRELOAD=${JEMALLOC}")
    fi
    perf record -q -F 99 -g -o "${data}" -- \
        "${runtime_env[@]}" "${EXE}" \
        --n_actins=320 --n_myosins=160 \
        --Lx=7.72 --Ly=1 --Lz=1 --dimension=3 --pbc_mask=111 \
        --initial_structure=myosin_x_noise \
        --pre_phase_steps=200 --interaction_ramp_steps=2000 \
        --physical_mobility=true --temperature_kelvin=310 \
        --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000 \
        --k_aa=3000 --kappa_aa=50 \
        --dt=1e-6 --nsteps="${NSTEPS}" --save_every=1000 \
        --seed=1731 --filename="${output}" >/dev/null 2>&1
    perf report -i "${data}" --stdio --percent-limit 0.5 >"${report}"
}

if command -v perf >/dev/null 2>&1; then
    profile_one baseline || true
    profile_one both || true
fi

echo "timings=${csv}"
echo "summary=${ROOT}/summary.csv"
