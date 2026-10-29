#!/bin/bash
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=prodbuild
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_steric080_20260930/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_steric080_20260930/logs/%x_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
export OMP_NUM_THREADS=1

repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_steric080_20260930}
mkdir -p "${root}"/{bin,logs,data,analysis,source_snapshot,src}

# Stage only the CMake inputs: cpp/ as a whole is 1.3 GB of test output.
echo "=== staging"
rm -rf "${root}/src"; mkdir -p "${root}/src/cpp"
for f in src include autodiff tests CMakeLists.txt; do cp -a "${repo}/cpp/${f}" "${root}/src/cpp/"; done
cp "${repo}/analysis/monitor_live_production_hdf5.py" "${root}/source_snapshot/"
cp "${repo}/analysis/check_live_production_hdf5.py" "${root}/source_snapshot/"
cp -a "${root}/src/cpp" "${root}/source_snapshot/"
( cd "${repo}" && git rev-parse HEAD && git status --short | grep -vE '^\?\?' ) \
    > "${root}/source_snapshot/GIT_STATE.txt" 2>&1 || true
( cd "${repo}" && git diff ) > "${root}/source_snapshot/WORKING_TREE.patch" 2>&1 || true
staged_mb=$(du -sm "${root}/src" | awk '{print $1}')
test "${staged_mb}" -lt 200 || { echo "FAIL staged ${staged_mb} MB"; exit 1; }
echo "  staged ${staged_mb} MB"

echo "=== build"
cmake -S "${root}/src/cpp" -B "${root}/build" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin/g++ \
      > "${root}/logs/cmake.log" 2>&1
cmake --build "${root}/build" --target sarcomere -j 8 > "${root}/logs/build.log" 2>&1
cp "${root}/build/sarcomere" "${root}/bin/sarcomere_steric080"
sha256sum "${root}/bin/sarcomere_steric080" | tee "${root}/source_snapshot/BINARY_SHA256.txt"

echo "=== tests that cover the changed code"
for t in actn2_am_occlusion_test zdisc_myosin_tip_contact_test catch_bond_kinetics_test \
         barbed_lateral_test pbc_round_test; do
    cmake --build "${root}/build" --target "$t" -j 8 >> "${root}/logs/build.log" 2>&1
    echo "--- $t"
    "${root}/build/$t"
done

echo "=== the new code path must be reachable and recorded"
"${root}/bin/sarcomere_steric080" --help 2>&1 | grep -c "aa_overlap_occlusion_veto" \
    | xargs -I{} sh -c 'test {} -ge 1 && echo "  flag aa_overlap_occlusion_veto present in --help"'
strings "${root}/bin/sarcomere_steric080" | grep -q "compute_barbed_tip_myosin_endstop\|aa_barbed_tip_endstop_force" \
    && echo "  barbed-tip end-stop present in the binary"

echo "ALL BUILD GATES PASSED"
