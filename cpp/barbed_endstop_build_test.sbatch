#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --job-name=barbed_endstop_test
#SBATCH --output=/scratch/users/xsli/sarcomere/barbed_endstop_20260930/build_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/barbed_endstop_20260930/build_%j.err
set -euo pipefail

export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
export OMP_NUM_THREADS=1
repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere_am_kmc/.codex_work/barbed_endstop
root=/scratch/users/xsli/sarcomere/barbed_endstop_20260930
stage=${root}/stage
build=${root}/build
mkdir -p "${root}"
rm -rf "${stage}" "${build}"
mkdir -p "${stage}/cpp"
for item in src include autodiff tests CMakeLists.txt; do
    cp -a "${repo}/cpp/${item}" "${stage}/cpp/"
done
cmake -S "${stage}/cpp" -B "${build}" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin/g++
cmake --build "${build}" -j 8 --target sarcomere \
  catch_bond_kinetics_test zdisc_myosin_tip_contact_test \
  actn2_am_occlusion_test
ctest --test-dir "${build}" --output-on-failure \
  -R 'catch_bond_kinetics|zdisc_myosin_tip_contact|actn2_am_occlusion'
cp "${build}/sarcomere" "${root}/sarcomere_barbed_endstop"
