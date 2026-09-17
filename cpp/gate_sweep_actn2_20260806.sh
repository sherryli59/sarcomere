#!/bin/bash
# Minimal-rig ACTN2 factorial gate.
#   4 actins + 2 myosins (motor_validation balanced), seconds per run.
#   {links 1,3} x {perp_all 0,1} x {rupture 0,0.12} x {max_bonds 2,3} x 6 seeds
#   = 96 runs. ref_links is tied to links so the pair-edge stiffness is held
#   fixed: N links at k/N and one link at k are the same bond.
set -uo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

BIN=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/motor_validation
ROOT=/scratch/users/xsli/sarcomere/gate_sweep_20260806
rm -rf "${ROOT}"; mkdir -p "${ROOT}"/{data,logs}
RESULTS=${ROOT}/results.tsv
printf "links\tperp_all\trupture\tmax_bonds\tseed\tstatus2_final\tstatus1_final\tgeo_frac\tdeaths\n" > "${RESULTS}"

MS=60          # 60 ms: past the 25-35 ms window where status-2 maturation starts
SAVE=200

for LINKS in 1 3; do
for PERP in 0 1; do
for RUP in 0 0.12; do
for BONDS in 2 3; do
for SEED in 20260805 20260806 20260807 20260808 20260809 20260810; do
  TAG="L${LINKS}_P${PERP}_R${RUP}_B${BONDS}_s${SEED}"
  MV_AA_LINKS=${LINKS} MV_AA_REF_LINKS=${LINKS} \
  MV_AA_PERP_ALL=${PERP} MV_AA_RUPTURE=${RUP} MV_AA_MAX_BONDS=${BONDS} \
  "${BIN}" balanced ${MS} ${SEED} 35 0 0.5 -1 0 0.16 0.010 100 3 1 \
      "${ROOT}/data/${TAG}.csv" "${ROOT}/data/${TAG}.h5" ${SAVE} \
      100 0.64 0.80 5000 0.05 3 0.95 1 100 120 0 \
      > "${ROOT}/logs/${TAG}.log" 2>&1
  echo "${TAG} rc=$?"
done; done; done; done; done

echo "sweep complete: $(ls ${ROOT}/data/*.h5 2>/dev/null | wc -l) h5 files"
