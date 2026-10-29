#!/bin/bash
#SBATCH --job-name=probehph
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --no-requeue
#SBATCH --array=0-5
#SBATCH --output=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1     # 2 filaments; 1 thread also makes runs reproducible

# ===========================================================================
# CONFINED ALIGNMENT PROBE -- does raising am_heads_per_half unlock the span
# couple at a MISALIGNED start?
#
# WHY. The span couple is established (confined probe v2, job 42764332, 144
# runs): same-crown pairs give +1.6+-4.3 kT (nothing), span 225 nm +13.5+-3.1,
# span 450 nm +31.0+-10.2, i.e. ~0.065 kT/nm of span, and extra heads at FIXED
# span add nothing.  But it CANNOT BOOTSTRAP: every tilt=60 arm failed, and the
# reason is geometric.  A head reaches 80 nm transversely, so the axial extent
# of actin within reach of one thick filament is W(theta) = 2R/tan(theta).  A
# couple needs two crowns inside W, which caps the angle at which any couple
# exists:
#     hph=9  (225.0 nm spacing) -> couple only up to theta = 35.4 deg
#     hph=27 ( 75.0 nm spacing) -> couple up to theta = 64.8 deg
# An isotropic start is 57.3 deg mean / 60 deg median, so at hph=9 the couple is
# switched OFF at the angles the system starts from.  hph=27 is the smallest
# budget whose couple reaches past an isotropic start.  THIS RUN TESTS THAT.
#
# WHAT IS NEW.  Only am_heads_per_half.  The hph=9 arms already exist at 12
# seeds from this same binary and runner as span3_tilt60 / span3_tilt30 in
# confined_probe_v2_20260910, and seeds are 1..12 here to match them
# seed-for-seed.  After the run, diff all /meta attrs against those files and
# require am_heads_per_half to be the ONLY difference.
#
# NO dt CONFOUND.  am_max_heads_per_actin_half = 6 in the runner, so hph changes
# WHERE heads can attach (crown spacing) and not HOW MANY bind.  Total AM force
# is bounded identically in both arms, so dt and the displacement cap stay put.
# Record cumulative_force_cap_events anyway and check it did not move.
#
# SCORING -- these arms are BISTABLE.  Score P(latch) over the 12 seeds, never
# pooled means: the historical "22 deg" was an average over latched
# (theta ~ 2-6 deg) and escaped (theta ~ 13-68 deg) runs.  Read the couple from
# myosin/am_align_drive_{minus,plus} / am_align_head_count_{minus,plus}, which
# excludes the am_alignment_k tether and is a force-field property at a
# configuration.  ALWAYS group by (actin, myosin, HALF) -- pooling the two
# halves admits spans above the 675 nm bearing length, which are opposite-half
# pairs and invalid.
#
# MINIMAL-RIG CAVEAT: Lx=4, 1 actin + 1 myosin.  This is a two-body mechanism
# probe.  It is INADMISSIBLE for any lattice-period, registration or
# assembly-outcome claim, and it cannot predict production per-pair occupancy
# (production has ~6 actin partners per myosin, this has 1).
# ===========================================================================

src=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910
root=/scratch/users/xsli/sarcomere/probe_hph_20261005
exe=${src}/bin/sarcomere_conf2_20260910
runner=${src}/run_conf2.sh
mkdir -p "${root}"/{data,params,logs}
test -x "${exe}"; test -x "${runner}"

# arm | hph | transverse_k | align_k | tilt     (transverse_k/align_k match span3)
configs=(
  "hph27_tilt60|27|100|0|60"  "hph27_tilt60|27|100|0|60"  "hph27_tilt60|27|100|0|60"
  "hph27_tilt30|27|100|0|30"  "hph27_tilt30|27|100|0|30"  "hph27_tilt30|27|100|0|30"
)
# 3 array tasks per arm, 4 seeds each, so a task finishes well inside 8 h even
# though hph=27 iterates 3x the head slots.
seed_chunks=( "1 2 3 4" "5 6 7 8" "9 10 11 12" "1 2 3 4" "5 6 7 8" "9 10 11 12" )

IFS='|' read -r label hph transverse_k align_k tilt <<< "${configs[$SLURM_ARRAY_TASK_ID]}"
seeds=${seed_chunks[$SLURM_ARRAY_TASK_ID]}
nsteps=${PROBE_NSTEPS:-200000}        # 2 s at dt=1e-5, same as confined probe v2
CONFINE_K=50

echo "task=$SLURM_ARRAY_TASK_ID arm=${label} hph=${hph} tilt=${tilt} seeds='${seeds}' host=$(hostname)"

for seed in ${seeds}; do
    tag=${label}_seed${seed}
    out=${root}/data/${tag}.h5
    if [ -s "${out}" ]; then echo "skip ${tag} (exists)"; continue; fi
    printf '%s\n' "label=${label}" "am_heads_per_half=${hph}" \
        "am_head_transverse_k=${transverse_k}" "am_alignment_k=${align_k}" \
        "probe_tilt_deg=${tilt}" "probe_confine_k=${CONFINE_K}" \
        "am_tether_occupancy_floor=0" "box=4x4x4" \
        "myosin=frozen_via_crowding_1e9" "seed=${seed}" \
        "nsteps=${nsteps}" "exe=${exe}" "runner=${runner}" \
        "control_arm=confined_probe_v2_20260910/data/span3_tilt${tilt}_seed${seed}.h5" \
        >"${root}/params/${tag}.params"
    "${runner}" "${out}" \
      --Lx=4 --Ly=4 --Lz=4 \
      --am_heads_per_half=${hph} --am_head_transverse_k=${transverse_k} \
      --am_alignment_k=${align_k} --probe_tilt_deg=${tilt} \
      --probe_confine_k=${CONFINE_K} \
      --nsteps=${nsteps} --seed=${seed} \
      >"${root}/logs/${tag}.log" 2>&1
    echo "done ${tag}"
done
echo "task $SLURM_ARRAY_TASK_ID complete"
