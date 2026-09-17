#!/bin/bash
#SBATCH --job-name=confv2
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --array=0-11
#SBATCH --output=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1     # 2 filaments; 1 thread also makes runs reproducible

root=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910
exe=${root}/bin/sarcomere_conf2_20260910
runner=${root}/run_conf2.sh
mkdir -p "${root}"/{data,params,logs}
test -x "${exe}"; test -x "${runner}"
case "$(basename "${exe}")" in *_20260910) ;; *) echo "exe/campaign date mismatch"; exit 3 ;; esac

# CONFINED ALIGNMENT PROBE v2 -- fixes both flaws of 42753933.
#
# FIX 1: box 4 x 4 x 4 instead of 4 x 1 x 1.  Half-box is now 2.0 um, far beyond
#   any observed excursion.  In v1 the half-box was 0.5 um in y/z and the actin
#   routinely sat at 0.95-1.00 of it, where two periodic images are equidistant,
#   so the minimum-image choice flipped and the drawn actin teleported by exactly
#   one box vector.  Worse, `am_tether_occupancy_floor_reach` was 1.0 um -- twice
#   the half-box -- so the confinement could pull the actin toward a periodic
#   GHOST of the myosin.  Every interaction reach must stay under half the box.
#
# FIX 2: a real centre trap (`probe_confine_k`) replaces the tether-floor hack,
#   and the floor is switched off (`am_tether_occupancy_floor=0`).  The pair
#   radial tether constrains only the segment-segment CLOSEST APPROACH, which a
#   crossed rod can satisfy with its centre 500 nm away -- so in v1 the low-head
#   control arms drifted off and read as false nulls.  The trap is a harmonic
#   spring on the actin's centre toward its initial position, applied AT the
#   centre so it contributes exactly zero torque and leaves orientation free.
#
# Tuned at k=50 pN/um: perp median 84 nm (target 80), perp max 146 nm, 88.8% of
# frames bound, 66.6% with >=2 heads, 62.1% with span>0, max drawn jump 56 nm
# (v1: up to 3942 nm), force-cap events 4.4% of steps.
#
# The myosin is frozen via crowding=1e9 (crowding divides the mobility and the
# integrator uses sqrt(2 D dt) noise, so drift and noise vanish together).
#
# ARMS -- now genuine controls, because escape is prevented for all of them:
#   span3   hph=9 -> 3 crowns, spans 0/225/450 nm.  Expect alignment.
#   span1   hph=3 -> crown_count=(3+2)/3=1, span identically 0.  If span is the
#           operative variable, expect NO alignment despite being engaged.
#           Note hph=3 also cuts the head budget, so condition on bound count.
#   notrans hph=9, am_head_transverse_k=0.  Tests whether the couple is the
#           radial head spring.
#   alignON hph=9, am_alignment_k=120.  Positive control.
#
# 12 seeds per (arm, tilt): these arms are BISTABLE and must be scored on
# P(latch), not on pooled means -- v1's "22 deg" was an average over a mixture
# of latched (theta ~ 2-6 deg) and escaped (theta ~ 13-68 deg) runs.

configs=(
  "span3_tilt10|9|100|0|10"     "span3_tilt30|9|100|0|30"     "span3_tilt60|9|100|0|60"
  "span1_tilt10|3|100|0|10"     "span1_tilt30|3|100|0|30"     "span1_tilt60|3|100|0|60"
  "notrans_tilt10|9|0|0|10"     "notrans_tilt30|9|0|0|30"     "notrans_tilt60|9|0|0|60"
  "alignON_tilt10|9|100|120|10" "alignON_tilt30|9|100|120|30" "alignON_tilt60|9|100|120|60"
)
IFS='|' read -r label hph transverse_k align_k tilt <<< "${configs[$SLURM_ARRAY_TASK_ID]}"

nsteps=${CONF2_NSTEPS:-200000}        # 2 s at dt=1e-5
CONFINE_K=50

for seed in $(seq 1 12); do
    tag=${label}_seed${seed}
    printf '%s\n' "label=${label}" "am_heads_per_half=${hph}" \
        "am_head_transverse_k=${transverse_k}" "am_alignment_k=${align_k}" \
        "probe_tilt_deg=${tilt}" "probe_confine_k=${CONFINE_K}" \
        "am_tether_occupancy_floor=0" "box=4x4x4" \
        "myosin=frozen_via_crowding_1e9" "seed=${seed}" \
        >"${root}/params/${tag}.params"
    "${runner}" "${root}/data/${tag}.h5" \
      --Lx=4 --Ly=4 --Lz=4 \
      --am_heads_per_half=${hph} --am_head_transverse_k=${transverse_k} \
      --am_alignment_k=${align_k} --probe_tilt_deg=${tilt} \
      --probe_confine_k=${CONFINE_K} \
      --nsteps=${nsteps} --seed=${seed} \
      >"${root}/logs/${tag}.log" 2>&1
    echo "done ${tag}"
done
