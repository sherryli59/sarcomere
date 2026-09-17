#!/bin/bash
set -euo pipefail

binary=/scratch/users/xsli/sarcomere/actn2_two_arm_pilot_A800_M320_100ms_dt10us_20260819/bin/sarcomere_actn2_two_arm_20260819
arm_model=${ACTN2_ARM_MODEL:?ACTN2_ARM_MODEL must be false or true}
case "${arm_model}" in
    false|true) ;;
    *) echo "ACTN2_ARM_MODEL must be false or true" >&2; exit 2 ;;
esac

echo "actn2_two_arm_wrapper_model=${arm_model}" >&2
exec "${binary}" "$@" \
    --aa_arm_state_model="${arm_model}" \
    --aa_arm_open_rate_weak=20 \
    --aa_arm_open_rate_mature=1 \
    --aa_arm_bell_distance=0.001 \
    --aa_tether_off_rate=5 \
    --aa_tether_rebind_rate=100 \
    --aa_tether_capture_distance=0.10 \
    --aa_arm_rebind_theta=0.5 \
    --aa_actn2_contact_radius=0.005 \
    --aa_myosin_contact_epsilon=0.00428 \
    --aa_myosin_contact_force_cap=40 \
    --aa_occlusion_rupture_penetration=0.010
