#!/bin/bash
set -euo pipefail

binary=/scratch/users/xsli/sarcomere/actn2_local_capture_pilot_A800_M320_100ms_dt10us_20260819/bin/sarcomere_actn2_local_capture_20260819
test -x "${binary}" || { echo "missing executable: ${binary}" >&2; exit 1; }

# These trailing arguments intentionally override the historical Cell-A pilot
# launcher.  Everything not listed here remains byte-for-byte the earlier
# A800/M320, 10-us, module40 setup.
exec "${binary}" "$@" \
    --myosin_radius=0.05 \
    --myomesin_sites_per_filament=3 \
    --aa_arm_state_model=true \
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
    --aa_occlusion_rupture_penetration=0.010 \
    --aa_total_link_pool=360
