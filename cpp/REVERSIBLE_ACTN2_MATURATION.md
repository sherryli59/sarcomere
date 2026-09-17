# Reversible ACTN2 maturation

## Purpose

Status 2 represents a mechanically stabilized ACTN2 conformation, not a
permanent historical label. A mature bridge now may relax to status 1 while
remaining attached to both actins. This lets the weak bridge resume local
annealing after sustained loss of productive motor support.

The default rate is zero, so old inputs and resumed trajectories keep the
previous irreversible behavior unless the feature is explicitly enabled.

## Kinetics

The mature-to-weak rate is

```
k_unlock = k0 exp[-F_filtered/F_scale
                  + max(0, delta-delta_off)/delta_scale].
```

- `F_filtered` is the existing motor-timescale support signal. In
  `motor_transmitted_register` it is the source-resolved support carried by the
  locked motor pair. In `productive_register` it is the aggregate productive
  motor-support EMA. Passive A-A springs and steric forces do not directly
  create either signal.
- `delta` is the barbed-end registry error and `delta_off` is the existing
  mature-state registry retention width.
- A transition changes status 2 to status 1 but does not detach ACTN2 or move
  either material attachment coordinate. Subsequent weak-state diffusion or
  partner exchange uses the already selected annealing model.
- In `motor_transmitted_register`, motor identities remain fixed while mature.
  Unlocking releases that identity lock; it does not unbind any myosin head.

The per-step probability is `1-exp(-k_unlock*dt)` and uses a dedicated
counter-based random stream keyed by link ID and timestep.

## Runtime controls

| option | units | default | suggested first pilot |
|---|---:|---:|---:|
| `aa_status2_unlock_rate` | 1/s | 0 | 20 |
| `aa_status2_unlock_force_scale` | pN | 0.1 | 0.1 |
| `aa_status2_unlock_registry_scale` | um | 0.02 | 0.02 |

`aa_status2_locked_until_detach=true` is incompatible with a positive unlock
rate and is rejected by the command-line driver.

## Contact/occlusion distinction

The pre-existing ACTN2-myosin corona is a conservative segment-to-segment
excluded-volume force on the ACTN2 connector. It is not an axial A-M overlap
end-stop. With `myosin_radius=0.05 um` and
`aa_actn2_contact_radius=0.005 um`, the hard-core distance is 0.055 um. The
previous 0.015-um shell therefore ends at 0.070 um, short of the usual
0.080-um A-M centerline. A 0.035-um test shell reaches 0.090 um and acts before
physical contact, but its force remains normal to the nearest connector/myosin
segments. It therefore can reduce bridge intrusion without guaranteeing that
axial actin sliding stops before the scalar A-M occlusion boundary.

The 85,329 `am_overlap_occlusion` records in the completed trajectory are
cumulative detachment events, not 85,329 simultaneously present bridges. They
show that this geometric failure mode recurred often enough to matter; the
corona range and the reversible maturation kinetics should be tested as
separate interventions.

## Minimal paired pilot

Keep the completed Cell-A-plus-myomesin setup and seed fixed. Compare:

1. reversible maturation with the existing 0.015-um corona; and
2. the same reversible kinetics with a 0.035-um corona.

Both arms use:

```
--aa_status2_locked_until_detach=false
--aa_status2_unlock_rate=20
--aa_status2_unlock_force_scale=0.1
--aa_status2_unlock_registry_scale=0.02
```

For the range arm, use

```
--aa_myosin_contact_shell_width=0.035
--aa_myosin_contact_shell_energy=0.1165
```

The larger energy scales as width squared from 0.0214 pN um, retaining the
same small-penetration shell stiffness rather than weakening the shoulder when
its range is enlarged. Score mature-to-weak transitions, mature occupancy,
AM-occlusion detachments per link-time, zero-crosslinkable-contour occupancy,
force-cap fractions, and Z/M order. This paired design separates the kinetic
change from the contact-range change.

## Diagnostics

- Status-change event code 3 records each 2-to-1 transition.
- Motor-pair event code 6 records release of a mature identity lock to the weak
  state.
- HDF5 metadata records all three unlock parameters.
- Occlusion detachment remains explicitly labeled
  `am_overlap_occlusion`; deep connector/thick-filament collision remains
  `actn2_myosin_collision`.
