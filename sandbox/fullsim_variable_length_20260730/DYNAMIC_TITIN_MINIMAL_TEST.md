# Dynamic tension-only titin: minimal-cell sandbox test

Date: 2026-07-31

This sandbox prototype was tested only in the periodic four-actin/two-myosin
balanced cell. It has not been copied into the production assembly driver.

## Implemented mechanism

- Each myosin half discovers its actin partner from the actual weak+strong head
  occupancy; no actin-to-myosin ownership table is used for capture.
- Capture also requires an active local ACTN2 link in the actin barbed-end
  zone, the correct myosin half, a maximum endpoint-to-myosin distance, and a
  one-millisecond occupancy dwell in the diagnostic configuration.
- There is at most one captured tether per myosin half and one tether owner per
  actin.
- The captured identity persists across fast head and ACTN2 turnover.
- An optional counter-RNG slip-bond off-rate permits slow release followed by
  local recapture. It is disabled in the first diagnostic so that the force
  mechanism can be tested without another timescale.
- The force law is strictly tension-only:

      extension = max(0, distance - rest_length)
      force = stiffness * extension

  Here `distance` is the absolute axial separation from the shared Z node
  (periodic centroid of the ACTN2-paired barbed ends) to the appropriate thick
  filament tip. The rest length is 0.375 um, obtained from the former 1.125-um
  Z-to-M-center target minus half the 1.5-um thick-filament length.
- The Z-node force is strictly axial, divided equally between its two actins,
  and applied reciprocally to the myosin. Titin applies no independent
  transverse force or director torque to either actin.
- Compression produces exactly zero force.

## Tests

The focused test covers topology discovery, dwell, ACTN2 support, identity
memory, one-link ownership, Z-partner identity, zero compression force, equal
pair-force allocation, reciprocal force balance, zero titin torque,
deterministic counter-RNG release, and recapture readiness.

All 23 sandbox tests passed before the final optional-turnover addition. After
that addition, the affected `counter_rng` and `axial_cell` tests passed again,
and a full-Brownian diagnostic smoke test captured all four intended links.

## 500 ms Brownian comparison (seed 1001)

Both runs use 35 heads per myosin half, 50 ms fixed-geometry chemistry
pre-equilibration, `Lx = 4.7 um`, initial `AA = 0.16 um`, and the
`aa_k100x_register_all_l3` ACTN2 mechanics.

| condition | final A-A overlap | final A-M overlap | final tether state | final total tether force |
|---|---:|---:|---|---:|
| axial Z-node--thick-tip tethers | 0.0995 um | 0.5255 um | 4 bound, 4 tensile, 0 compressed | 41.7 pN |
| no tether | 0.0953 um | 0.6307 um | none | 0 pN |

The corrected tethered trajectory was extended to 1 s. Over its final 500 ms,
A-A overlap averaged 0.0977 um and A-M overlap averaged 0.5262 um. Their fitted
slopes were -0.0047 um/s and +0.0059 um/s, respectively. At 1 s, all four
tethers remained tensile with a mean axial tip gap of 0.4218 um and total force
of 41.2 pN. No force, torque, or rotation cap fired. By contrast, the no-tether
control advanced at approximately +0.878 um/s over 0.4--0.5 s.

## Interpretation and unresolved point

The corrected minimal test supports a slide-to-stall effect without the former
3D endpoint-to-center brace. It simultaneously keeps A-A near its 0.10-um
target and holds A-M near 0.526 um for 1 s. This is still a one-seed mechanism
test with titin release disabled; rest-length/stiffness sensitivity, multiple
seeds, and an enabled slow release/recapture control remain necessary before
promotion to the large assembly.

Primary outputs:

- `diagnostics/dynamic_titin/znode_tip_full_500ms.csv`
- `diagnostics/dynamic_titin/znode_tip_full_500ms_titin.csv`
- `diagnostics/dynamic_titin/znode_tip_full_1s.csv`
- `diagnostics/dynamic_titin/znode_tip_full_1s_titin.csv`
- `diagnostics/dynamic_titin/none_full_500ms.csv`
