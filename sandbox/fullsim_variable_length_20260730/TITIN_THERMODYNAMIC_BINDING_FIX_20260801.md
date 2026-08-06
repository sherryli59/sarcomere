# Thermodynamic explicit-titin binding fix (2026-08-01)

Scope: sandbox `explicit_titin_hex_fixed` diagnostic only. No production run was launched.

## Changes

- Explicit-titin elastic energy is now
  `U = 0.5 k max(0, d - L_slack)^2`.
- Optional passive strain gating uses
  `k_on = k_on0 exp(-beta theta U)` and
  `k_off = k_off0 exp(beta (1-theta) U)`, with exponent limits for numerical safety.
- The sandbox hexagonal diagnostic enables strain gating with `theta = 0.5`.
- The effective coarse-grained titin parameters changed from
  `k = 220 pN/um, L_slack = 0.375 um` to
  `k = 40 pN/um, L_slack = 0.500 um`.
- Compression remains exactly force-free.
- Titin CSV output now records candidate distance, elastic energy, effective on-rate,
  and the distance, energy, and force at the most recent attachment.

## Seed-1001 regression

| Metric | Previous | Thermodynamic titin |
|---|---:|---:|
| A4 final offset from A3/A5 | 456.7 nm | -0.9 nm |
| Maximum central A4 excursion | 456.7 nm | 26.0 nm |
| A4 distinct ACTN2 neighbors | 0 | 2 |
| Distinct ACTN2 pairs | 10 | 12 |
| Maximum titin attachment force | about 28 pN | 2.00 pN |
| Final linked A--A overlap | 90.7 nm | 95.5 nm |
| Final A--M overlap | 566.6 nm | 528.2 nm |
| Final Z-radius CV | 0.46 | 0.28 |

All 24 sandbox tests pass. No force, torque, or rotation cap was reached.

## Four-seed robustness check

All four seeds avoided the previous hundreds-of-nanometers single-filament escape.
Seeds 1001--1003 closed all 12 ACTN2 neighbor pairs. Seed 1004 formed a disconnected
four-actin cycle plus a two-actin dimer at the central barbed-end region, leaving one
actin degree short and two axially registered subassemblies separated by about 156 nm.
That residual is an ACTN2 network-topology/annealing limitation, not a titin attachment
impulse. It should be addressed separately without adding a centroid or pre-existing
Z-disc constraint.
