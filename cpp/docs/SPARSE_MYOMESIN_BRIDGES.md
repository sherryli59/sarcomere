# Sparse local finite-valency myomesin modules

Date: 2026-08-18

## Why this replaces the old `k_mm`

The former `k_mm` term was not a molecular bond model. Every geometrically
near myosin pair received an always-on axial spring on every force evaluation.
It had no occupancy, no finite number of binding sites, no association or
dissociation, and no transverse rest length. Increasing `k_mm` could therefore
couple many neighbors simultaneously and help create a single overconnected
component.

With this change, `k_mm` is the axial stiffness of **one occupied local
module**. `k_mm=0` remains an exact opt-out. A positive value enables explicit
dynamic myomesin-module occupancy.

## Physical coarse graining

Each graph edge represents a local M-band myomesin module between two thick
filament bare zones. It is an effective coarse-grained link that may represent
several parallel myomesin-family molecules, not an assertion that one isolated
myomesin dimer spans the full lattice spacing. Consequently, `k_mm` and
`myomesin_force_max` are module-level parameters and must be scaled together
when changing the effective molecular occupancy represented by one edge.

The model has the following physical ingredients:

1. **Local capture.** A new bridge can form only if bare zones are within the
   axial capture window, close to the transverse lattice spacing, and their
   thick-filament axes are sufficiently parallel or antiparallel.
2. **Finite valency.** A bridge consumes one site on each thick filament.
   `myomesin_sites_per_filament` limits the number of simultaneous neighbors,
   and one myosin pair can carry at most one coarse-grained bridge.
3. **Reversible chemistry.** Local encounters associate with rate
   `myomesin_k_on`; occupied bridges dissociate with `myomesin_k_off`. An
   optional Bell transition distance accelerates dissociation under bridge
   load. Counter-based random streams make events independent of OpenMP
   scheduling.
4. **M-band register mechanics.** The axial rest offset is zero, so a bound
   bridge registers the two bare-zone centers. This is local pair mechanics,
   not attraction to a global M-band centroid.
5. **Finite transverse spacing.** The transverse spring has a nonzero rest
   distance, `myomesin_radial_rest`; it does not collapse the two thick
   filament centerlines onto each other.
6. **Reciprocal force.** The bridge applies equal and opposite forces to its
   two myosins and therefore cannot inject net linear momentum.
7. **Independent actomyosin chemistry.** Bridge association or release never
   deletes, releases, or changes the rates of any AM head. Ordinary individual
   head turnover continues unchanged while the M-band network anneals.

This is physically motivated by the observed localization of myomesin to the
M-band, its binding to myosin/titin, and the elastic response of its modular
chain. The topology and parameter values are still coarse-grained assumptions
that must be calibrated against M-band stoichiometry and mechanics.

## Parameters

| Option | Meaning | Default |
|---|---|---:|
| `k_mm` | effective axial stiffness of one occupied module (pN/um); zero disables | 0 |
| `myomesin_sites_per_filament` | maximum simultaneous bridges per myosin | 2 |
| `myomesin_k_on` | local association rate (1/s) | 5000 |
| `myomesin_k_off` | zero-load dissociation rate (1/s) | 1 |
| `myomesin_axial_capture` | new-bond axial capture window (um) | 0.25 |
| `myomesin_radial_rest` | transverse bridge rest spacing (um) | 0.11 |
| `myomesin_radial_capture` | allowed radial deviation at capture (um) | 0.04 |
| `myomesin_radial_k` | independently adjustable module transverse stiffness (pN/um); zero uses `k_mm` | 0 |
| `myomesin_alignment_cos` | minimum absolute axis dot product at capture | 0.94 |
| `myomesin_bell_distance` | load-sensitive off-rate distance (um) | 0.0005 |
| `myomesin_force_max` | module resultant-force cap (pN); zero disables | 10 |

Association gates apply only to formation. A bound bridge is not destroyed
merely because thermal motion carries it outside its original capture shell;
it persists until its own stochastic dissociation event.

## HDF5 diagnostics and resume

Metadata records the complete parameter set and rejects resume with a changed
bridge model. `/myomesin/bridge_state` is sparse and frame-ragged; each row is

```text
id_hi, id_lo, myosin_i, myosin_j, attach_step,
axial_offset, radial_distance, force
```

`/myomesin/frame_offsets` maps saved steps to their rows. `/myomesin/events`
records association and release events. `/resume` stores the exact occupied
graph, next molecule ID, diagnostics, and pending events. Arbitrary saved-frame
resume restores the sparse graph from the frame-ragged trajectory state.

## Deterministic tests

`myomesin_mband_test` verifies:

- restoring axial force and exact action-reaction symmetry;
- independent module stiffness and resultant-force capacity, including exact
  10 and 20 pN force clamping when `k_mm` and the cap are scaled together;
- zero transverse force at the nonzero rest spacing, attraction above it, and
  repulsion below it;
- rejection outside axial and radial capture windows;
- finite valency in a three-myosin local encounter;
- stochastic bridge turnover without any AM intervention;
- exact recovery of the `k_mm=0` behavior.

## Four-actin/four-myosin annealing challenge

The validation geometry is
`minimal_catch_balanced_double_myosin_anneal`: every actin is shared by two
myosins, and the two intended transverse M-band pairs start 150 nm out of axial
register. It retains the existing physical rod mobility, dynamic AM heads,
`motor_transmitted_register` ACTN2 state machine, and `k_aa=3000 pN/um`.

The matched no-myomesin control ended at 198.1 and 191.3 nm axial mismatch for
M0--M2 and M1--M3 (final-10-ms medians 195.6 and 190.6 nm), so the motor network
made the imposed 150-nm mismatch worse.

The acceptance metric is deliberately referenced to the exact association
event, not to the prescribed initial geometry or the first saved frame. For
each intended M--M pair, `analysis/check_myomesin_annealing.py` compares the
absolute axial mismatch in the association event with the median mismatch over
the final 10 ms. Both pairs must improve by at least 10%. This corrected an
earlier false-positive interpretation of the 60/10-pN run.

The module-level sweep keeps the transverse stiffness fixed at 300 pN/um while
scaling axial stiffness and force capacity together:

| axial / transverse stiffness (pN/um) | module force cap (pN) | M0--M2 formation -> final-10-ms median | M1--M3 formation -> final-10-ms median | final radial medians (nm) | result | wall time |
|---:|---:|---:|---:|---:|---:|---:|
| 60 / 300 | 10 | 145.16 -> 142.44 nm (1.9%) | 145.81 -> 139.89 nm (4.1%) | 93.27, 93.60 | fail | 43.57 s |
| 120 / 300 | 20 | 145.16 -> 106.06 nm (26.9%) | 145.82 -> 85.61 nm (41.3%) | 94.11, 94.28 | pass | 42.91 s |
| 240 / 300 | 40 | 145.16 -> 53.19 nm (63.4%) | 145.83 -> 53.60 nm (63.2%) | 94.08, 93.46 | pass | 42.49 s |

Every run formed exactly the two intended local modules, M0--M2 and M1--M3,
within 0.21 ms. The two new runs had no release or topology-change event. The
20-pN run reached 18.17 pN and the 40-pN run reached 30.26 pN, so neither new
run sat on its module cap. No global actin/myosin force, torque, or rotation cap
was invoked. Both ACTN2 links remained status 2; their final connector lengths
were 49--53 nm, close to the 50.6-nm rest length.

The 40-pN module gave the strongest and most symmetric registration in this
50-ms test. It is a calibration candidate, not yet a universal default: the
module stiffness and capacity must ultimately be related to the effective
number of myomesin-family molecules represented by one coarse edge.

Validation outputs:

```text
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/sparse-myomesin/k30_50ms.h5
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/sparse-myomesin/k60_rad300_50ms.h5
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/sparse-myomesin/module_k120_cap20_rad300_50ms.h5
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/sparse-myomesin/module_k240_cap40_rad300_50ms.h5
```

The complete CTest suite passed after the module calibration change: 28/28
tests. Exact latest-snapshot resume
restored two bridges with pairs `(0,2)` and `(1,3)` and next ID 3. An
arbitrary-frame resume at step 25,000 restored the same two-pair graph and
wrote a consistent new resume snapshot.

## Remaining calibration questions

- `myomesin_sites_per_filament` is an effective network valency and should be
  matched to the represented lattice resolution.
- Axial/radial stiffness and the module force capacity need calibration to the
  number of molecules represented by one graph edge.
- Constant zero-load turnover plus a Bell factor is the minimal reversible
  kinetics. It does not yet represent domain-by-domain unfolding/refolding.
- The bridge acts at the coarse-grained bare-zone center and therefore applies
  no separate orientational torque. Capture still requires local axis
  compatibility.

## Physical references

- Obermann et al., myomesin localization and M-band architecture:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC2121001/>
- Obermann et al., myomesin interactions with myosin and titin:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC1169628/>
- Pinotsis et al., elastic myomesin structure and mechanics:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3279516/>
- Lange et al., molecular structure of the M-band network:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC2814874/>
