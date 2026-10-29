# Standalone ACTN2 formation, repulsion, and failure

This is an independently runnable implementation of the current proposal, not a production simulator change or a calibrated molecular model.

## Findings

- Formation: full-connector maturation rejects a connector whose midpoint lies outside 130 nm; radial and axial conditions at different points do not cause a false rejection. Weak birth and maturation use distinct rules.
- With an explicit elastic scaffold and 20 pN drive, alpha=1 finishes at 252.6 nm penetration; alpha=3 finishes at -42.4 nm. Positive means the physical tip is past the original scaffold plane. This result is conditional on the scaffold.
- Alpha=3 does not eliminate sideways loading: at 20 nm burial the 40 nm outward connector carries approximately 19.9 pN axial and 30.3 pN radial. A point at its inner endpoint is almost purely axial.
- A 150 pN drive exceeds the four-pair axial capacity (at most 120 pN) and passes, ending at 1735.2 nm penetration. This is a bounded soft resistance, not an impenetrable wall.
- In 64 coupled 40 pN failure realizations, 79.7% lose at least one link and 1.6% lose all four by 2 ms. No connector is killed merely for entering the buffered region.
- 75.0% of the failure ensemble pass the original scaffold plane; 47 realizations do so with live connectors. Asymmetric loss permits lateral myosin escape, shown in lateral_escape.png.
- Weak sites detach at the barbed endpoint or cross the 200 nm structural edge and incur a higher hazard; the reflecting reference instead confines surviving sites. Failure and edge parameters are illustrative.

## Implemented rules

- Exact same-point full-connector maturation gate: radius 130 nm and physical half-length 750 nm; both weak sites must remain in the 200 nm barbed zone and the supplied registry condition must hold.
- Weak-birth physical-core veto using a 55 nm inflated myosin capsule; incremental maximum-pair contact-energy penalty exp(-positive delta U/kBT). This physical-core reference is an explicit standalone approximation, not a port of all production sterics.
- Smoothed weighted coordinate p=max(0,smin(alpha*(R-rho_epsilon),H-z_epsilon)); H=800 nm, epsilon=5 nm. Quadratic-then-linear Phi with k=1000 pN/um and saturation slope 30 pN. At alpha=3, radial/net force is bounded by 90 pN per connector–myosin interaction, axial force by 30 pN.
- Normalized integral over the entire connector, with endpoint derivatives weighted by (1-t) and t. Equal-and-opposite myosin reaction and exact physical torque from the same energy. Maximum-energy budget over redundant links is tested; distinct pairs add.
- Weak diffusion along fixed antiparallel actins, D=0.01 um^2/s and mobility D/kBT; actual contact-gradient drift, absorbing physical endpoints with local Brownian-bridge correction; soft inner-edge hazard 100/s * sum(([b-200 nm]+/20 nm)^2). Baseline weak off-rate is 10/s. Weak cohorts are held weak in these boundary experiments.
- Failure loading experiment: baseline 2/s plus contact hazard 50/s * [exp(T/5 pN)-1], where T is the larger endpoint contact-force magnitude. Contour rupture at 120 nm. Competing integrated-hazard clocks; no deterministic A-band deletion. These failure parameters are deliberately illustrative, not fitted.

## Numerical and stochastic checks

- Maximum endpoint finite-difference force error: 1.08e-08 pN; myosin translation error: 5.28e-09 pN.
- Maximum orientation-gradient error: 6.16e-08 pN um; total contact torque residual: 3.55e-14 pN um.
- Formation agrees with an independent 3,001-point segment reference on 1,200 random segments, plus explicit adversarial fixtures.
- Point-force bounds checked on 10,000 random points at both alpha values. Duplicate-link invariance, distinct-pair additivity and incremental redundant birth energy pass.
- Closed-loop work, including maximum-energy winner changes, at 128 / 512 / 2048 samples per side: [-0.12128757782914334, -0.007512424657907246, -6.806456731567323e-06] pN nm.
- Timestep/quadrature checks: {"myosin_center_dt_plus_quadrature_change_nm": 0.3885335651674926, "myosin_center_quadrature_only_change_nm": 0.36532989436680735, "free_relaxation_energy_monotone": true, "dt_s": 1e-07, "refined_dt_s": 5e-08, "fixed_Gauss_nodes": 16, "refined_nodes": 32, "trajectory_scope": "maximum center difference over all 201 common saved times", "scope": "deterministic supported 20 pN rig; does not certify stochastic production integration"}.
- Absorbing endpoint survival checked against the analytic first-passage law with 20,000 walkers. Weak-birth Bernoulli acceptance and fixed-contact-load lifetimes checked against their analytic probabilities (20,000 realizations per condition).

## Controlled mechanical rigs

Four distinct radial connectors have endpoints initially at 80 and 120 nm radius, on the same X plane, in four azimuthal directions. Myosin initially overshoots by 25 nm. Each endpoint translates with mobility 70 um/(pN s), myosin with 19 um/(pN s); connector extension stiffness 1000 pN/um and rest length 40 nm.

The supported cases attach each endpoint to an explicit harmonic external scaffold with stiffness 1000 pN/um. Those springs provide the missing surrounding-network support: a stopped myosin in this rig is not evidence that an unconstrained production Z-disc would hold. Body rotations and body Brownian noise are omitted. The free-relaxation case has no anchors or drive.

Fixed-topology cases isolate contact mechanics. The stochastic ensemble includes molecular death but no new births, remodeling or replacement. All realizations are included in statistics; the interactive view starts with realization zero and permits selecting every realization.

Penetration is measured relative to the original scaffold plane X=750 nm, not a moving connector. `connector_axial_overlap` separately records the current tip position minus each connector's nearest axial endpoint; radial overlap and the live mask must also be considered.

| Case | Drive (pN) | Alpha | Final scaffold penetration, median (nm) | Any link death | All links dead |
|---|---:|---:|---:|---:|---:|
| alpha1 | 20 | 1 | 252.55 | 0.000 | 0.000 |
| alpha3 | 20 | 3 | -42.38 | 0.000 | 0.000 |
| overdrive | 150 | 3 | 1735.15 | 0.000 | 0.000 |
| loaded_control | 40 | 3 | -34.86 | 0.000 | 0.000 |
| failure | 40 | 3 | 1490.85 | 0.797 | 0.016 |
| free | 0 | 3 | -23.83 | 0.000 | 0.000 |

## Weak-site boundary outcomes

| Initial condition | Rule | Alive at 10 ms | Endpoint loss | Soft-edge loss | Baseline loss |
|---|---|---:|---:|---:|---:|
| barbed | reflecting reference | 0.899 | 0 | 0 | 103 |
| barbed | absorbing + soft edge | 0.000 | 1019 | 0 | 5 |
| inner | reflecting reference | 0.899 | 0 | 0 | 103 |
| inner | absorbing + soft edge | 0.282 | 0 | 664 | 71 |

## Interpretation and limits

The unweighted geometry fails the intended near-tip direction test. Alpha=3 repairs point-level axial capture over a finite range, but an outward finite connector can still carry a larger radial than axial net force. Collective cancellation of lateral reactions on a myosin does not eliminate the sideways loading on each actin.

Energy-weighted weak association is effectively a hard veto through much of the 50 nm approach zone at these illustrative energies. The separate maturation gate can allow a location where weak-birth acceptance is extremely small. This is a substantive consequence to calibrate, not an implementation error to hide.

Failure tests verify the specified stochastic law and its interaction with this supported toy rig. They do not show that the chosen off-rate matches ACTN2, or that the production model will stop myosins, preserve connectivity, or share loads correctly.

This script is nonperiodic and uses one myosin. It does not reproduce production motor chemistry, catch-bond lifetimes, rigid-rod rotations, the full force-cap policy, multiple myosin images, site-registry kinetics, arm-state turnover, or network reformation. Those still require C++ kernel/integrator tests and the staged production comparisons in the implementation plan.

## Files and reproduction

- `standalone_actn2_contact.py`: copied standalone source, with no simulator imports.
- `summary.json`: all numerical checks, parameter values, seeds, timings, fixture outputs and experiment summaries.
- `trajectories.npz`: mechanical positions, contact endpoint forces, live/dead masks, local death codes and energies; weak-site histories/final cohorts; statistical clamp curves.
- `formation_repulsion.png`, `mechanics_failure.png`, `lateral_escape.png`, `weak_boundaries.png`: scientific figures.
- `standalone_contact.html`: self-contained rotatable 3D trajectory viewer and all figures; no WebGL, external dependency or internet connection needed.

From a compute allocation with numpy, scipy and matplotlib:

```bash
python standalone_actn2_contact.py --output /path/to/new_empty_results --replicates 64 --seed 20261028
```

The output directory must be empty; the script refuses to overwrite existing results. Read current Sherlock agent instructions before scheduling work. A one-CPU run aggregates the deterministic, statistical, stochastic-ensemble and rendering work.

Wall time: 51.69 s. Python 3.14.2, NumPy 2.4.1. Seed 20261028. No production source or trajectories were modified.
