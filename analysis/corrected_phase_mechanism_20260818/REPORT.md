# Mechanism analysis of the corrected ACTN2 lifetime sweep

## Dataset

- Campaign: `/scratch/users/xsli/sarcomere/data_nalam_s2fstab_base1e-3`
- Launcher: `cpp/nalam_s2fstab_base1e-3.sbatch`
- Grid: 6 actin counts (600--2100), 13 directional lifetime coefficients, 3 replicas
- Fixed parameters: 1000 myosins, 12 x 2 x 2 box, 300,000 steps, 1-ms status-1/base lifetime
- Complete trajectories: 233/234. The excluded file is
  `A1800_M1000_s2fstab_lt1e-5_r01_L12x2x2.h5` (checkpoint 196,601).

This is the corrected legacy active-matter phase sweep. It is not the later
explicit-head/myosin-motor branch.

## Main result

There is no density-independent structural-order threshold when the
directional status-2 lifetime contribution exceeds the 1-ms base lifetime.
For loaded status-2 links, `f_load` is approximately 0.997, so the crossing is
indeed very close to `lambda/base = 1`. The chemistry responds sharply, but the
structural order does not.

The nearest sampled boundary contrast is
`Phi3D(lambda/base=3.16) - Phi3D(lambda/base=1)`. Across 18 paired
density/replica blocks its mean is -0.0151 (95% CI -0.0710 to 0.0408;
paired t-test p=0.576; Wilcoxon p=0.702). The broader within-block contrast,
mean above the boundary minus mean at or below it, is -0.0380 (95% CI -0.0985
to 0.0224; p=0.202).

The boundary effect changes sign with density:

| actins | immediate boundary delta Phi3D | mean above - at/below |
|---:|---:|---:|
| 600 | +0.0147 | +0.0707 |
| 900 | +0.0564 | +0.1086 |
| 1200 | -0.0792 | -0.0073 |
| 1500 | +0.0485 | -0.1103 |
| 1800 | -0.0468 | -0.1638 |
| 2100 | -0.0842 | -0.1263 |

Using `Phi3D >= 0.5` as an assembled-state label gives the same conclusion.
Above the crossing, the assembled fraction changes from 0 to 0 (A=600), 0 to
0.111 (A=900), 0.429 to 0.444 (A=1200), 0.952 to 0.389 (A=1500), 1.0 to
0.833 (A=1800), and 1.0 to 1.0 (A=2100).

## What actually happens at the crossing

The lifetime implementation adds `lambda * f_load` to the base lifetime of a
status-2 contact. Completed-bond lifetimes in a representative A=1800 replica
confirm the expected kinetic response:

| lambda/base | median lifetime (s) | 90th percentile (s) | completed bonds |
|---:|---:|---:|---:|
| 0 | 0.00068 | 0.00227 | 1,759,578 |
| 1 | 0.00134 | 0.00449 | 1,428,610 |
| 3.16 | 0.00270 | 0.00916 | 1,028,712 |
| 100 | 0.02795 | 0.13424 | 110,056 |

Across all densities, moving from lambda/base=1 to 3.16 raises mean status-2
occupancy from 0.585 to 0.763 and lowers unbinding from 236 to 171 events per
actin per second. Yet mean Phi3D falls slightly, from 0.554 to 0.538. At
lambda/base=100, occupancy is 0.975 and turnover is only 12.4 events per actin
per second, while mean Phi3D is down to 0.469.

## Assembly mechanism

1. The initial condition is already orientationally aligned. Tail myosin
   nematic order is 0.94--0.997 and changes little across the sweep. The phase
   transition therefore measures axial density-band formation, not generic
   nematic alignment.
2. Actomyosin engagement and status-2 conversion happen first. In assembled
   trajectories, the median half-time is 0.02 s for actomyosin engagement and
   0.06 s for status-2 occupancy, versus 0.33 s for axial myosin order.
3. Myosin-rich axial bands then nucleate and reorganize over roughly
   0.2--1 s. Reversible ACTN2 links allow bands to move, merge, and select a
   strong periodic mode.
4. Directional stabilization helps sparse systems because longer-lived links
   make otherwise weak nuclei persistent. This explains the modest positive
   response at A=600--900.
5. In dense systems, bands already form with lambda=0. Further stabilization
   fills nearly every available status-2 contact and suppresses turnover.
   Crosslinks then preserve competing or misregistered bands instead of letting
   them anneal. At A=1500--2100, status-2 fraction is negatively correlated
   with Phi3D (Spearman rho -0.46, -0.57, and -0.44, respectively).

The empirical picture is thus a density-dependent optimum, not a universal
`directional lifetime > base lifetime` switch. The observed crossing is a
kinetic threshold for bond persistence; whether it helps structural order is
set by the balance between nucleus survival and the turnover needed for
annealing.

## Outputs and reproducibility

- `run_summary.csv`: refreshed exact Phi3D and tail mechanism metrics for all
  233 complete runs
- `cell_summary.csv`: 78 parameter-cell summaries
- `threshold_tests.csv`: paired boundary and broad contrasts
- `mechanism_correlations.csv`: density-resolved rank correlations
- `temporal_lags.csv`: per-run half-times
- `phase_mechanism.png`: phase/order, status-2 occupancy, and turnover
- `representative_timeseries.png`: order and status-2 time series
- `myosin_density_kymographs.png` and `cb2_density_kymographs.png`: spatial
  assembly pathways for representative sparse and dense trajectories
- Analysis scripts: `analysis/analyze_corrected_phase_mechanism.py`,
  `analysis/summarize_corrected_phase_mechanism.py`, and
  `analysis/plot_corrected_phase_kymographs.py`

The current sweep has only three replicas per parameter cell, so the precise
density-dependent optimum is not resolved tightly. A decisive follow-up would
densely sample lambda/base from 0.5 to 10 at A=900, 1500, and 1800 with more
replicas, and measure coarsening/mode-switch rates rather than only final order.
