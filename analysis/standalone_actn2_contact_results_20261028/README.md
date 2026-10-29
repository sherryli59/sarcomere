# Standalone ACTN2 contact results

[Interactive trajectory viewer](standalone_contact.html) · [Full experiment report](report.md) · [Multi-seed refinement report](refinement/convergence.md)

Completed with one CPU in Slurm allocation 50081778. This standalone model uses illustrative parameters and an explicit external scaffold. Production simulator files were not edited by this work.

- At 20 pN drive, alpha=1 ends 252.6 nm past the original scaffold plane; alpha=3 ends 42.4 nm before it. At 150 pN the myosin passes despite alpha=3.
- In the 64-realization viewer ensemble at 40 pN, 79.7% lose at least one connector and 75% pass the plane. Of the 48 passing realizations, 47 retain live connectors. The lateral-escape figure shows how asymmetric loss opens a sideways route.
- The independent three-seed audit uses 768 failure realizations at each resolution: any-link loss is 83.7% at 100 ns and 85.9% at 50 ns. The difference is 1.2 standard errors. The audit passes its stated 6-standard-error comparison, which detects large resolution effects but is not an equivalence test.
- The same audit uses 12,288 weak-site realizations per condition/resolution. Inner-edge survival is 30.6% at 10 us and 30.9% at 2.5 us; barbed-end survival is zero in both cases. The main experiment separately validates both physical endpoints against analytic first-passage survival.
- Mechanical work-balance residuals fall from 8.39 to 4.19 to 2.09 pN nm under timestep halving. The final experiment also compares complete deterministic center trajectories: maximum change is 0.389 nm under combined timestep/quadrature refinement and 0.365 nm under quadrature refinement alone.

## Figures and data

- [Formation, association, and force maps](formation_repulsion.png)
- [Mechanics and failure](mechanics_failure.png)
- [Lateral escape](lateral_escape.png)
- [Weak-site boundaries](weak_boundaries.png)
- `trajectories.npz` saves mechanical trajectories, forces, live masks, local death codes, weak-site histories and final populations, and statistical clamp curves.
- `summary.json` records checks, parameters, seeds, versions and per-case summaries. `refinement/convergence.json` records the larger audit.
- `manifest.json` records file hashes and verifies that the final numerical kernels, model parameters, weak-site mechanics and mechanical integration match the audited source snapshot. The audit snapshot is retained in `refinement/`.

The HTML is self-contained and includes all four figures plus six rotatable, playable mechanical cases. Every one of the 64 failure realizations is selectable. Its JavaScript/data wiring passed 207 frame checks across 69 case/realization combinations. Actual browser rendering was not tested; open the saved HTML in a modern browser.

## Reproduce

The standalone script needs NumPy, SciPy and Matplotlib. The scheduler launcher aggregates the primary experiment and larger audit; read current Sherlock instructions before submitting. From the repository root, use a new result directory:

```bash
sbatch --output="$SCRATCH/actn2-standalone-%j.log" analysis/run_standalone_actn2_contact.sbatch "$SCRATCH/actn2-new-results"
```

Copy results that must be retained out of purged scratch storage. This result directory is already preserved in the project under GROUP_HOME.
