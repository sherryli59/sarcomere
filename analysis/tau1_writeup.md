# Is status-1 a physical intermediate, or was the old assembly a model artifact?

**Question.** In the original ("s1s2") kinetics, singly-engaged catch bonds
(`cb_status == 1`) were stabilized by `lifetime = λ·f_load` — a *load-dependent*
lifetime. That boosted 3D sarcomeric order (legacy myosin φ₃D) by up to +0.33
over the model that stabilizes only the doubly-engaged, load-bearing crosslink
(`cb_status == 2`, "s2only"). The concern: load-stabilizing a bond that is *not
yet load-bearing* is unphysical, so that +0.33 might be an artifact.

**Test.** Replace the load-dependent status-1 lifetime with a **constant,
load-independent** `τ1` (`--base_lifetime=τ1`, `--stabilize_status1=false`,
`--directional_base_lifetime=0`), catch stabilization on status 2 only:

```
lifetime(status 1) = τ1                    # intrinsic off-rate, no load term
lifetime(status 2) = τ1 + λ·f_load         # catch reserved for the crosslink
```

Sweep τ1 and ask where φ₃D lands between the **s2only floor** (τ1→0) and the
**s1s2 ceiling** (old load-dependent law). If a finite constant τ1 reaches the
ceiling, the boost was just "single bonds live a while" — physical and keepable.
If none does, the boost required the load-dependence, and s2only is correct.

## Result

**A constant, load-independent τ1 fully recovers the s1s2 boost.** At λ=1e-3
the recovery is complete and the curve is fully resolved:

| τ1 | τ1/λ | φ₃D | recovered fraction |
|---:|---:|---:|---:|
| 0 | 0 | 0.596 | 0% (= floor) |
| 3.16e-5 | 0.03 | 0.601 | ~0% |
| 1e-4 | 0.1 | 0.720 | 44% |
| 3.16e-4 | 0.32 | 0.863 | 96% |
| **1e-3** | **1.0** | **0.882** | **103%** |
| 3.16e-3 | 3.16 | 0.539 | **−21% (below floor)** |

Floor (s2only) = 0.597, ceiling (s1s2) = 0.874, gap = +0.277.

See `phi3d_tau1.png` (left panel). Three regimes:

1. **τ1 ≲ 0.03λ — dead.** A few KMC steps of lifetime; the bond breaks before it
   can register anything. φ₃D sits on the floor.
2. **τ1 ≈ λ — peak, reaches the ceiling.** 103% recovery. This is the
   "matched-to-old" point: since `f_load ≈ 1` exactly whenever a pair is loaded
   (measured median 1.0, max 1.0), the old status-1 lifetime *was* effectively
   the constant λ. A literal constant τ1=λ reproduces s1s2 to within noise.
3. **τ1 ≳ 3λ — overshoot, degrades below the floor.** Bonds that essentially
   never release freeze actin in whatever mis-registered contact it first made:
   kinetic trapping / glassy jamming.

λ=1e-2 is **not yet fully resolved** (right panel): the base sweep topped out at
τ1 = 0.32λ (recovered ≈ 0.48, still climbing), and the τ1 = 0.1λ cell has no
complete replica yet. The two rungs that bracket the peak — τ1 = λ and τ1 = 3.16λ
— were launched as `cpp/a2000_tau1_ext_lt1e-2.sbatch` (job 34333283) and will
confirm whether the ceiling and the overshoot reproduce at λ=1e-2.

## Is it physical?

**The *form* is more physical than the old law; the *value* still owes a
justification.**

- **A constant τ1 for the nascent single bond is plain first-order unbinding** —
  the Bell-model zero-force off-rate of a single molecular contact. Reserving the
  catch (force-stabilization) term for status 2 is the honest split: *a contact
  exists* (τ1) vs *a load-bearing crosslink has formed and now shows catch*
  (τ1 + λ·f_load). Catch behavior is documented for the load-bearing object, not
  for a nascent contact that isn't yet transmitting force.

- **The old law was unphysical two ways.** (1) It applied catch stabilization to
  a bond that isn't yet bearing load. (2) `f_load` is effectively *binary*
  (median 1.0), so `λ·f_load` was a *switch* ("loaded → live λ, unloaded → die
  instantly"), not the smooth biphasic force-lifetime curve of a real catch bond.
  All of the physics lived in the constant λ — which is exactly why a literal
  constant τ1≈λ reproduces it.

- **The overshoot is itself a physical constraint.** Assembly needs bonds to be
  *reversible on the assembly timescale* so mis-registered contacts can dissolve
  and re-anneal. That puts a hard *upper* bound on physical τ1, below the ~3λ
  jamming point the sweep exposes. τ1 must be nonzero **and** sub-jamming.

- **The one number to be suspicious of is the magnitude.** τ1 ≈ λ works, but if λ
  itself isn't independently justified as a timescale, "τ1 = λ recovers the old
  result" may just mean both parameters absorb the same fitted number. The
  defensible move: pin τ1 to a *measured* single-contact off-rate (actomyosin
  weak-binding lifetime, or the crosslinker zero-force k_off) through the dt→s
  mapping, then check it lands in the nonzero, sub-jamming window this sweep
  defines.

**Recommendation.** Keep status-1 with a constant τ1 and drop the load
dependence — physically cleaner, and it loses no assembly. Set τ1's magnitude
from a measured off-rate rather than from matching s1s2, staying inside the
[nonzero, ≲3λ] window.

## Methods / provenance

- φ₃D = legacy myosin-only definition (max over 3D Fourier modes), identical to
  `compare_phi3d_s1s2_vs_s2only.py` and `legacy_status2_directional_summary.csv`,
  so floor/ceiling numbers are directly comparable across campaigns.
- Each cell averaged over its own last 300 frames (steady state), completed
  replicas only (`checkpoint ≥ 290k steps`). Seeds = 97000+rep, shared across
  the tau1 / s1s2 / s2only sweeps, so cells at a given (λ, rep) are the *same*
  trajectory up to the first divergent kinetics event.
- dt = 1e-5, `break_probability = dt/lifetime` ⇒ mean status-1 lifetime in steps
  = τ1/1e-5. τ1 = 1e-5 would equal τ1 = 0 (p_break = 1) and is omitted.
- Data: `/scratch/users/xsli/sarcomere/data_a2000_tau1_sweep`.
  Analysis: `analysis/phi3d_tau1_answer.py` → `phi3d_tau1_answer.csv`.
  Plot: `analysis/plot_phi3d_tau1.py` → `phi3d_tau1.png`.

**Caveat.** Most cells are single completed replicas (n=1, ringed orange in the
plot); treat individual points as ±0.05. The λ=1e-3 *shape* (dead / peak /
overshoot) is robust; individual λ=1e-2 points are provisional until the
extension and second replicas land.
