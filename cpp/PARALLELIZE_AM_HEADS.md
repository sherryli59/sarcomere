# Parallelization proposal: `_process_dynamic_am_heads()`

## Why it's currently serial
The dominant per-step cost runs single-threaded because the call site wraps it in
`#pragma omp single`:

```cpp
// src/sarcomere.cpp ~415
#pragma omp single
{
  if (dynamic_am) {
      _update_thin_filament_activation();
      _process_dynamic_am_heads();      // 45k head slots, one thread
  }
  _enforce_myosin_bond_limit();
}
```

Measured effect: ~2–3 effective cores out of 16–32 (Amdahl — this is ~90% of the step).

## Why it is safe to parallelize over myosins
The function's hot loop is `for (int j = 0; j < myosin.n; ++j)`, and each myosin's
head KMC is independent. The code is already thread-ready:

| concern | status | evidence |
|---|---|---|
| RNG | already per-thread | `rng_engines[thread_id]`, sized `omp_get_max_threads()`, seeded `seed+t` |
| transition log | already per-thread | `am_transition_events_temp[thread_id]` in `_record_am_transition` |
| scratch connections | already per-thread | `actinIndicesPerMyosin_temp[thread_id]` |
| capacity accounting | per-myosin-iteration local | `occupied_by_actin`/`strong_by_actin` declared inside the `for j` body |
| filament forces | not written here | forces applied later in the existing parallel `omp for` loops |
| head state writes | per-myosin | `am_heads_by_myosin[j]` indexed by the loop var |

Only one shared structure needs care: the `candidates` list built in phase 1.

## Proposed changes (3 edits)

### Edit 1 — call site (~line 414-423): stop serializing the head work
```cpp
if (dynamic_am) {
    #pragma omp single
    { _update_thin_filament_activation(); }     // cheap; keep serial
    _process_dynamic_am_heads();                 // now contains its own omp for
}
#pragma omp single
{ _enforce_myosin_bond_limit(); }
#pragma omp barrier
```
(The function is now entered by ALL threads in the existing `#pragma omp parallel`
region, so its internal work-sharing directives take effect.)

### Edit 2 — `candidates` must be shared, not a per-thread local
Currently `std::vector<std::vector<Candidate>> candidates(myosin.n);` is a function
local — if all threads enter, each gets its own empty copy. Promote it to a shared
member (`am_candidates_`, resized once) OR declare before the region. Build it under
`omp single`, then barrier:

```cpp
// phase 1: candidate discovery (scatter into candidates[j]) — keep serial, it is cheap
#pragma omp single
{
    for (auto& c : am_candidates_) c.clear();
    for (int i = 0; i < actin.n; ++i) { /* ... existing gate ... */
        am_candidates_[j].push_back({...});
    }
}   // implicit barrier
```

### Edit 3 — parallelize the myosin KMC loop (the expensive phase)
```cpp
#pragma omp for schedule(dynamic)
for (int j = 0; j < myosin.n; ++j) {
    /* ... existing per-myosin head KMC, unchanged ...
       reads am_candidates_[j]; uses rng_engines[thread_id];
       occupied_by_actin/strong_by_actin stay local to this iteration ... */
}
// phase 3 crosslink_ratio copy:
#pragma omp for
for (int i = 0; i < actin.n; ++i)
    actin["crosslink_ratio"][i] = actin_crosslink_ratio[i];
```

## Expected result
Phase 2 is embarrassingly parallel over 160+ myosins → near-linear scaling until
load imbalance (heads-per-myosin varies). Realistic ~10–20× on 16–32 cores, which
combined with `am_heads_per_half=24` turns a "months" full-scale run into "days".

## Correctness caveats to weigh BEFORE applying
1. **Reproducibility changes.** Parallel KMC consumes `rng_engines[t]` in a
   thread-dependent order, so trajectories become a function of thread count and are
   no longer bit-identical to the serial version (they remain *statistically*
   equivalent). Resume stays exact only if thread count is held fixed. If the serial
   `omp single` was a deliberate determinism choice, this is the tradeoff.
2. **Validation gate before trusting it:** run the prebound thermal lattice serial vs
   parallel (N threads) and require matched distributions — tail |psi6|, weak/strong
   head counts, status-2 occupancy, AM transition-rate histograms — within noise.
   Then a dt/thread-count sweep to confirm no new cap events.
3. Keep an env/flag escape hatch (e.g. `OMP_NUM_THREADS=1`) to recover the exact
   serial trajectory for debugging.

## Suggested rollout
Branch `parallel-am-heads`; implement Edits 1–3; add a serial-vs-parallel statistical
regression test; benchmark 320/160 and 1280/640; only merge if the distributions match.
