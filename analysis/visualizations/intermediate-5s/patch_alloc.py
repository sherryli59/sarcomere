#!/usr/bin/env python3
"""Allocate the per-frame byte budget across the three layers.

At 100 frames inside 1.6 MB the budget is ~1000 segments/frame, while actin and
myosin alone average ~1150, so ACTN2 cannot be the only thing thinned. Split the
budget explicitly -- 40% actin, 25% myosin, the rest ACTN2 -- and derive a
keep-fraction per layer. Subsets are fixed in time (a filament index set, and a
hash of the actin pair for bridges) so filaments and bridges persist across
frames instead of flickering.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/build_100.py")
s = p.read_text()
if "keep_actin" in s:
    sys.exit("already patched")

s = s.replace("def geometry(h, f, box, la, lm, keep_pair):",
              "def geometry(h, f, box, la, lm, keep_pair, keep_actin=None, keep_myo=None):")
s = s.replace("""    k2 = st == 2""",
              """    k2 = st == 2
    if keep_actin is not None:
        k2 = k2 & keep_actin""")
s = s.replace("""    my = np.asarray(pipeline.clip_periodic_segments(cm-mh, cm+mh, box),
                    np.float32).reshape(-1, 6)""",
              """    km = slice(None) if keep_myo is None else keep_myo
    my = np.asarray(pipeline.clip_periodic_segments((cm-mh)[km], (cm+mh)[km], box),
                    np.float32).reshape(-1, 6)""")

old = """    budget = TARGET_KB/len(picks)
    keep = 1000
    for candidate in (1000, 700, 500, 350, 250, 180, 120, 80, 55, 35, 20, 12, 6, 0):
        if cost(candidate) <= budget:
            keep = candidate; break
    print(f"{len(picks)} frames, budget {budget:.1f} KB/frame -> ACTN2 hash keep "
          f"{keep}/1000 (mean {cost(keep):.1f} KB/frame)")"""
new = """    budget_kb = TARGET_KB/len(picks)
    seg_budget = budget_kb*1024*3/4/12          # base64 -> bytes -> segments
    # mean segment counts at full fidelity
    n_a2 = np.mean([len(geometry(h, f, box, la, lm, 1000)[0]) for f in probe])
    n_my = np.mean([len(geometry(h, f, box, la, lm, 1000)[1]) for f in probe])
    n_aa = np.mean([len(geometry(h, f, box, la, lm, 1000)[2]) for f in probe])
    share_a2 = min(1.0, 0.40*seg_budget/max(n_a2, 1))
    share_my = min(1.0, 0.25*seg_budget/max(n_my, 1))
    left = seg_budget - share_a2*n_a2 - share_my*n_my
    keep = int(np.clip(1000*left/max(n_aa, 1), 0, 1000))
    rng = np.random.default_rng(20260805)
    n_actin = h["actin/center"].shape[1]; n_myosin = h["myosin/center"].shape[1]
    keep_actin = np.zeros(n_actin, bool)
    keep_actin[rng.choice(n_actin, max(1, int(share_a2*n_actin)), replace=False)] = True
    keep_myo = np.zeros(n_myosin, bool)
    keep_myo[rng.choice(n_myosin, max(1, int(share_my*n_myosin)), replace=False)] = True
    print(f"{len(picks)} frames, {budget_kb:.1f} KB/frame = {seg_budget:.0f} segments")
    print(f"  full fidelity would need {n_a2+n_my+n_aa:.0f} segments/frame")
    print(f"  keeping actin {100*share_a2:.0f}%, myosin {100*share_my:.0f}%, "
          f"ACTN2 {keep/10:.0f}%")"""
if old not in s:
    sys.exit("allocator anchor not found")
s = s.replace(old, new, 1)
s = s.replace("        a2, my, aa, st, cm, ca, da, counts = geometry(h, f, box, la, lm, keep)",
              "        a2, my, aa, st, cm, ca, da, counts = geometry(\n"
              "            h, f, box, la, lm, keep, keep_actin, keep_myo)")
s = s.replace('f"{keep/10:.0f}% of mature ACTN2 bridges shown"',
              'f"{100*share_a2:.0f}% of status-2 actin, {100*share_my:.0f}% of myosin, "\n'
              '                  f"{keep/10:.0f}% of mature ACTN2 shown"')
s = s.replace('"· nominal ON 0.75 · all status-2 actin and all myosin; "',
              '"· nominal ON 0.75 · showing "')
p.write_text(s)
print("allocator added")
