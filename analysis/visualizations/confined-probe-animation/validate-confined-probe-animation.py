#!/usr/bin/env python3
"""VISUALIZATION_PIPELINE.md checklist items 3-7 and 9 for the animation artifacts.

Item 9's warning applies: the payload is a bare JSON list injected at
`const runs = /*__RUNS__*/`, so this asserts a nonzero run count FIRST and
refuses to pass vacuously.
"""
import base64, glob, json, os, re, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np, h5py

DATA = "/scratch/users/xsli/sarcomere/confined_probe_20260910/data"
B64 = re.compile(r"^[A-Za-z0-9+/]*={0,2}$")
failures, checks = [], 0


def check(name, ok, detail=""):
    global checks
    checks += 1
    if not ok:
        failures.append(f"{name}: {detail}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  ' + detail if detail else ''}")


def extract(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"const runs = (\[.*?\]);\n", text, re.S)
    assert m, f"{path}: could not find the injected `const runs = [...]` list"
    return text, json.loads(m.group(1))


for path in sorted(glob.glob("confined-probe-*.html")):
    if "template" in path:
        continue
    print(f"\n=== {path}")
    text, runs = extract(path)
    # Refuse to pass vacuously.
    check("nonzero run count", len(runs) > 0, f"{len(runs)} runs")
    if not runs:
        continue

    check("item 3: no unreplaced marker",
          "__RUNS__" not in text and "__TITLE__" not in text)
    check("item 5: no atob(", "atob(" not in text)
    bare = re.findall(r"var\(--viz-series-[0-9]\)", text) + \
           re.findall(r"var\(--foreground\)", text) + \
           re.findall(r"var\(--muted-foreground\)", text)
    check("item 5b: no fallback-less theme variable", not bare, f"{len(bare)} bare refs")
    for name in ("--sarc-series-3", "--sarc-series-4", "--sarc-series-5",
                 "--sarc-series-6", "--sarc-foreground", "--sarc-muted-foreground"):
        check(f"item 5b: artifact defines {name}", f"{name}:" in text)

    for run in runs:
        tag = run["key"]
        fields = ("actin", "myosin", "barbed", "scalars", "arms", "armState", "armOffset")
        check(f"item 4: {tag} base64-clean",
              all(B64.match(run[f]) for f in fields))
        n = run["frames"]
        dec = {f: np.frombuffer(base64.b64decode(run[f]),
                                dtype=np.int32 if f.startswith("arm") and f != "arms"
                                else np.float32) for f in fields}
        check(f"item 6: {tag} segment layers divisible by 6",
              len(dec["actin"]) == 6 * n and len(dec["myosin"]) == 6 * n
              and len(dec["arms"]) % 6 == 0,
              f"actin {len(dec['actin'])}, myosin {len(dec['myosin'])}, arms {len(dec['arms'])}")
        check(f"item 6: {tag} point layer divisible by 3", len(dec["barbed"]) == 3 * n)
        check(f"item 6: {tag} all values finite",
              all(np.isfinite(dec[f]).all() for f in
                  ("actin", "myosin", "barbed", "scalars", "arms")))
        check(f"item 6: {tag} arm offsets monotone and consistent",
              len(dec["armOffset"]) == n + 1
              and bool(np.all(np.diff(dec["armOffset"]) >= 0))
              and int(dec["armOffset"][-1]) == len(dec["arms"]) // 6
              and int(dec["armOffset"][-1]) == len(dec["armState"]))

        # Items 7 and 9: recompute from HDF5 independently of the builder.
        h5 = os.path.join(DATA, f"{run["config"]}_seed{run["seed"]}.h5")
        with h5py.File(h5, "r", locking=False) as h:
            committed = int(h["meta/committed_frames"][0])
            stride = max(1, -(-committed // n))          # ceil(committed / n)
            sl = slice(0, committed, stride)
            steps_ref = h["state/current_step"][sl, 0].astype(int)
            last_step = int(steps_ref[-1])
            box = np.array([h["meta"].attrs[k] for k in ("box_x", "box_y", "box_z")], float)
            ad = h["actin/direction"][sl, 0, :]
            md = h["myosin/direction"][sl, 0, :]
            ac = h["actin/center"][sl, 0, :]
            mc = h["myosin/center"][sl, 0, :]
            wk = h["actin/am_weak_count"][sl, 0, 0]
            sg = h["actin/am_strong_count"][sl, 0, 0]
        sc = dec["scalars"].reshape(n, 8)
        check(f"item 7: {tag} decimated frame count consistent",
              n == len(steps_ref), f"payload {n} vs stride-{stride} reference {len(steps_ref)}")
        check(f"item 7: {tag} every step matches the strided reference",
              np.array_equal(sc[:, 7].astype(int), steps_ref),
              f"stride {stride}, first mismatch "
              f"{int(np.argmax(sc[:,7].astype(int) != steps_ref)) if not np.array_equal(sc[:,7].astype(int), steps_ref) else -1}")
        check(f"item 7: {tag} last displayed step == HDF5 step",
              int(sc[-1, 7]) == last_step, f"{int(sc[-1,7])} vs {last_step}")
        theta = np.degrees(np.arccos(np.clip(np.abs(np.einsum(
            "fk,fk->f", ad / np.linalg.norm(ad, axis=1, keepdims=True),
            md / np.linalg.norm(md, axis=1, keepdims=True))), 0, 1)))
        check(f"item 9: {tag} theta matches independent recompute",
              float(np.abs(sc[:, 0] - theta).max()) < 2e-3,
              f"max dev {float(np.abs(sc[:,0]-theta).max()):.2e} deg")
        umd = md / np.linalg.norm(md, axis=1, keepdims=True)
        d = ac - mc
        d -= box * np.round(d / box)
        perp = np.linalg.norm(d - np.einsum("fk,fk->f", d, umd)[:, None] * umd, axis=1)
        check(f"item 9: {tag} perp matches independent recompute",
              float(np.abs(sc[:, 1] - perp).max()) < 2e-6,
              f"max dev {float(np.abs(sc[:,1]-perp).max()):.2e} um")
        check(f"item 9: {tag} head counts match actin/am_*_count",
              np.array_equal(sc[:, 3].astype(int), wk.astype(int))
              and np.array_equal(sc[:, 4].astype(int), sg.astype(int)))
        # Myosin drawn centred at the origin with its lab orientation.
        ms = dec["myosin"].reshape(n, 6)
        check(f"item 9: {tag} myosin held at the origin",
              float(np.abs(ms[:, :3] + ms[:, 3:]).max()) < 2e-6)

print(f"\n{checks} checks run, {len(failures)} failed")
for f in failures:
    print("  FAILED:", f)
print("\nNote: M order and Z order (checklist item 7) are not applicable to this rig -- "
      "one myosin and one actin make the band order trivially 1.0. The per-frame readout "
      "carries theta, weak/strong head counts, rho_perp, axial offset and arm length instead.")
sys.exit(1 if failures else 0)
