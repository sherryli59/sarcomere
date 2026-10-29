"""Trajectory viewer for the production resume A/B (hph27, 2.0 -> 2.4 s).
Reuses build_minimal_wrap_core.build_run unchanged (head geometry validated
against state/am_head_strain_energy_current), then, to fit a production-size
system in one page:
  * drops the unbound-head anchor layers (~15k points per frame),
  * stores coordinates as int16 at QUANT counts per um (0.25 nm resolution),
  * replaces the strip scalar (index 8) with the median angle between BOUND
    actin-myosin pairs, the quantity used throughout the analysis.
Setup text is built from each arm's recorded argv diff, not typed."""
import os, sys, json, base64, numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_minimal_wrap_core as core
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023"
QUANT = 4000.0
RUNS = [("r27ctrl", "control"), ("r27nolat", "lateral repulsion off"), ("r27lat4", "lateral cap 40 to 4 pN")]
FIRST, LAST, STRIDE = 400, 480, 8
def b64(a): return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")
def unb64(s, dt): return np.frombuffer(base64.b64decode(s), dtype=dt)
def pair_angle(f, t):
    u = f["actin/direction"][t]; u = u / np.linalg.norm(u, axis=1, keepdims=True)
    m = f["myosin/direction"][t]; m = m / np.linalg.norm(m, axis=1, keepdims=True)
    hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero((hs == 1) | (hs == 2)); k = np.unique(hb[mi, hi].astype(np.int64) * 1000 + mi)
    return float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.sum(u[k // 1000] * m[k % 1000], 1)), 0, 1)))))
def argv_diff(arm):
    src = open("%s/logs/hph27m24_seed20260805.argv" % CAMP).read().split("\n")
    new = open("%s/logs/%s_seed20260805.argv" % (CAMP, arm)).read().split("\n")
    s = {x.split("=", 1)[0]: x for x in src[1:] if x.strip()}; n = {x.split("=", 1)[0]: x for x in new[1:] if x.strip()}
    return [(k, s.get(k, "(absent)"), n.get(k, "(absent)")) for k in sorted(set(s) | set(n)) if s.get(k) != n.get(k) and k != "--filename"]
runs = []; maxerr = 0.0
for arm, label in RUNS:
    path = "%s/data/%s_seed20260805.h5" % (CAMP, arm)
    r = core.build_run(path, stride=STRIDE, step_range=(FIRST, LAST), label=label)
    for key in ("anchorOff", "anchorDet"):
        r["layers"][key] = dict(values="", offsets=b64(np.zeros(r["frames"] + 1, np.int32)), kind="points")
    for key, spec in r["layers"].items():
        if not spec["values"]: continue
        v = unb64(spec["values"], np.float32).astype(np.float64)
        q = np.round(v * QUANT); assert np.abs(q).max() < 32767, (arm, key, np.abs(v).max())
        maxerr = max(maxerr, float(np.abs(q / QUANT - v).max())) if len(v) else maxerr
        spec["values"] = b64(q.astype(np.int16))
    r["quant"] = QUANT
    sc = unb64(r["scalars"], np.float32).reshape(r["frames"], -1).copy()
    with h5py.File(path, "r", locking=False) as f:
        ts = list(range(FIRST, LAST + 1, STRIDE)); assert len(ts) == r["frames"]
        sc[:, 8] = [pair_angle(f, t) for t in ts]
        assert np.array_equal(sc[:, 7], f["state/current_step"][ts, 0].astype(np.float32)), "step column mismatch"
    r["scalars"] = b64(sc.astype(np.float32).reshape(-1))
    d = argv_diff(arm)
    rows = "".join("<tr><td><code>%s</code></td><td>%s</td><td><b>%s</b></td></tr>" % (k, o.split("=", 1)[-1], n.split("=", 1)[-1]) for k, o, n in d)
    r["group"] = "27 heads per half, resumed from the 2.0 s production state"
    r["setupHtml"] = ("<b>%s</b>: a copy of the production run <code>hph27m24_seed20260805</code>, resumed at 2.0 s (frame 400) "
                      "and run to 2.4 s with its original argv. Frames 0&ndash;400 are the shared history; shown here are "
                      "2.0&ndash;2.4 s every %d ms. Flags that differ from the source run (read from the recorded argv):"
                      "<table>%s</table>" % (label, int(STRIDE * 5), rows))
    runs.append(r)
    print(arm, "frames", r["frames"], "pair angle first/last %.1f/%.1f" % (sc[0, 8], sc[-1, 8]), flush=True)
tpl = open(os.path.join(HERE, "resab-traj-template.html"), encoding="ascii").read()
out = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":"))).replace("__TITLE__", "Production resume: lateral repulsion on, off, capped")
assert out.isascii()
dst = os.path.join(HERE, "resab-traj.html"); open(dst, "w", encoding="ascii").write(out)
print("wrote %s (%.2f MB), max quantisation error %.3f nm" % (dst, os.path.getsize(dst) / 1e6, 1000 * maxerr))
