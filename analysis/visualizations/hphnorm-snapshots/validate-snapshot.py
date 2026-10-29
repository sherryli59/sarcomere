#!/usr/bin/env python3
"""Structural validation of the minimal-wrap viewer artifact + the physics answer."""
import base64, json, re, sys, os
import numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")

html = open(sys.argv[1], encoding="ascii").read()
checks = 0
def ck(cond, msg):
    global checks
    assert cond, "FAIL: " + msg
    checks += 1

ck(html.isascii(), "artifact is not ASCII-safe")
ck("/*__RUNS__*/" not in html, "payload placeholder not substituted")
ck("__TITLE__" not in html, "title placeholder not substituted")
m = re.search(r"const runs = (\[.*?\]);\n", html, re.S)
ck(m is not None, "cannot locate the runs payload")
runs = json.loads(m.group(1))
ck(len(runs) >= 1, "no runs in payload")

# every layer the template declares must exist in every run's payload
decl = re.findall(r"\{key:'([A-Za-z0-9]+)'", html[html.index("const layers = ["):])
decl = decl[:9]
ck(len(decl) == 9, f"expected 9 declared layers, got {len(decl)}")
kinds = dict(re.findall(r"\{key:'([A-Za-z0-9]+)'.*?kind:'(\w+)'", html[html.index("const layers = ["):]))

for r in runs:
    f = r["frames"]
    ck(f >= 1, f"{r["label"]}: frames={f}")  # snapshot: one frame is valid
    ck(r["validated"] is True, f"{r['label']}: head geometry NOT validated")
    for key in decl:
        ck(key in r["layers"], f"{r['label']}: template layer '{key}' missing from payload")
        spec = r["layers"][key]
        stride = 3 if kinds[key] == "points" else 6
        vals = np.frombuffer(base64.b64decode(spec["values"]), dtype=np.float32)
        if spec.get("perFrame"):
            ck(vals.size == f * spec["perFrame"] * stride,
               f"{r['label']}/{key}: values {vals.size} != frames*perFrame*stride")
        else:
            offs = np.frombuffer(base64.b64decode(spec["offsets"]), dtype=np.int32)
            ck(offs.size == f + 1, f"{r['label']}/{key}: offsets {offs.size} != frames+1")
            ck(offs[0] == 0 and np.all(np.diff(offs) >= 0),
               f"{r['label']}/{key}: offsets not non-decreasing from 0")
            ck(int(offs[-1]) * stride == vals.size,
               f"{r['label']}/{key}: offsets end {offs[-1]} inconsistent with {vals.size} floats")
        ck(np.all(np.isfinite(vals)), f"{r['label']}/{key}: non-finite coordinates")
    # Stride is read from the template, never hardcoded: it was hardcoded to 8
    # and a builder that grew a 9th scalar only failed because the reshape threw.
    m_stride = re.search(r"const NSCAL = (\d+)", html)
    ck(m_stride is not None, "template does not declare const NSCAL")
    nscal = int(m_stride.group(1))
    sc = np.frombuffer(base64.b64decode(r["scalars"]), dtype=np.float32).reshape(f, nscal)
    ck(np.all(np.isfinite(sc)), f"{r['label']}: non-finite scalars")
    # the strip series must be an angle in [0, 90], not a head count
    strip_ix = int(re.search(r"yAt\(values\.scalars\[i \* NSCAL \+ (\d+)\]\)", html).group(1))
    ck(strip_ix < nscal, f"strip series index {strip_ix} outside stride {nscal}")
    th = sc[:, strip_ix]
    ck(np.all((th >= -1e-6) & (th <= 90.0 + 1e-6)),
       f"{r['label']}: strip series ranges {th.min():.3f}-{th.max():.3f}, not a 0-90 deg angle")
    # the shading condition must key on weak+strong, not on an ACTN2 count
    shade = re.search(r"const bound = i < n \? values\.scalars\[i \* NSCAL \+ (\d+)\] "
                      r"\+ values\.scalars\[i \* NSCAL \+ (\d+)\]", html)
    ck(shade is not None and {int(shade.group(1)), int(shade.group(2))} == {2, 3},
       "strip shading must use scalars[2]+scalars[3] (weak+strong), not an ACTN2 count")
    tot = sc[:, 0] + sc[:, 1] + sc[:, 2] + sc[:, 3]
    ck(np.allclose(tot, r["nHeadTotal"]),
       f"{r['label']}: head states sum to {tot.min()}-{tot.max()}, not {r['nHeadTotal']}")

# THE CHECK THAT WAS MISSING: every run.<field> the template reads must exist
# in the payload.  The first version of this validator passed 264 checks on an
# artifact that rendered nothing, because run.boundsEngaged was read by the
# draw path and never emitted -- undefined threw on the first draw.
read = set(m.split(".", 1)[1] for m in re.findall(r"run\.[A-Za-z][A-Za-z0-9]*", html))
for r in runs:
    missing = sorted(read - set(r.keys()))
    ck(not missing, f"{r['label']}: template reads run.{{{','.join(missing)}}} "
                    f"but the payload does not provide it")
ck("boundsEngaged" in read, "sanity: expected the template to read boundsEngaged")


# Every CSS variable a layer references must be DEFINED in the artifact's own
# :root, not merely mentioned.  An undefined var resolves to the initial colour
# and the layer renders near-black -- which is what happened to the weak/strong
# head arms on the first build, and a whole-file grep for the name gave a false
# pass because the name appeared only in the layer declaration.
root_block = html[html.index(":root"):html.index("}", html.index(":root"))]
for var in sorted(set(re.findall(r"--sarc-series-\d+", html))):
    ck(f"{var}:" in root_block,
       f"layer references {var} but it is not defined in :root (would render black)")

print(f"PASS {checks} checks on {os.path.basename(sys.argv[1])}, {len(runs)} runs")

# ---------------- the physics answer ----------------
print()
print("Did the minimal wrap-around catch bond stabilize?  (5 s, current parameters)")
print(f"{'run':34s} {'off':>6s} {'det':>6s} {'weak':>6s} {'strong':>7s} "
      f"{'AA s1':>6s} {'AA s2':>6s} {'%frames s2>0':>13s} {'Morder':>7s}")
for r in runs:
    f = r["frames"]
    sc = np.frombuffer(base64.b64decode(r["scalars"]), dtype=np.float32).reshape(f, nscal)
    tail = sc[int(0.6 * f):]          # last 40% = steady state
    print(f"{r['label'][:34]:34s} {tail[:,0].mean():6.1f} {tail[:,1].mean():6.1f} "
          f"{tail[:,2].mean():6.2f} {tail[:,3].mean():7.2f} {tail[:,4].mean():6.2f} "
          f"{tail[:,5].mean():6.2f} {100*np.mean(tail[:,5]>0):13.1f} {tail[:,6].mean():7.3f}")
