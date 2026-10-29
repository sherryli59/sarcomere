"""MP4 of any small rig via build_minimal_wrap_core.build_run (same geometry as its viewer).
tilted 30 deg about their barbed ends, 2 free thick filaments, 5 s, seed 20260805.
Geometry from build_minimal_wrap_core.build_run (the builder this rig's viewer uses);
unbound-head anchors are not drawn (not in the legend)."""
import os, sys, base64, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_minimal_wrap_core as core, render_mp4 as R


STRIDE = int(sys.argv[2]) if len(sys.argv) > 2 else 5
FPS = int(sys.argv[3]) if len(sys.argv) > 3 else 24
path = sys.argv[1]
r = core.build_run(path, stride=STRIDE)
dec = lambda s: np.frombuffer(base64.b64decode(s), np.float32) if s else np.zeros(0, np.float32)
frames = []
L = {k: v for k, v in r["layers"].items() if k in ("actin", "myosin", "barbed", "armWeak", "armStrong", "aa1", "aa2")}
store = {}
for k, spec in L.items():
    st = 3 if spec["kind"] == "points" else 6
    vals = dec(spec["values"]).reshape(-1, st)
    offs = np.frombuffer(base64.b64decode(spec["offsets"]), np.int32) if spec.get("offsets") else None
    store[k] = (vals, offs, spec.get("perFrame", 0))
for i in range(r["frames"]):
    d = {}
    for k, (vals, offs, pf) in store.items():
        d[k] = vals[i * pf:(i + 1) * pf] if pf else vals[offs[i]:offs[i + 1]]
    frames.append(d)
out = os.path.join(HERE, "videos", os.path.basename(path).replace(".h5", ".mp4")); os.makedirs(os.path.dirname(out), exist_ok=True)
R.render(frames, out, fps=FPS)
print("wrote", out, "%.2f MB" % (os.path.getsize(out) / 1e6), "frames", len(frames), "frameMs", r["frameMs"], "validated", r["validated"])
