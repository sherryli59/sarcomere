import base64, json, os, re
import numpy as np
for name in ("intermediate-5s-traj-aa-800kb.html", "intermediate-5s-traj-aa-1600kb.html"):
    h = open(name).read()
    fr = json.loads(re.search(r"const runs = (.+?);\n", h, re.S).group(1))[0]["frames"]
    keys = re.findall(r"\{key:'([a-z0-9]+)'", h)
    cols = re.findall(r"color:'(#[0-9a-fA-F]{6})'", h)
    missing = sum(1 for f in fr for k in keys if k not in f["layers"])
    bad = sum(1 for f in fr for k in keys if k in f["layers"]
              and (lambda a: len(a) % 6 or not np.isfinite(a).all())(
                  np.frombuffer(base64.b64decode(f["layers"][k]), np.float32)))
    print(f"  {name}")
    print(f"    {os.path.getsize(name)/1024:6.0f} KB  frames {len(fr):3d}  "
          f"t {fr[0]['timeMs']:.0f}-{fr[-1]['timeMs']:.0f} ms")
    print(f"    layers {keys}  colours {cols}")
    print(f"    missing keys {missing}   malformed {bad}   "
          f"bare theme vars {len(re.findall(r'var[(]--viz-series-[0-9][)]', h))}")
