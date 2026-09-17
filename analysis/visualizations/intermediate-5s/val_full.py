import base64, json, os, re
import numpy as np
name = "intermediate-5s-full-100frames.html"
h = open(name).read()
run = json.loads(re.search(r"const runs = (.+?);\n", h, re.S).group(1))[0]
fr, box = run["frames"], np.asarray(run["box"], float)
keys = re.findall(r"\{key:'([a-z0-9]+)'", h)
defaults = dict(re.findall(r"\{key:'([a-z0-9]+)'.*?on:(true|false)\}", h))
missing = bad = 0
mx = np.zeros(3)
for f in fr:
    for k in keys:
        if k not in f["layers"]:
            missing += 1; continue
        q = np.frombuffer(base64.b64decode(f["layers"][k]), np.int16)
        if len(q) % 6:
            bad += 1
        if len(q):
            mx = np.maximum(mx, np.abs(q.reshape(-1, 3)*(box/2)/32767.0).max(axis=0))
print(f"  {name}")
print(f"    {os.path.getsize(name)/1048576:5.1f} MB   frames {len(fr)}   "
      f"t {fr[0]['timeMs']:.0f}-{fr[-1]['timeMs']:.0f} ms")
print(f"    layers {keys}")
print(f"    default visible {defaults}")
print(f"    missing keys {missing}   malformed {bad}   "
      f"bare theme vars {len(re.findall(r'var[(]--viz-series-[0-9][)]', h))}   "
      f"markers {len(re.findall(r'__[A-Z_]+__', h))}")
print(f"    int16 extremes {mx.round(3)} um vs box/2 {(box/2).round(2)}")
# cross-check drawn counts against the true totals in the readout
lf = fr[-1]
print(f"    final frame drawn: " + ", ".join(f"{k} {lf['layerCounts'][k]}" for k in keys))
print(f"    final frame真 totals: actin {lf['actinStatusCounts']}, "
      f"ACTN2 {lf['linkStatusCounts']}, myosin {lf['myosinCount']}")
