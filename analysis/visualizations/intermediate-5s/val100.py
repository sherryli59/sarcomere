import base64, json, os, re
import numpy as np
name = "intermediate-5s-traj-100frames.html"
h = open(name).read()
run = json.loads(re.search(r"const runs = (.+?);\n", h, re.S).group(1))[0]
fr, box = run["frames"], np.asarray(run["box"], float)
keys = re.findall(r"\{key:'([a-z0-9]+)'", h)
missing = sum(1 for f in fr for k in keys if k not in f["layers"])
bad = 0
extremes = []
for f in fr:
    for k in keys:
        q = np.frombuffer(base64.b64decode(f["layers"][k]), np.int16)
        if len(q) % 6:
            bad += 1
        if len(q):
            xyz = q.reshape(-1, 3) * (box/2) / 32767.0
            extremes.append(np.abs(xyz).max(axis=0))
ex = np.max(extremes, axis=0)
print(f"  {name}")
print(f"    {os.path.getsize(name)/1024:6.0f} KB   frames {len(fr)}   "
      f"t {fr[0]['timeMs']:.0f}-{fr[-1]['timeMs']:.0f} ms")
print(f"    layers {keys}   missing {missing}   malformed {bad}")
print(f"    int16 decode -> max |x,y,z| = {ex.round(3)} um  (box/2 = {(box/2).round(2)})")
print(f"    quantisation step = {(box/32767*1000).round(4)} nm per axis")
print(f"    decodeSegments present: {'decodeSegments' in h}   "
      f"bare theme vars {len(re.findall(r'var[(]--viz-series-[0-9][)]', h))}")
print(f"    final frame: {fr[-1]['actinStatusCounts']['2']} status-2 filaments exist, "
      f"{fr[-1]['linkStatusCounts']['2']} mature bridges exist")
