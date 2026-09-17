#!/usr/bin/env python3
import base64, json, os, re, sys
import numpy as np

for name in ("intermediate-5s-status2-20frames.html",
             "intermediate-5s-status2-40frames.html"):
    if not os.path.exists(name):
        print(f"  {name}: MISSING"); continue
    h = open(name).read()
    runs = json.loads(re.search(r"const runs = (.+?);\n", h, re.S).group(1))
    fr = runs[0]["frames"]
    keys = re.findall(r"\{key:'([a-z0-9]+)'", h)
    problems = 0
    for f in fr:
        for k in keys:
            if k not in f["layers"]:
                problems += 1
                continue
            a = np.frombuffer(base64.b64decode(f["layers"][k]), np.float32)
            if len(a) == 0 or len(a) % 6 or not np.isfinite(a).all():
                problems += 1
    bare = len(re.findall(r"var\(--viz-series-[0-9]\)", h))
    marker = len(re.findall(r"__[A-Z_]+__", h))
    print(f"  {name}")
    print(f"    {os.path.getsize(name)/1024:6.0f} KB   frames {len(fr)}   "
          f"t {fr[0]['timeMs']:.0f} to {fr[-1]['timeMs']:.0f} ms")
    print(f"    layers {keys}   bare theme vars {bare}   markers {marker}   "
          f"payload problems {problems}")
    print(f"    last frame: {fr[-1]['actinStatusCounts']['2']} status-2 filaments, "
          f"M {fr[-1]['mOrder']:.3f}  Z {fr[-1]['zOrder']:.3f}  S {fr[-1]['nematic']:.3f}")
