import base64, json, re
import numpy as np
h = open("intermediate-5s-status2-20frames.html").read()
fr = json.loads(re.search(r"const runs = (.+?);\n", h, re.S).group(1))[0]["frames"]
for i, f in enumerate(fr):
    for k, v in f["layers"].items():
        a = np.frombuffer(base64.b64decode(v), np.float32)
        if len(a) == 0:
            print(f"  empty layer: frame index {i} (traj frame {f['frame']}, "
                  f"t={f['timeMs']:.0f} ms), layer '{k}', "
                  f"status-2 filaments = {f['actinStatusCounts']['2']}")
        elif len(a) % 6:
            print(f"  MALFORMED: frame {i} layer {k} len {len(a)}")
print("  key presence:", all(set(f["layers"]) == {"actin2", "myosin"} for f in fr))
print("  EMPTY_LAYER guard in file:", "EMPTY_LAYER" in h)
