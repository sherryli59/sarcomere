#!/usr/bin/env python3
"""Interactive trajectory viewer: status-2 actin and myosin, no subsampling.

Dropping the other layers is what makes full fidelity affordable. The previous
animation had to subsample filaments because aa1 (status-1 ACTN2) alone was
~68% of the payload. With only actin2 + myosin every filament can be drawn.

Status-2 is cb_status == 2, the same definition the snapshot pipeline uses for
its actin2 layer and the same one visualize_traj.py's --actin_display cb uses.
"""
import base64, importlib.util, json, os, sys
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import h5py, numpy as np

D = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
         "visualizations/intermediate-5s")
SRC = ("/scratch/users/xsli/sarcomere/intermediate_rest74_lateral_5s_20260902/"
       "data/A1600_M320_Lx6_intermediate_seed20260805_5000ms.h5")
N_FRAMES = int(sys.argv[1]) if len(sys.argv) > 1 else 40

spec = importlib.util.spec_from_file_location("pipeline", D / "build-simulations.py")
pipeline = importlib.util.module_from_spec(spec); spec.loader.exec_module(pipeline)


def band_order(v, box_x):
    w = np.asarray(v, float) % box_x
    k = np.arange(1, 5, dtype=float)
    return float(np.max(np.abs(np.exp(2j*np.pi*w[:, None]*k[None, :]/box_x).mean(0))))


def encode(a):
    return base64.b64encode(np.ascontiguousarray(a, np.float32).tobytes()).decode("ascii")


with h5py.File(SRC, "r", locking=False) as h:
    usable = pipeline.committed_frames(h)
    m = h["meta"].attrs
    box = np.asarray([m["box_x"], m["box_y"], m["box_z"]], float)
    dt, la, lm = float(m["dt"]), float(m["actin_length"]), float(m["myosin_length"])
    picks = np.unique(np.linspace(0, usable - 1, N_FRAMES).round().astype(int))
    frames = []
    for f in picks:
        f = int(f)
        ca = np.asarray(h["actin/center"][f], float)
        da = np.asarray(h["actin/direction"][f], float)
        cm = np.asarray(h["myosin/center"][f], float)
        dm = np.asarray(h["myosin/direction"][f], float)
        st = np.asarray(h["actin/cb_status"][f, :, 0], int)
        ah, mh = 0.5*la*da, 0.5*lm*dm
        keep = st == 2
        a2 = np.asarray(pipeline.clip_periodic_segments(
            (ca - ah)[keep], (ca + ah)[keep], box), np.float32).reshape(-1, 6)
        my = np.asarray(pipeline.clip_periodic_segments(
            cm - mh, cm + mh, box), np.float32).reshape(-1, 6)
        step = int(np.asarray(h["state/current_step"][f]).ravel()[0])
        frames.append({
            "frame": f, "step": step, "timeMs": step*dt*1e3,
            "mOrder": band_order(cm[:, 0], box[0]),
            "zOrder": max(band_order(ca[:, 0]+0.5*la*da[:, 0], box[0]),
                          band_order(ca[:, 0]-0.5*la*da[:, 0], box[0])),
            "nematic": float(1.5*np.mean(da[:, 0]**2)-0.5),
            "actinStatusCounts": {str(s): int(np.count_nonzero(st == s)) for s in (0, 1, 2)},
            "linkStatusCounts": {"1": 0, "2": 0},
            "myosinCount": int(len(cm)),
            "layers": {"actin2": encode(a2), "myosin": encode(my)},
            "layerCounts": {"actin2": int(len(a2)), "myosin": int(len(my))},
        })

kb = sum(len(v) for fr in frames for v in fr["layers"].values())/1024
print(f"{usable} committed frames, {len(frames)} animated, payload {kb:.0f} KB")
print(f"  final frame: actin2 {frames[-1]['layerCounts']['actin2']} segments "
      f"({frames[-1]['actinStatusCounts']['2']} filaments), "
      f"myosin {frames[-1]['layerCounts']['myosin']}")

run = {"key": "intermediate-5s", "label": "Intermediate · seed 20260805",
       "config": ("A1600 M320 Lx=6 · ACTN2 rest 73.6 nm + barbed-end lateral "
                  "repulsion · nominal ON 0.75 · ALL status-2 actin and ALL "
                  "myosin, no subsampling"),
       "jobId": "41728503", "stateLabel": "Completed", "runTime": "", "reason": "",
       "box": box.tolist(), "targetTimeMs": 5000.0, "frames": frames}

tpl = (D / "intermediate-5s-animation-template.html").read_text()
# this viewer ships two layers; drop the toggles for the rest
import re
tpl = re.sub(r"\s*\{key:'(?:actin0|actin1|aa1|aa2|barbed)'[^\n]*\n", "\n", tpl)
tpl = tpl.replace("{key:'actin2',label:'status-2 actin'",
                  "{key:'actin2',label:'actin (status-2)'")
out = D / "intermediate-5s-status2-trajectory.html"
out.write_text(tpl.replace("/*__RUNS__*/", json.dumps([run], separators=(",", ":"))))
print(f"  wrote {out.name}  {out.stat().st_size/1024:.0f} KB")
print(f"  layers in template: {re.findall(chr(123)+chr(107)+'ey:.([a-z0-9]+).', tpl)}")
