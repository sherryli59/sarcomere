#!/usr/bin/env python3
"""Standalone 100-frame viewer: every filament, every ACTN2 bridge.

No subsampling. Coordinates are int16 fixed-point over [-L/2, L/2] per axis,
which halves the payload at a quantisation step of 0.18 nm along x -- three
orders of magnitude below the 15 nm filament radius, so nothing is lost at
display resolution. The file is meant to be downloaded and opened locally,
where payload size is not a constraint.

Actin is split by cb_status and ACTN2 by link status so each can be toggled;
status-1 layers start hidden so the opening view is the mature network.
"""
import base64, importlib.util, json, os, re, sys
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import h5py, numpy as np

D = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
         "visualizations/intermediate-5s")
SRC = ("/scratch/users/xsli/sarcomere/intermediate_rest74_lateral_5s_20260902/"
       "data/A1600_M320_Lx6_intermediate_seed20260805_5000ms.h5")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
OUT = sys.argv[2] if len(sys.argv) > 2 else "intermediate-5s-full-100frames.html"

spec = importlib.util.spec_from_file_location("pipeline", D / "build-simulations.py")
pipeline = importlib.util.module_from_spec(spec); spec.loader.exec_module(pipeline)

def band(v, bx):
    w = np.asarray(v, float) % bx
    k = np.arange(1, 5, dtype=float)
    return float(np.max(np.abs(np.exp(2j*np.pi*w[:, None]*k[None, :]/bx).mean(0))))

def enc(seg, box):
    if not len(seg):
        return ""
    half = (box/2.0)[None, None, :]
    q = np.clip(np.round(np.asarray(seg).reshape(-1, 2, 3)/half*32767.0),
                -32767, 32767).astype(np.int16)
    return base64.b64encode(np.ascontiguousarray(q).tobytes()).decode("ascii")

with h5py.File(SRC, "r", locking=False) as h:
    usable = pipeline.committed_frames(h)
    m = h["meta"].attrs
    box = np.asarray([m["box_x"], m["box_y"], m["box_z"]], float)
    dt, la, lm = float(m["dt"]), float(m["actin_length"]), float(m["myosin_length"])
    picks = [int(x) for x in np.unique(np.linspace(0, usable-1, N).round().astype(int))]
    frames = []
    for n, f in enumerate(picks):
        ca = np.asarray(h["actin/center"][f], float)
        da = np.asarray(h["actin/direction"][f], float)
        cm = np.asarray(h["myosin/center"][f], float)
        dm = np.asarray(h["myosin/direction"][f], float)
        st = np.asarray(h["actin/cb_status"][f, :, 0], int)
        ah, mh = 0.5*la*da, 0.5*lm*dm
        g = {}
        for s, key in ((0, "actin0"), (1, "actin1"), (2, "actin2")):
            k = st == s
            g[key] = np.asarray(pipeline.clip_periodic_segments(
                (ca-ah)[k], (ca+ah)[k], box), np.float32).reshape(-1, 6)
        g["myosin"] = np.asarray(pipeline.clip_periodic_segments(
            cm-mh, cm+mh, box), np.float32).reshape(-1, 6)
        counts = {"1": 0, "2": 0}
        g["aa1"] = np.empty((0, 6), np.float32); g["aa2"] = np.empty((0, 6), np.float32)
        op, ip, cp = ("catch_bond/local_link_frame_offsets",
                      "catch_bond/local_link_identity",
                      "catch_bond/local_link_coordinates")
        if all(x in h for x in (op, ip, cp)):
            _, off, cnt = np.asarray(h[op][f], int)
            ident = np.asarray(h[ip][off:off+cnt], int)
            coord = np.asarray(h[cp][off:off+cnt], float)
            for s, key in ((1, "aa1"), (2, "aa2")):
                k = ident[:, 4] == s
                counts[str(s)] = int(np.count_nonzero(k))
                si, sc = ident[k], coord[k]
                if len(si):
                    s0 = ca[si[:, 2]] + sc[:, :1]*da[si[:, 2]]
                    raw = ca[si[:, 3]] + sc[:, 1:]*da[si[:, 3]]
                    e0 = s0 + (raw - s0 - box[None, :]*np.rint((raw-s0)/box[None, :]))
                    g[key] = np.asarray(pipeline.clip_periodic_segments(s0, e0, box),
                                        np.float32).reshape(-1, 6)
        step = int(np.asarray(h["state/current_step"][f]).ravel()[0])
        frames.append({
            "frame": f, "step": step, "timeMs": step*dt*1e3,
            "mOrder": band(cm[:, 0], box[0]),
            "zOrder": max(band(ca[:, 0]+0.5*la*da[:, 0], box[0]),
                          band(ca[:, 0]-0.5*la*da[:, 0], box[0])),
            "nematic": float(1.5*np.mean(da[:, 0]**2)-0.5),
            "actinStatusCounts": {str(s): int(np.count_nonzero(st == s)) for s in (0,1,2)},
            "linkStatusCounts": counts, "myosinCount": int(len(cm)),
            "layers": {k: enc(v, box) for k, v in g.items()},
            "layerCounts": {k: int(len(v)) for k, v in g.items()}})
        if n % 20 == 0:
            print(f"  frame {n+1}/{len(picks)}  t={frames[-1]['timeMs']:.0f} ms  "
                  f"{sum(len(v) for v in g.values())} segments", flush=True)

kb = sum(len(v) for fr in frames for v in fr["layers"].values())/1024
last = frames[-1]["layerCounts"]
print(f"{len(frames)} frames, payload {kb/1024:.1f} MB")
print(f"  final frame: " + ", ".join(f"{k} {v}" for k, v in last.items()))

run = {"key": "intermediate-5s", "label": "Intermediate · seed 20260805",
       "config": ("A1600 M320 Lx=6 · ACTN2 rest 73.6 nm + barbed-end lateral repulsion "
                  "· nominal ON 0.75 · every filament and every ACTN2 bridge, "
                  "no subsampling"),
       "jobId": "41728503", "stateLabel": "Completed", "runTime": "", "reason": "",
       "box": box.tolist(), "targetTimeMs": 5000.0, "frames": frames}

LAYERS = [("actin0", "actin (unlinked)", "#b8bfbc", 0.9, 0.30, "false"),
          ("actin1", "actin (status-1)", "#6baed6", 1.2, 0.50, "false"),
          ("actin2", "actin (status-2)", "#08306b", 1.7, 0.80, "true"),
          ("myosin", "myosin", "#f5a45b", 6.6, 0.92, "true"),
          ("aa1", "ACTN2 (status-1)", "#f2b8cd", 0.8, 0.40, "false"),
          ("aa2", "ACTN2 (mature)", "#e87ba4", 1.35, 0.85, "true")]
tpl = (D / "intermediate-5s-animation-template.html").read_text()
tpl = re.sub(r"\s*\{key:'[a-z0-9]+'[^\n]*\n", "\n", tpl)
body = "".join(f"      {{key:'{k}',label:'{lab}',color:'{c}',width:{w},"
               f"alpha:{a},kind:'segments',on:{on}}},\n"
               for k, lab, c, w, a, on in LAYERS).rstrip(",\n") + "\n"
tpl = re.sub(r"    const layers = \[\s*\n\s*\];\n",
             "    const layers = [\n" + body + "    ];\n", tpl)
# honour the per-layer default visibility
tpl = tpl.replace("const visible = new Map(layers.map(layer => [layer.key,true]));",
                  "const visible = new Map(layers.map(layer => "
                  "[layer.key,layer.on !== false]));")
tpl = tpl.replace('type="checkbox" checked>',
                  'type="checkbox" ${layer.on === false ? "" : "checked"}>')
out = D / OUT
out.write_text(tpl.replace("/*__RUNS__*/", json.dumps([run], separators=(",", ":"))))
print(f"  wrote {out.name}  {out.stat().st_size/1048576:.1f} MB")
