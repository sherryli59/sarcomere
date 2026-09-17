#!/usr/bin/env python3
"""Interactive trajectory viewer: status-2 actin, myosin, ACTN2 bridges.

Colours match the PyVista render so the two are comparable:
  status-2 actin  #08306b  (top of plt.cm.Blues, what --actin_color uniform gives)
  myosin          #f5a45b  (the script's bonded-myosin orange)
  ACTN2           pink
Nothing is subsampled; the frame count is chosen to hit a payload target,
because ACTN2 roughly triples the bytes per frame.
"""
import base64, importlib.util, json, os, re, sys
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import h5py, numpy as np

D = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
         "visualizations/intermediate-5s")
SRC = ("/scratch/users/xsli/sarcomere/intermediate_rest74_lateral_5s_20260902/"
       "data/A1600_M320_Lx6_intermediate_seed20260805_5000ms.h5")
TARGET_KB = float(sys.argv[1]) if len(sys.argv) > 1 else 800.0
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "intermediate-5s-traj-aa.html"

ACTIN2, MYOSIN, AA = "#08306b", "#f5a45b", "#e87ba4"

spec = importlib.util.spec_from_file_location("pipeline", D / "build-simulations.py")
pipeline = importlib.util.module_from_spec(spec); spec.loader.exec_module(pipeline)

def band(v, bx):
    w = np.asarray(v, float) % bx
    k = np.arange(1, 5, dtype=float)
    return float(np.max(np.abs(np.exp(2j*np.pi*w[:, None]*k[None, :]/bx).mean(0))))

def enc(a):
    return base64.b64encode(np.ascontiguousarray(a, np.float32).tobytes()).decode("ascii")

def build(h, f, box, la, lm):
    ca = np.asarray(h["actin/center"][f], float)
    da = np.asarray(h["actin/direction"][f], float)
    cm = np.asarray(h["myosin/center"][f], float)
    dm = np.asarray(h["myosin/direction"][f], float)
    st = np.asarray(h["actin/cb_status"][f, :, 0], int)
    ah, mh = 0.5*la*da, 0.5*lm*dm
    keep = st == 2
    a2 = np.asarray(pipeline.clip_periodic_segments((ca-ah)[keep], (ca+ah)[keep], box),
                    np.float32).reshape(-1, 6)
    my = np.asarray(pipeline.clip_periodic_segments(cm-mh, cm+mh, box),
                    np.float32).reshape(-1, 6)
    aa = np.empty((0, 6), np.float32); counts = {"1": 0, "2": 0}
    op, ip, cp = ("catch_bond/local_link_frame_offsets",
                  "catch_bond/local_link_identity", "catch_bond/local_link_coordinates")
    if all(n in h for n in (op, ip, cp)):
        _, off, cnt = np.asarray(h[op][f], int)
        ident = np.asarray(h[ip][off:off+cnt], int)
        coord = np.asarray(h[cp][off:off+cnt], float)
        for s in (1, 2):
            counts[str(s)] = int(np.count_nonzero(ident[:, 4] == s))
        sel = ident[:, 4] == 2
        si, sc = ident[sel], coord[sel]
        if len(si):
            s0 = ca[si[:, 2]] + sc[:, :1]*da[si[:, 2]]
            raw = ca[si[:, 3]] + sc[:, 1:]*da[si[:, 3]]
            e0 = s0 + (raw - s0 - box[None, :]*np.rint((raw-s0)/box[None, :]))
            aa = np.asarray(pipeline.clip_periodic_segments(s0, e0, box),
                            np.float32).reshape(-1, 6)
    step = int(np.asarray(h["state/current_step"][f]).ravel()[0])
    return {"frame": f, "step": step, "timeMs": step*float(h["meta"].attrs["dt"])*1e3,
            "mOrder": band(cm[:, 0], box[0]),
            "zOrder": max(band(ca[:, 0]+0.5*la*da[:, 0], box[0]),
                          band(ca[:, 0]-0.5*la*da[:, 0], box[0])),
            "nematic": float(1.5*np.mean(da[:, 0]**2)-0.5),
            "actinStatusCounts": {str(s): int(np.count_nonzero(st == s)) for s in (0,1,2)},
            "linkStatusCounts": counts, "myosinCount": int(len(cm)),
            "layers": {"actin2": enc(a2), "myosin": enc(my), "aa2": enc(aa)},
            "layerCounts": {"actin2": len(a2), "myosin": len(my), "aa2": len(aa)}}

with h5py.File(SRC, "r", locking=False) as h:
    usable = pipeline.committed_frames(h)
    m = h["meta"].attrs
    box = np.asarray([m["box_x"], m["box_y"], m["box_z"]], float)
    la, lm = float(m["actin_length"]), float(m["myosin_length"])
    # cost of the densest frame sets the frame budget
    # ACTN2 bridges accumulate over the run, so the last frame costs several
    # times the first. Sizing off the densest frame under-counts badly; sample
    # across the trajectory and use the mean.
    sample = [int(x) for x in np.linspace(0, usable-1, 5).round()]
    costs = [sum(len(v) for v in build(h, f, box, la, lm)["layers"].values())/1024
             for f in sample]
    per_kb = float(np.mean(costs))
    n = max(4, int(TARGET_KB/per_kb))
    print(f"frame cost {min(costs):.0f}-{max(costs):.0f} KB, mean {per_kb:.0f} "
          f"-> {n} frames for a {TARGET_KB:.0f} KB target")
    picks = np.unique(np.linspace(0, usable-1, n).round().astype(int))
    frames = [build(h, int(f), box, la, lm) for f in picks]

kb = sum(len(v) for fr in frames for v in fr["layers"].values())/1024
last = frames[-1]["layerCounts"]
print(f"  {len(frames)} frames, payload {kb:.0f} KB")
print(f"  final frame: actin2 {last['actin2']}, myosin {last['myosin']}, aa2 {last['aa2']}"
      f"  ({frames[-1]['linkStatusCounts']['2']} mature bridges, all drawn)")

run = {"key": "intermediate-5s", "label": "Intermediate · seed 20260805",
       "config": ("A1600 M320 Lx=6 · ACTN2 rest 73.6 nm + barbed-end lateral repulsion "
                  "· nominal ON 0.75 · all status-2 actin, all myosin, all mature ACTN2"),
       "jobId": "41728503", "stateLabel": "Completed", "runTime": "", "reason": "",
       "box": box.tolist(), "targetTimeMs": 5000.0, "frames": frames}

tpl = (D / "intermediate-5s-animation-template.html").read_text()
tpl = re.sub(r"\s*\{key:'(?:actin0|actin1|aa1|aa2|barbed|actin2|myosin)'[^\n]*\n", "\n", tpl)
layers = ("    const layers = [\n"
          f"      {{key:'actin2',label:'actin (status-2)',color:'{ACTIN2}',width:1.7,alpha:0.95,kind:'segments'}},\n"
          f"      {{key:'myosin',label:'myosin',color:'{MYOSIN}',width:4.8,alpha:0.90,kind:'segments'}},\n"
          f"      {{key:'aa2',label:'ACTN2 (mature)',color:'{AA}',width:1.35,alpha:0.85,kind:'segments'}}\n"
          "    ];\n")
tpl = re.sub(r"    const layers = \[\s*\n\s*\];\n", layers, tpl)
out = D / OUTNAME
out.write_text(tpl.replace("/*__RUNS__*/", json.dumps([run], separators=(",", ":"))))
print(f"  wrote {out.name}  {out.stat().st_size/1024:.0f} KB")
print(f"  template layers: {re.findall(chr(123)+'key:.([a-z0-9]+).', tpl)}")
