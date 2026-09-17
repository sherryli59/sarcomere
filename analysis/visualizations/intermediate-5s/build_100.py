#!/usr/bin/env python3
"""100-frame interactive viewer inside a payload budget.

Two levers, applied in this order:
  1. int16 fixed-point coordinates -- exactly 2x, and lossless at display
     resolution (0.09 nm in a 6 um box).
  2. subsample ACTN2 bridges -- they are the largest layer by far (3828 of
     5602 segments at t=5 s) and the only one that can be thinned without
     losing a species. Bridges are chosen by a hash of the actin PAIR, so the
     choice is stable in time and a persistent bridge does not flicker.
Actin and myosin are kept complete; ACTN2 is thinned only as far as the budget
demands, and the caption always reports the true totals.
"""
import base64, importlib.util, json, os, re, sys
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import h5py, numpy as np

D = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
         "visualizations/intermediate-5s")
SRC = ("/scratch/users/xsli/sarcomere/intermediate_rest74_lateral_5s_20260902/"
       "data/A1600_M320_Lx6_intermediate_seed20260805_5000ms.h5")
N_FRAMES = int(sys.argv[1]); TARGET_KB = float(sys.argv[2]); OUT = sys.argv[3]
ACTIN2, MYOSIN, AA = "#08306b", "#f5a45b", "#e87ba4"

spec = importlib.util.spec_from_file_location("pipeline", D / "build-simulations.py")
pipeline = importlib.util.module_from_spec(spec); spec.loader.exec_module(pipeline)

def band(v, bx):
    w = np.asarray(v, float) % bx
    k = np.arange(1, 5, dtype=float)
    return float(np.max(np.abs(np.exp(2j*np.pi*w[:, None]*k[None, :]/bx).mean(0))))

def quant(seg, box):
    """6 floats -> 6 int16, fixed point over [-L/2, L/2] per axis."""
    if not len(seg):
        return np.empty((0, 6), np.int16)
    half = (box / 2.0)[None, :]
    q = np.round(seg.reshape(-1, 2, 3) / half[:, None, :] * 32767.0)
    return np.clip(q, -32767, 32767).astype(np.int16).reshape(-1, 6)

def enc(a):
    return base64.b64encode(np.ascontiguousarray(a, np.int16).tobytes()).decode("ascii")

def geometry(h, f, box, la, lm, keep_pair, keep_actin=None, keep_myo=None):
    ca = np.asarray(h["actin/center"][f], float)
    da = np.asarray(h["actin/direction"][f], float)
    cm = np.asarray(h["myosin/center"][f], float)
    dm = np.asarray(h["myosin/direction"][f], float)
    st = np.asarray(h["actin/cb_status"][f, :, 0], int)
    ah, mh = 0.5*la*da, 0.5*lm*dm
    k2 = st == 2
    if keep_actin is not None:
        k2 = k2 & keep_actin
    a2 = np.asarray(pipeline.clip_periodic_segments((ca-ah)[k2], (ca+ah)[k2], box),
                    np.float32).reshape(-1, 6)
    km = slice(None) if keep_myo is None else keep_myo
    my = np.asarray(pipeline.clip_periodic_segments((cm-mh)[km], (cm+mh)[km], box),
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
            lo = np.minimum(si[:, 2], si[:, 3]).astype(np.int64)
            hi = np.maximum(si[:, 2], si[:, 3]).astype(np.int64)
            take = (((lo*73856093) ^ (hi*19349663)) % 1000) < keep_pair
            si, sc = si[take], sc[take]
        if len(si):
            s0 = ca[si[:, 2]] + sc[:, :1]*da[si[:, 2]]
            raw = ca[si[:, 3]] + sc[:, 1:]*da[si[:, 3]]
            e0 = s0 + (raw - s0 - box[None, :]*np.rint((raw-s0)/box[None, :]))
            aa = np.asarray(pipeline.clip_periodic_segments(s0, e0, box),
                            np.float32).reshape(-1, 6)
    return a2, my, aa, st, cm, ca, da, counts

with h5py.File(SRC, "r", locking=False) as h:
    usable = pipeline.committed_frames(h)
    m = h["meta"].attrs
    box = np.asarray([m["box_x"], m["box_y"], m["box_z"]], float)
    dt, la, lm = float(m["dt"]), float(m["actin_length"]), float(m["myosin_length"])
    picks = [int(x) for x in np.unique(np.linspace(0, usable-1, N_FRAMES).round().astype(int))]

    probe = [int(x) for x in np.linspace(0, usable-1, 5).round()]
    def cost(keep):
        tot = 0
        for f in probe:
            a2, my, aa, *_ = geometry(h, f, box, la, lm, keep)
            tot += (len(a2)+len(my)+len(aa)) * 12
        return tot/len(probe)/1024*4/3      # bytes -> KB of base64

    budget_kb = TARGET_KB/len(picks)
    seg_budget = budget_kb*1024*3/4/12          # base64 -> bytes -> segments
    # mean segment counts at full fidelity
    n_a2 = np.mean([len(geometry(h, f, box, la, lm, 1000)[0]) for f in probe])
    n_my = np.mean([len(geometry(h, f, box, la, lm, 1000)[1]) for f in probe])
    n_aa = np.mean([len(geometry(h, f, box, la, lm, 1000)[2]) for f in probe])
    share_a2 = min(1.0, 0.40*seg_budget/max(n_a2, 1))
    share_my = min(1.0, 0.25*seg_budget/max(n_my, 1))
    left = seg_budget - share_a2*n_a2 - share_my*n_my
    keep = int(np.clip(1000*left/max(n_aa, 1), 0, 1000))
    rng = np.random.default_rng(20260805)
    n_actin = h["actin/center"].shape[1]; n_myosin = h["myosin/center"].shape[1]
    keep_actin = np.zeros(n_actin, bool)
    keep_actin[rng.choice(n_actin, max(1, int(share_a2*n_actin)), replace=False)] = True
    keep_myo = np.zeros(n_myosin, bool)
    keep_myo[rng.choice(n_myosin, max(1, int(share_my*n_myosin)), replace=False)] = True
    print(f"{len(picks)} frames, {budget_kb:.1f} KB/frame = {seg_budget:.0f} segments")
    print(f"  full fidelity would need {n_a2+n_my+n_aa:.0f} segments/frame")
    print(f"  keeping actin {100*share_a2:.0f}%, myosin {100*share_my:.0f}%, "
          f"ACTN2 {keep/10:.0f}%")

    frames = []
    for f in picks:
        a2, my, aa, st, cm, ca, da, counts = geometry(
            h, f, box, la, lm, keep, keep_actin, keep_myo)
        step = int(np.asarray(h["state/current_step"][f]).ravel()[0])
        frames.append({
            "frame": f, "step": step, "timeMs": step*dt*1e3,
            "mOrder": band(cm[:, 0], box[0]),
            "zOrder": max(band(ca[:, 0]+0.5*la*da[:, 0], box[0]),
                          band(ca[:, 0]-0.5*la*da[:, 0], box[0])),
            "nematic": float(1.5*np.mean(da[:, 0]**2)-0.5),
            "actinStatusCounts": {str(s): int(np.count_nonzero(st == s)) for s in (0,1,2)},
            "linkStatusCounts": counts, "myosinCount": int(len(cm)),
            "layers": {"actin2": enc(quant(a2, box)), "myosin": enc(quant(my, box)),
                       "aa2": enc(quant(aa, box))},
            "layerCounts": {"actin2": len(a2), "myosin": len(my), "aa2": len(aa)}})

kb = sum(len(v) for fr in frames for v in fr["layers"].values())/1024
last = frames[-1]
print(f"  payload {kb:.0f} KB;  final frame actin2 {last['layerCounts']['actin2']} "
      f"myosin {last['layerCounts']['myosin']} aa2 {last['layerCounts']['aa2']} drawn "
      f"of {last['linkStatusCounts']['2']} mature")

run = {"key": "intermediate-5s", "label": "Intermediate · seed 20260805",
       "config": (f"A1600 M320 Lx=6 · ACTN2 rest 73.6 nm + barbed-end lateral repulsion "
                  f"· nominal ON 0.75 · showing "
                  f"{100*share_a2:.0f}% of status-2 actin, {100*share_my:.0f}% of myosin, "
                  f"{keep/10:.0f}% of mature ACTN2 shown"),
       "jobId": "41728503", "stateLabel": "Completed", "runTime": "", "reason": "",
       "box": box.tolist(), "targetTimeMs": 5000.0, "frames": frames}

tpl = (D / "intermediate-5s-animation-template.html").read_text()
tpl = re.sub(r"\s*\{key:'[a-z0-9]+'[^\n]*\n", "\n", tpl)
tpl = re.sub(r"    const layers = \[\s*\n\s*\];\n",
             "    const layers = [\n"
             f"      {{key:'actin2',label:'actin (status-2)',color:'{ACTIN2}',width:1.7,alpha:0.95,kind:'segments'}},\n"
             f"      {{key:'myosin',label:'myosin',color:'{MYOSIN}',width:4.8,alpha:0.90,kind:'segments'}},\n"
             f"      {{key:'aa2',label:'ACTN2 (mature)',color:'{AA}',width:1.35,alpha:0.85,kind:'segments'}}\n"
             "    ];\n", tpl)
out = D / OUT
out.write_text(tpl.replace("/*__RUNS__*/", json.dumps([run], separators=(",", ":"))))
print(f"  wrote {out.name}  {out.stat().st_size/1024:.0f} KB")
