"""Does a strong head ALIGN the actin, or only pin it?

A single attachment point cannot exert an aligning couple. Aligning an actin
toward the thick-filament axis needs two grips separated ALONG the filament.
Heads within one crown share an identical anchor_axial (and, with
am_head_azimuth=false, an identical 3-D anchor), so only a MULTI-CROWN grip can
make a couple. Crowns are 225 nm apart at hph=9 while the reach budget is 80 nm.

So: count distinct crowns per strong-head pair, and condition theta_AM on it.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
BARE, BEAR = 0.10, 0.5*(1.0-0.10)

def anchor_axial_table(hph, Lm, spokes=3):
    spokes = max(1, int(spokes)); cc = max(1, (hph+spokes-1)//spokes)
    out = np.empty(2*hph)
    for hi in (0, 1):
        sign = -1.0 if hi == 0 else 1.0
        for h in range(hph):
            out[hi*hph+h] = sign*(0.5*BARE*Lm + ((h//spokes + 0.5)/cc)*BEAR*Lm)
    return out

FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
K27  = "/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/data"
ARMS = [("kaa3_ksoff160", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("kaa3_ksoff27",  f"{K27}/*kaa3_seed*_5000ms_ksoff27.h5")]
TARGET, NFR, STRONG = 2000, 5, 2

for arm, pat in ARMS:
    nheads = collections.Counter(); ncrowns = collections.Counter()
    th_by_crowns = collections.defaultdict(list); spans = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            A = f["meta"].attrs
            hph = int(A["am_heads_per_half"]); Lm = float(A["myosin_length"])
            spokes = int(float(A.get("am_spokes_per_crown", 3)))
            anc = anchor_axial_table(hph, Lm, spokes)
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            k = int(np.searchsorted(st, TARGET*100, side="right"))-1
            if k < NFR: continue
            fr = list(range(k-NFR+1, k+1))
            ad = np.asarray(f["actin/direction"][fr], float)
            u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
            md = np.asarray(f["myosin/direction"][fr], float)
            w = md/np.linalg.norm(md, axis=2, keepdims=True)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr])
            for t in range(len(fr)):
                mi, hi = np.nonzero((hs[t] == STRONG) & (ha[t] >= 0))
                if mi.size == 0: continue
                ai = ha[t][mi, hi]
                grip = collections.defaultdict(list)
                for a_, m_, h_ in zip(ai.tolist(), mi.tolist(), hi.tolist()):
                    grip[(a_, m_)].append(anc[h_])
                for (a_, m_), ancs in grip.items():
                    nh = len(ancs); nc = len(set(np.round(ancs, 9)))
                    nheads[nh] += 1; ncrowns[nc] += 1
                    dot = abs(float(np.dot(u[t][a_], w[t][m_])))
                    th = float(np.degrees(np.arccos(min(1.0, dot))))
                    th_by_crowns[min(nc, 3)].append(th)
                    if nc >= 2: spans.append(max(ancs)-min(ancs))
    tot = sum(ncrowns.values())
    print("\n=== %s   (%d strong-head pairs)" % (arm, tot))
    print("  heads per pair : " + "  ".join("%d:%5.2f%%" % (k_, 100*v/tot) for k_, v in sorted(nheads.items())))
    print("  DISTINCT crowns: " + "  ".join("%d:%5.2f%%" % (k_, 100*v/tot) for k_, v in sorted(ncrowns.items())))
    for nc in sorted(th_by_crowns):
        v = np.array(th_by_crowns[nc])
        print("  theta_AM | %d crown(s): n=%7d  mean=%6.2f  median=%6.2f" % (nc, v.size, v.mean(), np.median(v)))
    if spans:
        sp = np.array(spans)
        print("  multi-crown axial span: n=%d  mean=%.4f um  min=%.4f  max=%.4f" % (sp.size, sp.mean(), sp.min(), sp.max()))
    else:
        print("  multi-crown axial span: NONE - no pair is gripped by two crowns")
