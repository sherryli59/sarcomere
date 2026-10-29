"""Early read of the production head-geometry ladder at a matched time every
run has reached. Three arms, paired seed-for-seed:
  A baseline  hph=9  spokes=3 -> 3 crowns/half, 225 nm, 18 heads/myosin
  B spokes1   hph=9  spokes=1 -> 9 crowns/half,  75 nm, 18 heads/myosin
  C hph27     hph=27 spokes=3 -> 9 crowns/half,  75 nm, 54 heads/myosin
B-A isolates crown spacing at fixed drive; C-B isolates drive at fixed spacing.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
BARE, BEAR = 0.10, 0.5*(1.0-0.10)
TARGET_MS = 1500
NFR, STRONG = 3, 2

def anchor_axial_table(hph, Lm, spokes):
    spokes = max(1, int(spokes)); cc = max(1, (hph+spokes-1)//spokes)
    out = np.empty(2*hph)
    for hi in (0, 1):
        sgn = -1.0 if hi == 0 else 1.0
        for h in range(hph):
            out[hi*hph+h] = sgn*(0.5*BARE*Lm + ((h//spokes + 0.5)/cc)*BEAR*Lm)
    return out

FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A baseline hph9/spokes3", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B spokes1  hph9/spokes1", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27    hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]

for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        seed = p.split("seed")[1][:8]
        try:
            with h5py.File(p, "r", locking=False) as f:
                A = f["meta"].attrs
                hph = int(A["am_heads_per_half"]); spk = int(float(A.get("am_spokes_per_crown", 3)))
                Lm = float(A["myosin_length"]); anc = anchor_axial_table(hph, Lm, spk)
                cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
                cf = min(cf, f["actin/direction"].shape[0], f["state/am_head_state_current"].shape[0])
                st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
                k = int(np.searchsorted(st, TARGET_MS*100, side="right"))-1
                if k < NFR: print("  SKIP %s (only %d frames)" % (seed, cf)); continue
                fr = list(range(k-NFR+1, k+1))
                ad = np.asarray(f["actin/direction"][fr], float)
                u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
                md = np.asarray(f["myosin/direction"][fr], float)
                w = md/np.linalg.norm(md, axis=2, keepdims=True)
                hs = np.asarray(f["state/am_head_state_current"][fr])
                ha = np.asarray(f["state/am_head_bound_actin_current"][fr])
                Q = 1.5*np.einsum("ia,ib->ab", u[-1], u[-1])/u.shape[1] - 0.5*np.eye(3)
                S = float(np.linalg.eigvalsh(Q).max())
                links = float(np.asarray(f["state/aa_pairs_current_count"][k]).reshape(-1)[-1])
                nstr = int((hs == STRONG).sum()/len(fr))
                nbound = int(((hs == 1) | (hs == STRONG)).sum()/len(fr))
                th = []; ncrowns = collections.Counter(); npairs = 0
                for t in range(len(fr)):
                    mi, hi = np.nonzero((hs[t] == STRONG) & (ha[t] >= 0))
                    if mi.size == 0: continue
                    ai = ha[t][mi, hi]; stack = np.stack([ai, mi], 1)
                    grip = collections.defaultdict(list)
                    for a_, m_, h_ in zip(ai.tolist(), mi.tolist(), hi.tolist()):
                        grip[(a_, m_)].append(anc[h_])
                    for (a_, m_), ancs in grip.items():
                        npairs += 1
                        ncrowns[len(set(np.round(ancs, 9)))] += 1
                        th.append(np.degrees(np.arccos(min(1.0, abs(float(np.dot(u[t][a_], w[t][m_])))))))
                multi = 100.0*sum(v for kk, v in ncrowns.items() if kk >= 2)/max(1, npairs)
                rows.append(dict(seed=seed, hph=hph, spk=spk, S=S, links=links, nstr=nstr,
                                 nbound=nbound, th=float(np.median(th)) if th else float("nan"),
                                 multi=multi, npairs=npairs//max(1, len(fr)), steps=int(st[k])))
        except Exception as e:
            print("  ERR %s: %s" % (seed, e))
    if not rows: print("\n=== %-26s NO DATA" % tag); continue
    g = lambda key: np.array([r[key] for r in rows], float)
    print("\n=== %-26s n=%d  (hph=%d spokes=%d)" % (tag, len(rows), rows[0]["hph"], rows[0]["spk"]))
    print("    seed      S     ACTN2  bound  strong  pairs  med_th  multi-crown%%")
    for r in rows:
        print("    %-8s %.4f  %5.0f  %5d  %5d  %5d  %6.2f  %6.2f" %
              (r["seed"], r["S"], r["links"], r["nbound"], r["nstr"], r["npairs"], r["th"], r["multi"]))
    print("    MEAN     %.4f  %5.0f  %5d  %5d  %5d  %6.2f  %6.2f" %
          (g("S").mean(), g("links").mean(), g("nbound").mean(), g("nstr").mean(),
           g("npairs").mean(), g("th").mean(), g("multi").mean()))
    print("    sd(S)=%.4f" % (g("S").std(ddof=1) if len(rows) > 1 else 0.0))
