"""WHY does alignment reduce the head-count crosslink penalty by ~30%?

Repeat the formation decomposition in all four cells at matched time:
   formation_rate = eligible_pairs x gate_window_factor x residual
The k=0 half already gave 0.633 = 0.842 x 0.880 x 0.855 at 1000-1500 ms. If
alignment leaves eligible pairs unimproved (measured: interaction ~0), then its
crosslink gain must come from the gate window or the residual -- this says which.

Window 40-80 ms, matched across all four cells. TRANSIENT, so treat as a
mechanism comparison at matched time, not as steady-state rates.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 40, 80
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
AL   = "/scratch/users/xsli/sarcomere/amalign_20261005/data"
CELL = [("A   hph9  k=0 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("A9  hph9  k=30", f"{AL}/*kaa3_seed*_5000ms_align9.h5"),
        ("C   hph27 k=0 ", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5"),
        ("C27 hph27 k=30", f"{AL}/*kaa3_seed*_5000ms_align27.h5")]
res = {}
for tag, pat in CELL:
    rows = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            dt = float(a["dt"]); L = float(a["actin_length"]); bz = float(a["aa_barbed_zone"])
            cut = float(a["aa_preferred_centerline"]) + float(a["aa_capture_margin"])
            maxdot = float(a["aa_association_max_direction_dot"])
            box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
            cols = a["aa_local_link_events_v3_columns"].decode().split(",")
            ic = {c: i for i, c in enumerate(cols)}
            ev = f["catch_bond/local_link_events_v3"]
            step = np.asarray(ev[:, ic["step"]], float)
            lo, hi = T0/(dt*1e3), T1/(dt*1e3)
            m = (step >= lo) & (step < hi)
            code = np.asarray(ev[:, ic["event_code"]], int)[m]
            form = int((code == 1).sum())/((T1-T0)/1000.0)
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, max(i0+1, i1), 2))
            u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
            ac = np.asarray(f["actin/center"][fr], float)
            cr = np.asarray(f["actin/crosslink_ratio"][fr], float)
            win = np.minimum(bz, np.clip(cr, 0, 1)*L).mean()
            el = []
            for t in range(len(fr)):
                b = ac[t] - 0.5*L*u[t]
                d = b[:, None, :] - b[None, :, :]; d -= box*np.rint(d/box)
                dist = np.linalg.norm(d, axis=2); np.fill_diagonal(dist, 1e9)
                dot = np.einsum("ik,jk->ij", u[t], u[t])
                el.append(int(((dist < cut) & (dot < maxdot)).sum())//2)
            rows.append(dict(form=form, elig=float(np.mean(el)), win=float(win)))
    g = lambda k: float(np.mean([r[k] for r in rows]))
    res[tag.split()[0]] = dict(form=g("form"), elig=g("elig"), win=g("win"))
    print("=== %s  n=%d   formation %8.1f /s   eligible %7.1f   mean window %.4f um"
          % (tag, len(rows), g("form"), g("elig"), g("win")))
print("\n=== head effect decomposed, at each alignment setting")
for k0, k1, lab in (("A", "C", "k=0 "), ("A9", "C27", "k=30")):
    a, c = res[k0], res[k1]
    fr = c["form"]/a["form"]; er = c["elig"]/a["elig"]; wr = c["win"]/a["win"]
    rr = fr/(er*wr)
    print("  %s  formation %.3f = eligible %.3f x window %.3f x residual %.3f" % (lab, fr, er, wr, rr))
print("\n=== what alignment changed, within each head count")
for a_, b_, lab in (("A", "A9", "hph9 "), ("C", "C27", "hph27")):
    a, c = res[a_], res[b_]
    print("  %s  formation x%.3f   eligible x%.3f   window x%.3f   residual x%.3f"
          % (lab, c["form"]/a["form"], c["elig"]/a["elig"], c["win"]/a["win"],
             (c["form"]/a["form"])/((c["elig"]/a["elig"])*(c["win"]/a["win"]))))
