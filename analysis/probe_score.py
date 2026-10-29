"""Pre-registered scoring for the confined-probe hph test.

LATCH CRITERION, fixed before seeing the hph=27 data: a run is latched if the
MEDIAN theta over the second half of the trajectory is < 10 deg. The confined
probe v2 campaign found a clean bimodal gap -- latched runs sit at theta median
2.7-6.0 deg, escaped runs at 13-68 deg -- so 10 deg is inside the gap. Per-seed
medians are printed so the threshold can be re-judged.

COUPLE: myosin/am_align_drive_{half} / am_align_head_count_{half}, which
excludes the am_alignment_k tether. Grouped BY HALF, never pooled: a span above
the 675 nm bearing length is an opposite-half pair and invalid.
"""
import os, sys, glob, json
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
KT = 4.114e-3          # pN.um
LATCH_DEG = 10.0
THETA_BINS = [(0, 10), (10, 30), (30, 50), (50, 70), (70, 90)]

def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v/np.where(n == 0, 1, n)

def score(p):
    with h5py.File(p, "r", locking=False) as h:
        if "committed_frames" not in h["meta"]:
            raise RuntimeError("in flight (no committed_frames yet)")
        cf = int(np.asarray(h["meta/committed_frames"]).reshape(-1)[-1])
        if cf < 200:
            raise RuntimeError("in flight (only %d frames)" % cf)
        cf = min(cf, h["actin/direction"].shape[0])
        u = unit(np.asarray(h["actin/direction"][:cf], float))[:, 0, :]
        w = unit(np.asarray(h["myosin/direction"][:cf], float))[:, 0, :]
        th = np.degrees(np.arccos(np.clip(np.abs(np.einsum("fk,fk->f", u, w)), 0, 1)))
        hs = np.asarray(h["state/am_head_state_current"][:cf])
        bound_frac = float(((hs == 1) | (hs == 2)).any(axis=(1, 2)).mean())
        half = cf//2
        out = dict(theta_med=float(np.median(th[half:])),
                   theta_med_all=float(np.median(th)),
                   bound_frac=bound_frac,
                   cap=float(np.asarray(h["meta"].attrs.get("cumulative_force_cap_events", np.nan)).reshape(-1)[-1])
                       if "cumulative_force_cap_events" in h["meta"].attrs else float("nan"))
        # couple, per half, BINNED BY THETA and over ALL frames. The drive is a
        # force-field property at a configuration, so using all frames is valid
        # and is the only way to get coverage at high theta in an arm that
        # latches quickly.
        for tag in ("minus", "plus"):
            d = np.asarray(h["myosin/am_align_drive_%s" % tag][:cf], float).reshape(cf, -1)[:, 0]
            n = np.asarray(h["myosin/am_align_head_count_%s" % tag][:cf], float).reshape(cf, -1)[:, 0]
            m = (n > 0) & (np.arange(cf) >= half)
            out["drive_%s_kT" % tag] = float((d[m]/n[m]).mean()/KT) if m.any() else float("nan")
            out["headfrac_%s" % tag] = float((n[half:] > 0).mean())
            per = d/np.where(n > 0, n, 1)/KT
            for lo, hi in THETA_BINS:
                sel = (n > 0) & (th >= lo) & (th < hi)
                out["bin_%s_%d_%d" % (tag, lo, hi)] = (float(per[sel].mean()), int(sel.sum())) \
                    if sel.any() else (float("nan"), 0)
        return out

ARMS = json.loads(sys.argv[1])
for label, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat), key=lambda q: int(q.split("seed")[-1].split(".")[0])):
        try: rows.append((int(p.split("seed")[-1].split(".")[0]), score(p)))
        except Exception as e: print("  ERR %s: %s" % (os.path.basename(p), e))
    if not rows:
        print("\n=== %-16s NO DATA YET" % label); continue
    th = np.array([r[1]["theta_med"] for r in rows])
    latch = th < LATCH_DEG
    print("\n=== %-16s n=%d   P(latch) = %d/%d = %.2f   [latch := median theta(2nd half) < %.0f deg]"
          % (label, len(rows), latch.sum(), len(rows), latch.mean(), LATCH_DEG))
    print("    per-seed median theta: " + " ".join("%d:%.1f%s" % (s, r["theta_med"], "*" if r["theta_med"] < LATCH_DEG else "") for s, r in rows))
    print("    bound-frame fraction : mean %.3f   min %.3f" % (
        np.mean([r[1]["bound_frac"] for r in rows]), np.min([r[1]["bound_frac"] for r in rows])))
    for tag in ("minus", "plus"):
        d = np.array([r[1]["drive_%s_kT" % tag] for r in rows]); d = d[np.isfinite(d)]
        hf = np.mean([r[1]["headfrac_%s" % tag] for r in rows])
        print("    drive/head %-5s     : mean %+.2f kT (n=%d seeds with heads)  frames-with-heads %.2f"
              % (tag, d.mean() if d.size else float("nan"), d.size, hf))
    caps = np.array([r[1]["cap"] for r in rows]); caps = caps[np.isfinite(caps)]
    if caps.size: print("    force-cap events     : mean %.0f" % caps.mean())
    print("    drive/head vs theta (minus half, kT; seed-mean, frames pooled across seeds):")
    for lo, hi in THETA_BINS:
        vals = []; nfr = 0
        for _, r in rows:
            v, c = r["bin_minus_%d_%d" % (lo, hi)]
            nfr += c
            if c > 0 and np.isfinite(v): vals.append(v)
        if vals:
            print("      theta [%2d,%2d): %+7.2f kT  (sd %5.2f over %2d seeds, %6d frames)"
                  % (lo, hi, np.mean(vals), np.std(vals, ddof=1) if len(vals) > 1 else 0.0, len(vals), nfr))
        else:
            print("      theta [%2d,%2d):     --     (no frames)" % (lo, hi))
