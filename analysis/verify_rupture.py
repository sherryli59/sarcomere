"""Does mechanical rupture QUANTITATIVELY account for arm C's 33% ACTN2 deficit?

A birth-death balance: at steady state  N = formation_rate / per_link_death_rate.
So if formation is unchanged and the per-link death rate rises by factor f, the
steady-state link count falls by 1/f. Measured N fell to 0.629/0.939 = 0.670,
i.e. 1/1.49 -- so the prediction is a death rate ~1.49x higher in C, and it must
come from reason 1 (connector_length > rupture_length), not from reason 2.

event_code 1 = formation, 2 = detachment, 3 = status change.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500          # ms
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B hph9/spokes1 ", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
out = {}
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            dt = float(a["dt"])
            cols = a["aa_local_link_events_v3_columns"].decode().split(",")
            ic = {c: i for i, c in enumerate(cols)}
            ev = f["catch_bond/local_link_events_v3"]
            step = np.asarray(ev[:, ic["step"]], float)
            lo, hi = T0/(dt*1e3), T1/(dt*1e3)     # ms -> steps
            m = (step >= lo) & (step < hi)
            if m.sum() == 0: continue
            code = np.asarray(ev[:, ic["event_code"]], int)[m]
            reas = np.asarray(ev[:, ic["detachment_reason"]], int)[m]
            dur_s = (T1-T0)/1000.0
            form = int((code == 1).sum())
            det = collections.Counter(reas[code == 2].tolist())
            ndet = int((code == 2).sum())
            # mean live link count over the same window
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            N = float(np.asarray(f["state/aa_pairs_current_count"][i0:i1]).reshape(i1-i0, -1)[:, -1].mean())
            rows.append(dict(form=form/dur_s, det=ndet/dur_s, N=N,
                             r1=det.get(1, 0)/dur_s, r2=det.get(2, 0)/dur_s,
                             r3=det.get(3, 0)/dur_s, r4=det.get(4, 0)/dur_s))
    if not rows: continue
    g = lambda k: float(np.mean([r[k] for r in rows]))
    N, form, det = g("N"), g("form"), g("det")
    out[tag[0]] = dict(N=N, form=form, det=det, r1=g("r1"), r2=g("r2"), r3=g("r3"),
                       k_death=det/N, k1=g("r1")/N, k2=g("r2")/N)
    print("\n=== %s  (n=%d seeds, %d-%d ms)" % (tag, len(rows), T0, T1))
    print("    mean live links N                 %8.1f" % N)
    print("    formation rate                    %8.1f /s" % form)
    print("    detachment rate (all reasons)     %8.1f /s" % det)
    print("      reason 1 mechanical rupture     %8.1f /s  (%.1f%% of detachments)" % (g("r1"), 100*g("r1")/det))
    print("      reason 2 stochastic             %8.1f /s  (%.1f%%)" % (g("r2"), 100*g("r2")/det))
    print("      reason 3 tip penetration        %8.1f /s  (%.1f%%)" % (g("r3"), 100*g("r3")/det))
    print("    PER-LINK death rate  k = det/N    %8.4f /s" % (det/N))
    print("      from reason 1                   %8.4f /s" % (g("r1")/N))
    print("      from reason 2                   %8.4f /s" % (g("r2")/N))
    print("    balance check  form/k = %8.1f  vs measured N = %8.1f  (ratio %.3f)"
          % (form/(det/N), N, form/(det/N)/N))
if "A" in out and "C" in out:
    A, C = out["A"], out["C"]
    print("\n=== C vs A")
    print("  link count ratio            N_C/N_A        = %.3f   (deficit %.1f%%)" % (C["N"]/A["N"], 100*(1-C["N"]/A["N"])))
    print("  formation rate ratio                       = %.3f" % (C["form"]/A["form"]))
    print("  PER-LINK death rate ratio   k_C/k_A        = %.3f" % (C["k_death"]/A["k_death"]))
    print("    reason-1 component        k1_C/k1_A      = %.3f" % (C["k1"]/A["k1"]))
    print("    reason-2 component        k2_C/k2_A      = %.3f" % (C["k2"]/A["k2"]))
    print("  PREDICTION if death-rate-driven: N ratio should be 1/%.3f = %.3f"
          % (C["k_death"]/A["k_death"], A["k_death"]/C["k_death"]))
    print("  times formation ratio %.3f -> predicted N ratio %.3f  vs OBSERVED %.3f"
          % (C["form"]/A["form"], (C["form"]/A["form"])*(A["k_death"]/C["k_death"]), C["N"]/A["N"]))
