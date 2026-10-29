"""1 actin + 1 free myosin from the tq30 pose: angle vs time and head occupancy by head kinetics."""
import h5py, numpy as np, os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from torque_an import anchors, unit
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
ARMS = ["tq30", "k01trap_30", "k16trap_30", "k160trap_30", "k160h81trap_30", "k160_30"]
SEEDS = ["20260805", "20260815", "20260825"]
out = {}
print("%-16s %-9s %6s %6s %6s %6s %6s %7s %7s %7s %8s %8s" % ("arm", "seed", "th0", "50ms", "100ms", "300ms", "1000ms", "t<10deg", "bound%", "heads", "span>=200", "perp_nm"))
for arm in ARMS:
    for s in SEEDS:
        p = D + "%s_seed%s.h5" % (arm, s)
        if not os.path.exists(p): continue
        f = h5py.File(p, "r", locking=False); a = dict(f["meta"].attrs)
        n = int(min(f["meta/committed_frames"][0], f["actin/center"].shape[0]))
        ua = unit(f["actin/direction"][:n, 0]); um = unit(f["myosin/direction"][:n, 0])
        th = np.degrees(np.arccos(np.clip(np.abs(np.sum(ua * um, 1)), 0, 1)))
        hs = f["state/am_head_state_current"][:n, 0]; b = (hs == 1) | (hs == 2)
        anc = anchors(int(a["am_heads_per_half"]), float(a["myosin_length"]))
        span = np.array([1000 * (anc[r].max() - anc[r].min()) if r.any() else -1 for r in b])
        box = np.array([a["box_x"], a["box_y"], a["box_z"]])
        d = f["actin/center"][:n, 0] - f["myosin/center"][:n, 0]; d -= box * np.round(d / box)
        perp = 1000 * np.linalg.norm(d - np.sum(d * um, 1)[:, None] * um, axis=1)
        dtf = float(a["dt"]) * float(np.median(np.diff(f["state/current_step"][:n, 0])))
        tt = np.arange(n) * dtf
        def at(ms): i = min(n - 1, int(round(ms / 1000 / dtf))); return th[i]
        below = np.nonzero(th < 10)[0]
        r = dict(th0=float(th[0]), th50=float(at(50)), th100=float(at(100)), th300=float(at(300)), th1000=float(at(1000)),
                 t10=float(tt[below[0]] * 1000) if len(below) else None, bound=float(b.any(1).mean()),
                 heads=float(b.sum(1)[b.any(1)].mean()) if b.any() else 0.0, span200=float((span >= 200).mean()),
                 perp=float(np.median(perp)), th_med_last=float(np.median(th[n // 2:])), n=n)
        out["%s_%s" % (arm, s)] = r
        print("%-16s %-9s %6.1f %6.1f %6.1f %6.1f %6.1f %7s %7.0f %7.1f %8.2f %8.0f" % (arm, s[-3:], r["th0"], r["th50"], r["th100"], r["th300"], r["th1000"],
              "%.0f" % r["t10"] if r["t10"] is not None else "never", 100 * r["bound"], r["heads"], r["span200"], r["perp"]))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/kin_an.json", "w"), indent=1)
