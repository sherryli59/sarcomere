"""Distribution of barbed-end spacing within a Z-disc, per resume arm.
Spacing = transverse (y-z, min-image) distance from each barbed end to its
nearest barbed end assigned to the SAME Z-disc (Z-discs from the dominant
harmonic of the barbed-end x distribution, exactly as zdisc_capture in the
2026-08-31 factorial script).  Same metric as zcheck.py, but every value kept.
Frames: 400 (2.0 s, the shared resume point) and 460/470/480 (2.30-2.40 s, pooled)."""
import h5py, numpy as np, json, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
ARMS = [("r27ctrl", 27, "control"), ("r27nolat", 27, "lateral repulsion off"), ("r27lat4", 27, "lateral cap 40 to 4 pN"),
        ("r27soff16", 27, "strong dwell x10"), ("r81ctrl", 81, "control"), ("r81nolat", 81, "lateral repulsion off")]
def nn_values(f, t, box, L):
    c = f["actin/center"][t]; d = f["actin/direction"][t]; d = d / np.linalg.norm(d, axis=1, keepdims=True)
    barbed = c - 0.5 * L * d; bx = barbed[:, 0] % box[0]
    h = np.arange(1, 9, dtype=float); amps = np.abs(np.exp(2j * np.pi * bx[:, None] * h[None, :] / box[0]).mean(0))
    k = int(h[int(np.argmax(amps))]); phase = np.angle(np.exp(2j * np.pi * k * bx / box[0]).mean())
    zid = np.floor(((bx / box[0] * k) - phase / (2 * np.pi) + 0.5) % k).astype(int)
    out = []
    for z in range(k):
        idx = np.flatnonzero(zid == z)
        if len(idx) < 10: continue
        yz = barbed[idx][:, 1:]; dd = yz[:, None, :] - yz[None, :, :]; dd -= box[1:] * np.round(dd / box[1:])
        r = np.sqrt((dd ** 2).sum(-1)); np.fill_diagonal(r, np.inf); out.append(r.min(1))
    return k, 1000 * np.concatenate(out)
res = dict(meta=dict(lateral_range_nm=80.0, lateral_cap_saturates_below_nm=60.0, actin_diameter_nm=30.0,
                     actn2_rest_nm=None), arms=[])
for arm, hph, lab in ARMS:
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, arm), "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]]); L = float(a["actin_length"])
    if res["meta"]["actn2_rest_nm"] is None:
        for key in ("aa_rest_length", "actn2_rest_length", "aa_link_rest_length"):
            if key in a: res["meta"]["actn2_rest_nm"] = 1000 * float(a[key]); res["meta"]["actn2_rest_key"] = key; break
    k0, v0 = nn_values(f, 400, box, L)
    late = [nn_values(f, t, box, L) for t in (460, 470, 480)]
    v1 = np.concatenate([v for _, v in late])
    res["arms"].append(dict(arm=arm, hph=hph, label=lab, zdiscs=k0, t0_s=2.0, t1_s="2.30-2.40",
                            at_resume=[round(float(x), 1) for x in v0], late=[round(float(x), 1) for x in v1]))
    q = lambda v: "median %.1f, p10 %.1f, p90 %.1f, <30 nm %.0f%%, <60 nm %.0f%%" % (np.median(v), np.percentile(v, 10), np.percentile(v, 90), 100 * np.mean(v < 30), 100 * np.mean(v < 60))
    print("%-10s Z-discs %d | 2.0 s: %s | 2.3-2.4 s: %s" % (arm, k0, q(v0), q(v1)), flush=True)
print(res["meta"])
json.dump(res, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/spacing_dist.json", "w"))
