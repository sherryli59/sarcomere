"""What does removing the barbed-end lateral repulsion cost the Z-disc?  For each
resume arm at frames 400/440/480 (2.0/2.2/2.4 s): Z-disc capture (zdisc_capture,
imported unchanged from the 2026-08-31 factorial's analysis script), the median
transverse (y-z, min-image) nearest-neighbour distance between barbed ends
assigned to the same Z-disc, Z-band order of the barbed ends, and status-2 links."""
import h5py, numpy as np, json, os, sys, importlib.util
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
spec = importlib.util.spec_from_file_location("acp", "/scratch/users/xsli/sarcomere/actn2_rest74_lateral_20260831/source_snapshot/analyze_assembly_control_pilot.py")
acp = importlib.util.module_from_spec(spec); spec.loader.exec_module(acp)
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
out = {}
ARMS_Z = sys.argv[1].split(",") if len(sys.argv) > 1 else ["r27ctrl", "r27nolat", "r27lat4", "r27soff16", "r81ctrl", "r81nolat"]
OUTZ = sys.argv[2] if len(sys.argv) > 2 else "/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/zcheck.json"
for arm in ARMS_Z:
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, arm), "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]]); L = float(a["actin_length"])
    rows = []
    for t in (400, 440, 480):
        c = f["actin/center"][t]; d = f["actin/direction"][t]; d = d / np.linalg.norm(d, axis=1, keepdims=True)
        raw, pairs, mature, s1, s2 = acp.local_links(f, t)
        roots, counts, deg = acp.components(len(c), pairs)
        k, cap_mean, cap_max, sy, sz = acp.zdisc_capture(c, d, box, L, pairs, roots, counts)
        barbed = c - 0.5 * L * d
        bx = barbed[:, 0] % box[0]
        phase = np.angle(np.exp(2j * np.pi * k * bx / box[0]).mean())
        zid = np.floor(((bx / box[0] * k) - phase / (2 * np.pi) + 0.5) % k).astype(int)
        nn = []
        for z in range(k):
            idx = np.flatnonzero(zid == z)
            if len(idx) < 10: continue
            yz = barbed[idx][:, 1:]; dd = yz[:, None, :] - yz[None, :, :]; dd -= box[1:] * np.round(dd / box[1:])
            r = np.sqrt((dd ** 2).sum(-1)); np.fill_diagonal(r, np.inf); nn.append(r.min(1))
        nn = np.concatenate(nn)
        rows.append(dict(t_ms=t * 5.0, zdiscs=int(k), zcap=cap_mean, nn_median_nm=1000 * float(np.median(nn)),
                         nn_p10_nm=1000 * float(np.percentile(nn, 10)), z_order=acp.band_order(bx, box[0]),
                         status2=int(s2), links=int(s1 + s2)))
    out[arm] = rows
    print(arm, " | ".join("t %.0f zcap %.3f nnYZ %.1f (p10 %.1f) nm Zorder %.3f status2 %d links %d" % (r["t_ms"], r["zcap"], r["nn_median_nm"], r["nn_p10_nm"], r["z_order"], r["status2"], r["links"]) for r in rows), flush=True)
json.dump(out, open(OUTZ, "w"), indent=1)
