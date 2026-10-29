"""Cluster E with its wider attachments.  On top of the cluster viewer (members fixed at the
shared 2.0 s state), each frame also draws
  * myoOther / oArmWeak / oArmStrong: myosins OUTSIDE the cluster that hold >= 1 bound head on a
    member actin, and those heads;
  * myoPartner / pArmStrong: myosins holding >= 1 STRONG head on an ACTN2 partner of a member
    (the actin across the Z-disc), and those heads (hidden by default).
Every drawn object is placed at the periodic image next to the member it attaches to.
Runs: control 2.0 -> 2.4 s, and kappa_aa 3 -> 30 from 2.0 s to its latest committed frame."""
import os, sys, json, base64, collections, numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_cluster_traj as B
import build_minimal_wrap_core as core
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def mimg(d, box): return d - box * np.round(d / box)
def q16(v):
    v = np.asarray(v, float).reshape(-1); q = np.round(v * B.QUANT); assert np.abs(q).max(initial=0) < 32767
    return base64.b64encode(q.astype(np.int16).tobytes()).decode("ascii")
ang = lambda a, b: float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b)))))))

def extras(arm, c, raw, first, last):
    f = h5py.File("%s/data/%s_seed20260805.h5" % (B.CAMP, arm), "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]]); La = float(a["actin_length"]); Lm = float(a["myosin_length"])
    hph = int(a["am_heads_per_half"]); anc = core.anchor_axial_table(hph, Lm, int(float(a.get("am_spokes_per_crown", 3))))
    A = list(c["actins"]); M = set(c["myosins"]); aidx = {x: k for k, x in enumerate(A)}
    lay = {k: ([], []) for k in ("myoOther", "oArmStrong", "myoPartner", "pArmStrong")}
    stats = []; fix = 0
    for fp, t in enumerate(range(first, last + 1)):
        seg = raw[fp]["actin"].reshape(-1, 6); cA = 0.5 * (seg[:, :3] + seg[:, 3:])          # unwrapped member centres
        ca = f["actin/center"][t]; ua = U(f["actin/direction"][t]); cm = f["myosin/center"][t]; um = U(f["myosin/direction"][t])
        hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int); tg = f["state/am_head_target_axial_current"][t]
        Mdir = U(np.sum([um[m] * np.sign(um[m] @ um[list(M)[0]]) for m in M], axis=0))
        def head_seg(jc, j, h, actin_c, i):
            m_att = jc + anc[h] * um[j]; a_att = actin_c + float(np.clip(tg[j, h], -0.5 * La, 0.5 * La)) * ua[i]
            if np.linalg.norm(a_att - m_att) > 0.3:
                a_att = m_att + mimg(a_att - m_att, box); return np.concatenate([m_att, a_att]), 1
            return np.concatenate([m_att, a_att]), 0
        # other myosins on member actins
        mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]          # STRONG heads only
        placed = {}; o_ang = []; n_strong = 0
        for j, h, i in zip(mi, hi, ai):
            if int(j) in M or int(i) not in aidx: continue
            il = aidx[int(i)]
            if j not in placed:
                placed[j] = cA[il] + mimg(cm[j] - ca[i], box)
                lay["myoOther"][0].append(np.concatenate([placed[j] - 0.5 * Lm * um[j], placed[j] + 0.5 * Lm * um[j]])); lay["myoOther"][1].append(fp)
                o_ang.append(ang(um[j], ua[i]))
            s_, k = head_seg(placed[j], j, h, cA[il], int(i)); fix += k; n_strong += 1
            lay["oArmStrong"][0].append(s_); lay["oArmStrong"][1].append(fp)
        # myosins with >= 1 STRONG head on an ACTN2 partner actin (the faint actins that are drawn)
        off = f["catch_bond/local_link_frame_offsets"][t]; o, cnt = int(off[1]), int(off[2]); ids = f["catch_bond/local_link_identity"][o:o + cnt]
        partners = {}
        for r in ids:
            x, y = int(r[2]), int(r[3])
            for i, j in ((x, y), (y, x)):
                if i in aidx and j not in aidx and j not in partners: partners[j] = cA[aidx[i]] + mimg(ca[j] - ca[i], box)
        pplaced = {}; n_pstrong = 0
        for j, h, i in zip(mi, hi, ai):
            if int(i) not in partners: continue
            pc = partners[int(i)]
            if int(j) in M and j not in pplaced:
                pass
            if j not in pplaced and j not in placed and int(j) not in M:
                pplaced[j] = pc + mimg(cm[j] - ca[i], box)
                lay["myoPartner"][0].append(np.concatenate([pplaced[j] - 0.5 * Lm * um[j], pplaced[j] + 0.5 * Lm * um[j]])); lay["myoPartner"][1].append(fp)
            jc = pplaced.get(j, placed.get(j))
            if jc is None:      # a cluster myosin gripping a partner actin: draw from the cluster myosin's unwrapped image
                jl = sorted(M).index(int(j)); mseg = raw[fp]["myosin"].reshape(-1, 6)[jl]; jc = 0.5 * (mseg[:3] + mseg[3:])
            s_, k = head_seg(jc, j, h, pc, int(i)); fix += k; n_pstrong += 1
            lay["pArmStrong"][0].append(s_); lay["pArmStrong"][1].append(fp)
        stats.append(dict(t_ms=t * 5.0, other_myosins=len(placed), other_strong_heads=n_strong,
                          partner_actins=len(partners), partner_myosins=len(pplaced), partner_strong_heads=n_pstrong,
                          other_myosin_to_its_actin=float(np.median(o_ang)) if o_ang else None))
    nF = last - first + 1; layers = {}
    for key, (rows, fo) in lay.items():
        fo = np.array(fo, int); offs = np.zeros(nF + 1, np.int32)
        if len(fo): offs[1:] = np.cumsum(np.bincount(fo, minlength=nF))
        order = np.argsort(fo, kind="stable") if len(fo) else []
        layers[key] = dict(values=q16(np.array(rows)[order]) if len(rows) else "",
                           offsets=base64.b64encode(offs.tobytes()).decode("ascii"), kind="segments")
    pts = [np.array(rows).reshape(-1, 3) for rows, _ in lay.values() if rows]
    return layers, stats, fix, pts

f0 = h5py.File("%s/data/r27ctrl_seed20260805.h5" % B.CAMP, "r", locking=False)
clusters, adj = B.find_clusters(f0, f0["actin/center"].shape[1], f0["myosin/center"].shape[1])
c = B.pick(clusters)[4]; assert (len(c["myosins"]), len(c["actins"])) == (2, 8)
g = h5py.File("%s/data/r27kaa30_seed20260805.h5" % B.CAMP, "r", locking=False); last_k = int(g["meta/committed_frames"][0]) - 1; g.close()
runs = []; allstats = {}
for arm, lab, last in (("r27ctrl", "control", 480), ("r27kaa30", "kappa_aa 3 -> 30", last_k)):
    B.LAST = last
    run, raw = B.build(arm, lab, c, "E", adj, None, return_raw=True)
    ext, st, fix, pts = extras(arm, c, raw, B.FIRST, last)
    run["layers"].update(ext)
    # camera: fit EVERYTHING drawn (members, ACTN2-partner actins, attached myosins, heads), not just members
    for fr in raw:
        for k in ("actin", "myosin", "ctx", "armWeak", "armStrong", "aa1", "aa2"):
            if len(fr[k]): pts.append(np.asarray(fr[k]).reshape(-1, 3))
    P = np.concatenate(pts); lo, hi = P.min(0), P.max(0); pad = 0.04 * np.maximum(hi - lo, 0.1)
    run["bounds"] = [(lo - pad).tolist(), (hi + pad).tolist()]; run["boundsEngaged"] = run["bounds"]
    run["setupHtml"] = ("<b>Cluster E</b> (%s), 2.000&ndash;%.3f s every 5 ms: the 2 cluster myosins and 8 actins (fixed at 2.0 s), "
                        "plus every OTHER myosin holding at least one STRONG head on one of those 8 actins, and every myosin holding a STRONG head on one of "
                        "their ACTN2 partners (the faint actins), each with those strong heads." % (lab, last * 0.005))
    runs.append(run); allstats[arm] = st
    om = [x["other_myosins"] for x in st]; oh = [x["other_strong_heads"] for x in st]; oa = [x["other_myosin_to_its_actin"] for x in st if x["other_myosin_to_its_actin"] is not None]
    pm = [x["partner_myosins"] for x in st]; pa = [x["partner_actins"] for x in st]
    print("%-18s frames %d | other myosins strongly on a cluster actin: median %d (range %d-%d), %d strong heads, %.1f deg to the actin | ACTN2-partner actins %d; other myosins strongly on them median %d (range %d-%d) | image fixups %d"
          % (lab, len(st), np.median(om), min(om), max(om), np.median(oh), np.median(oa), np.median(pa), np.median(pm), min(pm), max(pm), fix), flush=True)
tpl = open(os.path.join(HERE, "cluster-ext-template.html"), encoding="ascii").read()
out = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":"))).replace("__TITLE__", "Cluster E and every myosin strongly attached to it")
i0 = out.index("<script"); out = out[:i0] + "".join(ch if ord(ch) < 128 else "\\u%04x" % ord(ch) for ch in out[i0:])
assert out.isascii()
dst = os.path.join(HERE, "cluster-e-ext.html"); open(dst, "w", encoding="ascii").write(out)
json.dump(allstats, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/cluster_e_ext_stats.json", "w"))
print("wrote %s (%.2f MB)" % (dst, os.path.getsize(dst) / 1e6))
