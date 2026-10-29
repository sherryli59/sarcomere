"""Example actin-myosin clusters from the production resume trajectories, one by one.

Cluster = connected component of the actin-myosin graph at the shared 2.0 s state
(frame 400) where an edge is an actin-myosin PAIR held by >= K_STRONG strong
heads.  (K=1 percolates: one component holds ~1264 of 1600 actins and 318 of 320
myosins, so it cannot isolate anything.)  Membership is fixed at frame 400 and
the same filaments are followed 2.0 -> 2.4 s every 5 ms in two arms (control,
lateral repulsion off), so the same cluster can be compared across arms.

Drawn per frame: member actins / myosins / barbed ends (unwrapped across the
periodic box along the cluster's own edges, then tracked frame to frame), the
weak and strong heads between members (geometry of compute_am_two_spring_head,
checked against state/am_head_strain_energy_current), ACTN2 links with at least
one member end, and the non-member actins on the far end of those links (faint).
Strip: median angle between bound member pairs.  Scalar 6 = nematic S of the
member actins."""
import os, sys, json, base64, collections, numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_minimal_wrap_core as core
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023"
ARMS = [("r27ctrl", "control"), ("r27nolat", "lateral repulsion off")]
REF, FIRST, LAST = 400, 400, 480
K_STRONG = 2
TARGETS = [30, 12, 8, 4, 2, 1]          # myosins per example cluster
QUANT = 4000.0
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")
def mimg(d, box): return d - box * np.round(d / box)
def nem(d):
    if len(d) < 2: return float("nan")
    d = U(d); Q = 1.5 * (d[:, :, None] * d[:, None, :]).mean(0) - 0.5 * np.eye(3); return float(np.linalg.eigvalsh(Q)[-1])

def find_clusters(f, nA, nM):
    hs = f["state/am_head_state_current"][REF]; hb = f["state/am_head_bound_actin_current"][REF].astype(int)
    mi, hi = np.nonzero(hs == 2); pc = collections.Counter(zip(hb[mi, hi].tolist(), mi.tolist()))
    edges = [(a, m) for (a, m), v in pc.items() if v >= K_STRONG]
    adj = collections.defaultdict(list)
    for a, m in edges: adj[("a", a)].append(("m", m)); adj[("m", m)].append(("a", a))
    seen, comps = set(), []
    for node in list(adj):
        if node in seen: continue
        stack, comp = [node], []
        seen.add(node)
        while stack:
            x = stack.pop(); comp.append(x)
            for y in adj[x]:
                if y not in seen: seen.add(y); stack.append(y)
        comps.append(comp)
    out = []
    for c in comps:
        A = sorted(i for s, i in c if s == "a"); M = sorted(i for s, i in c if s == "m")
        out.append(dict(actins=A, myosins=M, edges=[(a, m) for a, m in edges if a in set(A)]))
    return out, adj

def pick(clusters):
    chosen, used = [], set()
    for tgt in TARGETS:
        cand = [c for i, c in enumerate(clusters) if i not in used]
        best = min(range(len(cand)), key=lambda j: (abs(len(cand[j]["myosins"]) - tgt), -len(cand[j]["actins"])))
        idx = clusters.index(cand[best]); used.add(idx); chosen.append(clusters[idx])
    return chosen

def unwrap_ref(f, c, adj, box):
    """BFS along the cluster's own edges at frame REF; returns unwrapped centres."""
    ca = f["actin/center"][REF]; cm = f["myosin/center"][REF]
    raw = lambda n: ca[n[1]] if n[0] == "a" else cm[n[1]]
    start = ("m", c["myosins"][0]); pos = {start: raw(start).astype(float)}; q = [start]
    members = set([("a", a) for a in c["actins"]] + [("m", m) for m in c["myosins"]])
    while q:
        x = q.pop()
        for y in adj[x]:
            if y in members and y not in pos:
                pos[y] = pos[x] + mimg(raw(y) - raw(x), box); q.append(y)
    assert len(pos) == len(members), "BFS did not reach every member"
    return np.array([pos[("a", a)] for a in c["actins"]]), np.array([pos[("m", m)] for m in c["myosins"]])

def build(arm, label, c, cid, adj, group, return_raw=False):
    path = "%s/data/%s_seed20260805.h5" % (CAMP, arm)
    f = h5py.File(path, "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]]); La = float(a["actin_length"]); Lm = float(a["myosin_length"])
    hph = int(a["am_heads_per_half"]); anc = core.anchor_axial_table(hph, Lm, int(float(a.get("am_spokes_per_crown", 3))))
    rho_rest = float(a["am_preferred_centerline"]); slack = float(a.get("am_tether_slack", 0.0)); cap = float(a["am_head_radial_force_max"])
    k_rho = float(a["am_head_transverse_k"]) + float(a["am_tether_k"]) / max(1.0, float(a["am_tether_saturation_heads"]))
    k_ax_weak = float(a["am_motor_k"]) * 0.1
    A = np.array(c["actins"]); M = np.array(c["myosins"]); aidx = {int(x): i for i, x in enumerate(A)}; midx = {int(x): i for i, x in enumerate(M)}
    Au, Mu = unwrap_ref(f, c, adj, box)
    frames = list(range(FIRST, LAST + 1)); nF = len(frames)
    act_seg, myo_seg, barb = [], [], []
    var = {k: ([], []) for k in ("armWeak", "armStrong", "aa1", "aa2", "ctx")}
    scal = np.zeros((nF, core.NSCAL)); fixups = 0; val_err = 0.0; nval = 0; raw_err = 0.0; emax = 0.0
    prevA, prevM = None, None
    step = f["state/current_step"][:, 0]
    off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]; coord = f["catch_bond/local_link_coordinates"]
    for fp, t in enumerate(frames):
        ca = f["actin/center"][t]; cm = f["myosin/center"][t]
        ua = U(f["actin/direction"][t]); um = U(f["myosin/direction"][t])
        if prevA is None:
            assert t == REF; cA, cM = Au, Mu
        else:
            cA = prevA + mimg(ca[A] - rawA_prev, box); cM = prevM + mimg(cm[M] - rawM_prev, box)
        prevA, prevM, rawA_prev, rawM_prev = cA, cM, ca[A].copy(), cm[M].copy()
        act_seg.append(np.concatenate([cA - 0.5 * La * ua[A], cA + 0.5 * La * ua[A]], 1))
        myo_seg.append(np.concatenate([cM - 0.5 * Lm * um[M], cM + 0.5 * Lm * um[M]], 1))
        barb.append(cA - 0.5 * La * ua[A])
        # heads between members
        hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
        tg = f["state/am_head_target_axial_current"][t]; sp = f["state/am_head_stroke_progress_current"][t]; en = f["state/am_head_strain_energy_current"][t]
        nw = ns = 0; angs = []
        for j_local, j in enumerate(M):
            for h in np.flatnonzero((hs[j] == 1) | (hs[j] == 2)):
                i = int(hb[j, h])
                if i not in aidx: continue
                il = aidx[i]
                m_att = cM[j_local] + anc[h] * um[j]
                s = float(np.clip(tg[j, h], -0.5 * La, 0.5 * La))
                a_att = cA[il] + s * ua[i]
                conn = a_att - m_att
                if np.linalg.norm(conn) > 0.3:          # image mismatch: draw the min-imaged connector
                    conn = mimg(conn, box); a_att = m_att + conn; fixups += 1
                z = float(conn @ ua[i]); rho = float(np.linalg.norm(conn - z * ua[i]))
                ld = rho - rho_rest; ld = ld - slack if ld > slack else (ld + slack if ld < -slack else 0.0)
                sc = cap / k_rho; x = abs(ld / sc); lc = x - np.log(2.0) if x > 20 else np.log(np.cosh(x))
                e = cap * sc * lc + (0.5 * k_ax_weak * z * z if hs[j, h] == 1 else 0.5 * k_ax_weak * max(0.0, z - sp[j, h]) ** 2)
                val_err = max(val_err, abs(e - float(en[j, h])) - 1e-6 * abs(float(en[j, h]))); raw_err = max(raw_err, abs(e - float(en[j, h]))); emax = max(emax, abs(float(en[j, h]))); nval += 1
                key = "armStrong" if hs[j, h] == 2 else "armWeak"
                var[key][0].append(np.concatenate([m_att, a_att])); var[key][1].append(fp)
                if hs[j, h] == 2: ns += 1
                else: nw += 1
        # pair angles for bound member pairs
        bp = set()
        for j_local, j in enumerate(M):
            for h in np.flatnonzero((hs[j] == 1) | (hs[j] == 2)):
                if int(hb[j, h]) in aidx: bp.add((int(hb[j, h]), int(j)))
        angs = [np.degrees(np.arccos(min(1.0, abs(float(ua[i] @ um[j]))))) for i, j in bp]
        # ACTN2 links with >= 1 member end
        row = off[t]; o, cnt = int(row[1]), int(row[2]); n1 = n2 = 0
        if cnt > 0:
            ids = ident[o:o + cnt]; cs = coord[o:o + cnt]
            for (i, j, st), (si, sj) in zip(ids[:, [2, 3, 4]].astype(int), cs[:, :2]):
                im, jm = i in aidx, j in aidx
                if not (im or jm): continue
                if not im: i, j, si, sj = j, i, sj, si; im, jm = jm, im
                p0 = cA[aidx[i]] + si * ua[i]
                if jm: p1 = cA[aidx[j]] + sj * ua[j]
                else:
                    cj = cA[aidx[i]] + mimg(ca[j] - ca[i], box); p1 = cj + sj * ua[j]
                    var["ctx"][0].append(np.concatenate([cj - 0.5 * La * ua[j], cj + 0.5 * La * ua[j]])); var["ctx"][1].append(fp)
                if np.linalg.norm(p1 - p0) > 0.3: p1 = p0 + mimg(p1 - p0, box); fixups += 1
                key = "aa2" if st == 2 else "aa1"
                var[key][0].append(np.concatenate([p0, p1])); var[key][1].append(fp)
                if st == 2: n2 += 1
                else: n1 += 1
        scal[fp] = [0, 0, nw, ns, n1, n2, nem(ua[A]), step[t], np.median(angs) if angs else np.nan]
    assert val_err <= 1e-11 * 0.2, "head energy check failed: %.3e (tolerance as in build_minimal_wrap_core: 1e-11*scale + 1e-6*|ref|)" % val_err
    if return_raw:
        raw = []
        for fp in range(nF):
            d = dict(actin=act_seg[fp], myosin=myo_seg[fp], barbed=barb[fp])
            for key, (rows, fo) in var.items():
                sel = [r for r, x in zip(rows, fo) if x == fp]
                d[key] = np.array(sel) if sel else np.zeros((0, 6))
            raw.append(d)
    layers = {}
    def q16(v):
        v = np.asarray(v, float).reshape(-1); qv = np.round(v * QUANT); assert np.abs(qv).max(initial=0) < 32767; return b64(qv.astype(np.int16))
    layers["actin"] = dict(values=q16(np.array(act_seg)), perFrame=len(A), kind="segments")
    layers["myosin"] = dict(values=q16(np.array(myo_seg)), perFrame=len(M), kind="segments")
    layers["barbed"] = dict(values=q16(np.array(barb)), perFrame=len(A), kind="points")
    for key, (rows, fo) in var.items():
        fo = np.array(fo, int); offs = np.zeros(nF + 1, np.int32); offs[1:] = np.cumsum(np.bincount(fo, minlength=nF)) if len(fo) else 0
        order = np.argsort(fo, kind="stable") if len(fo) else []
        vals = np.array(rows)[order] if len(rows) else np.zeros((0, 6))
        layers[key] = dict(values=q16(vals) if len(rows) else "", offsets=b64(offs), kind="segments")
    for key in ("anchorOff", "anchorDet"):
        layers[key] = dict(values="", offsets=b64(np.zeros(nF + 1, np.int32)), kind="points")
    allp = np.concatenate([np.array(act_seg).reshape(-1, 3), np.array(myo_seg).reshape(-1, 3)])
    lo, hi = allp.min(0), allp.max(0); pad = 0.06 * np.maximum(hi - lo, 0.1)
    bounds = [(lo - pad).tolist(), (hi + pad).tolist()]
    sc = np.where(np.isnan(scal), 0.0, scal).astype(np.float32)
    n_edges = len(c["edges"])
    setup = ("<b>Cluster %s</b> (%s): %d myosins and %d actins joined at 2.0 s by %d actin–myosin pairs that each hold "
             "≥%d strong heads. The same filaments are followed from 2.0 to 2.4 s every 5 ms; heads and ACTN2 links "
             "come and go, so later frames can show members no longer gripped. Faint actins are ACTN2 partners outside the cluster."
             % (cid, label, len(M), len(A), n_edges, K_STRONG))
    run = dict(key="%s_%s" % (cid, arm), label=label, group="Cluster %s · %d myosins, %d actins" % (cid, len(M), len(A)),
               setupHtml=setup, frames=nF, nActin=len(A), nMyosin=len(M), headsPerHalf=hph, nHeadTotal=int(2 * hph * len(M)),
               frameMs=5.0, box=box.tolist(), bounds=bounds, boundsEngaged=bounds, layers=layers, scalars=b64(sc.reshape(-1)),
               quant=QUANT, validated=True)
    print("  cluster %s %-22s myosins %d actins %d | strip first/last %.1f/%.1f deg | heads checked %d (raw max err %.1e of energies up to %.2e) | image fixups %d"
          % (cid, label, len(M), len(A), scal[0, 8], scal[-1, 8], nval, raw_err, emax, fixups), flush=True)
    return (run, raw) if return_raw else run

if __name__ == "__main__":
    f0 = h5py.File("%s/data/r27ctrl_seed20260805.h5" % CAMP, "r", locking=False)
    nA = f0["actin/center"].shape[1]; nM = f0["myosin/center"].shape[1]
    clusters, adj = find_clusters(f0, nA, nM)
    print("k>=%d: %d clusters at 2.0 s; sizes (myosins) %s" % (K_STRONG, len(clusters), sorted((len(c["myosins"]) for c in clusters), reverse=True)[:12]))
    chosen = pick(clusters)
    runs = []
    for n, c in zip("ABCDEF", chosen):
        for arm, label in ARMS:
            runs.append(build(arm, label, c, n, adj, None))
    tpl = open(os.path.join(HERE, "cluster-traj-template.html"), encoding="ascii").read()
    out = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":"))).replace("__TITLE__", "Actin&ndash;myosin clusters, one at a time")
    out = "".join(ch if ord(ch) < 128 else "&#%d;" % ord(ch) for ch in out[:out.index("<script")]) + out[out.index("<script"):]
    i0 = out.index("<script"); out = out[:i0] + "".join(ch if ord(ch) < 128 else "\\u%04x" % ord(ch) for ch in out[i0:])
    assert out.isascii()
    dst = os.path.join(HERE, "cluster-traj.html"); open(dst, "w", encoding="ascii").write(out)
    print("wrote %s (%.2f MB)" % (dst, os.path.getsize(dst) / 1e6))
