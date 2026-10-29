#!/usr/bin/env python3
"""Animated viewer for every gap-hunt minimal-rig trajectory.

Reuses build_minimal_wrap_core.build_run (the canonical minimal-wrap builder,
copied unchanged) for all geometry, layers and the head-energy validation.
Adds three things per run, none of them typed by hand:
  setupHtml.params  -- every flag that differs from the production argv,
                       parsed from the run's recorded logs/<arm>_seed<s>.argv
  setupHtml.result  -- tilt, heads, ACTN2, myosin spread and escape distance,
                       computed from the HDF5 itself
  setupHtml.prose   -- the one thing that is written: what the rung is FOR.
"""
import os, sys, json, math, subprocess, argparse
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_minimal_wrap_core import build_run

CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023"
PRODARGS = "/scratch/users/xsli/sarcomere/burialzone_20261020/prodargs.sh"
TARGET_FRAMES = 200
SKIP = {"filename", "seed", "nsteps", "save_every", "kill-on-bad-exit"}

G1 = "1 actin + 1 myosin (pose = standalone)"
G2 = "Partner ladder: N frozen parallel myosins"
G3 = "Scaffold disorder: 6 frozen myosins, translation off"
G4 = "Idealisation ladder at sigma = 0"
G5 = "Dwell sweep, FREE myosins"
G6 = "ACTN2 attempt v1 (geometry failed)"

RUNS = [
 (G1, "d_base", 20260805, "diag: production gates",
  "One actin at 36&deg; to +x, centre 80 nm from one myosin and over the nine crowns of its +x half. "
  "Myosin frozen. 200 ms diagnostic. Heads bind and the actin turns parallel."),
 (G1, "d_flip", 20260805, "diag: polarity flipped",
  "Same pose with the actin direction vector reversed. Heads in a myosin half bind only one actin "
  "polarity, so not one head binds and the actin simply diffuses."),
 (G1, "d_nothick", 20260805, "diag: thick feedback off",
  "Same as the production-gates diagnostic with <code>am_thick_feedback=false</code>. The thick-filament "
  "feedback turns out to RECRUIT heads: switching it off leaves fewer bound and alignment worse."),
 (G1, "d_noazi", 20260805, "diag: azimuth gate off",
  "Same with <code>am_binding_azimuth_selection=false</code>. Identical outcome: the binding azimuth gate is "
  "inert in this geometry."),
 (G1, "base", 20260805, "1 s, seed 805 (stays bound)",
  "One actin, one frozen myosin, 1 s, actin free to translate. This seed is caught by the heads and "
  "aligns."),
 (G1, "base", 20260815, "1 s, seed 815 (escapes)",
  "Identical configuration, different seed. The actin diffuses out of the 130 nm binding radius "
  "before the heads hold it and never comes back -- three of the four seeds did this. A lone actin beside "
  "a lone myosin is not held by anything else."),
 (G1, "nosteric", 20260805, "no A-M sterics",
  "As the 1 s run with <code>am_steric_k=0</code>, which the standalone script also lacks."),
 (G1, "freemyo", 20260805, "myosin free",
  "As the 1 s run but the myosin is not frozen."),
 (G2, "lat1par", 20260805, "1 myosin",
  "One actin with N myosins on a ring of radius 80 nm around it, all along +x and frozen; actin "
  "translation free. N = 1 here."),
 (G2, "lat2par", 20260805, "2 myosins", "N = 2 parallel frozen myosins around the actin."),
 (G2, "lat3par", 20260805, "3 myosins", "N = 3 parallel frozen myosins around the actin."),
 (G2, "lat4par", 20260805, "4 myosins", "N = 4 parallel frozen myosins around the actin."),
 (G2, "lat6par", 20260805, "6 myosins", "N = 6 parallel frozen myosins around the actin."),
 (G2, "lat3tilt", 20260805, "3 myosins, tilted",
  "Three frozen myosins with production direction noise (sigma 0.40, pairwise 21-45&deg;) instead of "
  "parallel."),
 (G3, "cf_par", 20260805, "sigma 0 (parallel)",
  "One actin inside six FROZEN myosins at 80 nm. Actin translation suppressed "
  "(<code>actin_crowding_parallel/perpendicular</code> 10 &rarr; 10000) so it cannot drift away; "
  "tumbling left at the production 160, so rotation is unchanged. Only the scaffold's direction spread "
  "varies across this group. Here the six are exactly parallel."),
 (G3, "cf_s0.20r11", 20260805, "sigma 0.20, realisation 11",
  "Same, six myosins drawn with direction noise sigma 0.20 (median pairwise ~19&deg;)."),
 (G3, "cf_s0.20r12", 20260805, "sigma 0.20, realisation 12",
  "Same noise level, a different random scaffold."),
 (G3, "cf_s0.40r11", 20260805, "sigma 0.40, realisation 11",
  "Six myosins at production direction noise, sigma 0.40 (median pairwise ~36&deg;)."),
 (G3, "cf_s0.40r12", 20260805, "sigma 0.40, realisation 12",
  "Same noise level, a different random scaffold."),
 (G4, "id_A", 20260805, "A: frozen + no translation",
  "Six PARALLEL myosins, one actin. A is the fully idealised baseline: myosins frozen and actin "
  "translation suppressed. Each later rung removes one idealisation."),
 (G4, "id_B", 20260805, "B: translation free",
  "As A but actin translation is free (production crowding)."),
 (G4, "id_C", 20260805, "C: myosins FREE",
  "As A but the six myosins are free to move. Watch the scaffold: it starts exactly parallel and the "
  "single actin drags it apart."),
 (G4, "id_D", 20260805, "D: myosins + translation free",
  "Both idealisations removed."),
 (G4, "id_E", 20260805, "E: free + 2 actins",
  "As D with a second, antiparallel actin, intended to let ACTN2 form. It did not: the two actins were "
  "placed on the same axis, so no link was ever possible. Read as a 2-actin run, not an ACTN2 test."),
 (G4, "id_F", 20260805, "F: frozen + 2 actins",
  "As A with the same second actin and the same pose bug -- no ACTN2 link formed."),
 (G5, "lf_ks40", 20260805, "k_strong_off 40 (25 ms dwell)",
  "id_A-type geometry but with the six myosins FREE and actin translation suppressed, sweeping only "
  "<code>am_k_strong_off</code> to test whether head dwell matters. 25 ms mean strong dwell."),
 (G5, "lf_ks160", 20260805, "k_strong_off 160 (6.25 ms, production)",
  "Same, production dwell."),
 (G5, "lf_ks640", 20260805, "k_strong_off 640 (1.6 ms)",
  "Same, four times faster turnover than production."),
 (G6, "aa_ks160", 20260805, "2 antiparallel actins, k_soff 160",
  "Two antiparallel actins in six FREE myosins, meant to sit 79 nm apart so ACTN2 could capture them. "
  "The 36&deg; tilt turned the axial offset into a 348 nm line separation, they drifted to ~900 nm, and "
  "no link ever formed. Kept so the failure is visible; the corrected pose is running."),
 (G6, "aa_ks27", 20260805, "2 antiparallel actins, k_soff 27",
  "Same failed pose with long-lived heads."),
]


G7 = "Tilt scan: 1 actin + 1 FREE myosin, strong heads"
_TS_PROSE = ("One actin and one FREE myosin. The actin starts at %d&deg; to the myosin, its centre 80 nm "
  "away along y and tilted in the plane perpendicular to that offset, so it never touches the myosin "
  "(closest approach 80 nm against the 65 nm contact distance). Heads bind for 100 ms with the filaments "
  "held at the pose (<code>fixed_lattice_preequil_steps=10000</code>), then both are released. Strong heads "
  "persist (<code>am_k_strong_off=0.1</code>, a 10 s dwell; no strain-dependent detachment) unless the "
  "model's geometric-rupture rule removes them, and detached heads can rebind. Production argv otherwise, "
  "<code>am_alignment_k=0</code>, 27 heads per half, motor saturation 24. The strip is the actin-myosin "
  "angle.")
for _T in (0, 15, 30, 45, 60, 75, 85):
    for _S in (20260805, 20260815, 20260825):
        RUNS.append((G7, "ts%d" % _T, _S, "%d deg, seed %s" % (_T, str(_S)[-3:]), _TS_PROSE % _T))

G8 = "Tilt scan, PRODUCTIVE polarity, pointed end outside the mid zone"
_TP_PROSE = ("One actin and one FREE myosin, actin starting at %d&deg; and pointing -x over the +x half, the "
  "productive polarity for that half: barbed end on the outer side, pointed end toward the M-line but outside "
  "the 0.075 um bare zone (pointed end at x = 0.65 - 0.5 cos&theta;). 80 nm off the myosin axis, tilted "
  "perpendicular to the offset, so it never touches it. Heads bind for 100 ms with the filaments held, then "
  "release; strong heads persist (10 s dwell, no strain-dependent detachment) unless the geometric-retention "
  "rule removes them. Production settings otherwise, <code>am_alignment_k=0</code>, motor saturation 24, 1 s.")
for _T in (0, 15, 30, 45, 60, 75, 85):
    for _S in (20260805, 20260815, 20260825):
        RUNS.append((G8, "tp%d" % _T, _S, "%d deg, seed %s" % (_T, str(_S)[-3:]), _TP_PROSE % _T))

G9 = "Tilt scan, productive polarity, Z-disc stop gated on ACTN2"
_TQ_PROSE = ("Identical to the productive-polarity tilt scan at %d&deg; except "
  "<code>aa_zdisc_tip_contact_requires_link=true</code>: the myosin capsule's Z-disc stop now acts only on an "
  "actin carrying a live ACTN2 link. This actin has none, so its bare barbed end can enter the A-band freely, "
  "as a thin filament should. Production leaves the flag off, so before this run the stop shoved uncrosslinked "
  "actins as if they were Z-discs.")
for _T in (0, 15, 30, 45, 60, 75, 85):
    for _S in (20260805, 20260815, 20260825):
        RUNS.append((G9, "tq%d" % _T, _S, "%d deg, seed %s" % (_T, str(_S)[-3:]), _TQ_PROSE % _T))

G10 = "One actin shared by 3 parallel myosins (vs 1)"
_TRI_PROSE = ("One actin shared by THREE parallel, FREE myosins in a trigonal cell (each myosin 80 nm from the "
  "actin axis, at azimuths 30, 150 and 270 deg). The actin starts at %d&deg;, pointing toward the M-line "
  "(productive), tilted along the gap between two myosins so its pointed half threads that gap (&ge; 69 nm "
  "clearance, contact 65 nm); its centre sits at x = 0.77 um so the barbed half lies beyond the myosins' ends. "
  "Pointed end outside the bare zone. Same settings as the gated tilt scan: heads bound 100 ms at the fixed "
  "pose, strong heads persist, Z-disc stop gated on ACTN2, 1 s.")
_MONO_PROSE = ("Matched control for the three-myosin cell: the SAME actin pose (centre x = 0.77 um, %d&deg;) "
  "with ONE free myosin, offset perpendicular to the tilt plane. Same settings, 1 s.")
for _T in (0, 15, 30, 45):
    RUNS.append((G10, "tri%d" % _T, 20260805, "3 myosins, %d deg" % _T, _TRI_PROSE % _T))
for _T in (0, 15, 30, 45):
    RUNS.append((G10, "mono%d" % _T, 20260805, "1 myosin, %d deg" % _T, _MONO_PROSE % _T))

G11 = "One actin shared by 3 TILTED myosins (sigma 0.40)"
_TRIS_PROSE = ("One actin shared by three FREE myosins drawn with the production direction-noise law at sigma "
  "0.40 (draw 241: 23.9, 26.3 and 4.9 deg to the axis, 25-50 deg apart pairwise), in the same trigonal cell as "
  "the parallel case. The actin starts at %d&deg;, productive polarity, pointed end outside the bare zone. A "
  "1.5 um myosin tilted that far crosses its neighbours, so as in production a sterics-only pre-phase "
  "(<code>pre_phase_steps=500</code>) resolves the overlaps before heads bind; the geometry at release is in "
  "the result line. Then heads bind for 100 ms with filaments held, strong heads persist, Z-disc stop gated "
  "on ACTN2, 1 s.")
for _T in (0, 15, 30, 45):
    RUNS.append((G11, "tris%d" % _T, 20260805, "sigma 0.4, %d deg" % _T, _TRIS_PROSE % _T))

def flags_of(lines):
    out = {}
    for a in lines:
        a = a.strip()
        if a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            out[k] = v
    return out


def production_flags():
    env = dict(os.environ, EXE="X", OUT="O", SEED="1", HPH="27", NSTEPS="1", SAVE_EVERY="1")
    raw = subprocess.run(["bash", "-c", "source %s; build_prod_cmd; printf '%%s\\n' \"${PROD_CMD[@]}\""
                          % PRODARGS], env=env, capture_output=True, text=True, check=True).stdout
    return flags_of(raw.split("\n"))


def result_of(path):
    with h5py.File(path, "r", locking=False) as f:
        n = int(min(f["meta/committed_frames"][0], f["actin/direction"].shape[0]))
        d = f["actin/direction"][:n]; c = f["actin/center"][:n]
        md = f["myosin/direction"][:n]; mc = f["myosin/center"][:n]
        h = (f["actin/am_strong_count"][:n, :, 0] + f["actin/am_weak_count"][:n, :, 0]).sum(axis=1)
        lk = f["actin/aa_local_link_count"][:n, :, 0].mean(axis=1)
        box = np.array([f["meta"].attrs["box_x"], f["meta"].attrs["box_y"], f["meta"].attrs["box_z"]], float)
    t0 = int(n * 0.6)
    u = d / np.linalg.norm(d, axis=2, keepdims=True)
    um = md / np.linalg.norm(md, axis=2, keepdims=True)
    tilt_x = np.degrees(np.arccos(np.clip(np.abs(u[t0:, :, 0]), 0, 1))).mean()
    # nearest myosin by centre distance, per frame per actin
    near = []
    for t in range(t0, n):
        for a in range(u.shape[1]):
            j = int(np.linalg.norm(mc[t] - c[t, a], axis=1).argmin())
            near.append(math.degrees(math.acos(min(1.0, abs(float(u[t, a] @ um[t, j]))))))
    last = n - 1
    sp = [math.degrees(math.acos(min(1.0, abs(float(um[last, i] @ um[last, k])))))
          for i in range(um.shape[1]) for k in range(i + 1, um.shape[1])]
    perp = []
    for a in range(u.shape[1]):
        best = 9e9
        for j in range(um.shape[1]):
            r = c[last, a] - mc[last, j]
            r = r - box * np.round(r / box)      # minimum image: the box is periodic
            best = min(best, float(np.linalg.norm(r - (r @ um[last, j]) * um[last, j])))
        perp.append(best)
    return dict(tilt_x=float(tilt_x), tilt_m=float(np.mean(near)), heads=float(h[t0:].mean()),
                links=float(lk[t0:].mean()), spread=float(np.median(sp)) if sp else 0.0,
                perp_nm=1000 * float(np.max(perp)))


def setup_html(prose, params, res, arm, seed):
    rows = "".join("<tr><td><code>%s</code></td><td>%s</td><td>%s</td></tr>" % (k, a, b)
                   for k, a, b in params)
    esc = res["perp_nm"] > 200
    # Result first and the flag table folded, so the canvas stays above the fold.
    return ("<strong>%s</strong> &middot; seed %d &middot; <code>%s</code>"
            "<div class='gh-result'><strong>Result</strong> (last 40%% of frames): tilt to +x %.1f&deg; "
            "&middot; to nearest myosin %.1f&deg; &middot; bound heads %.1f &middot; ACTN2 links %.2f "
            "&middot; myosin spread %.1f&deg; &middot; actin to myosin %.0f nm%s</div>"
            "<div>%s</div>"
            "<details><summary>%d flags differ from production (read from this run's argv)</summary>"
            "<table><tr><td><em>flag</em></td><td><em>production</em></td><td><em>this run</em></td></tr>%s"
            "</table></details>"
            % (arm, seed, "%s_seed%d.h5" % (arm, seed),
               res["tilt_x"], res["tilt_m"], res["heads"], res["links"], res["spread"], res["perp_nm"],
               " <strong>(escaped)</strong>" if esc else "", prose, len(params), rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", default=os.path.join(HERE, "gaphunt-rig-template.html"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups", nargs="*", default=None)
    ap.add_argument("--title", default="Gap-hunt minimal rigs")
    a = ap.parse_args()
    prod = production_flags()
    runs, summary = [], []
    for grp, arm, seed, label, prose in RUNS:
        if a.groups and grp not in a.groups:
            continue
        path = "%s/data/%s_seed%d.h5" % (CAMP, arm, seed)
        argv = flags_of(open("%s/logs/%s_seed%d.argv" % (CAMP, arm, seed)).read().split("\n"))
        params = []
        for k in sorted(set(prod) | set(argv)):
            if k in SKIP: continue
            pv, rv = prod.get(k, "(absent)"), argv.get(k, "(absent)")
            if pv != rv:
                if k == "initial_geometry_file": rv = os.path.basename(rv)
                params.append((k, pv, rv))
        with h5py.File(path, "r", locking=False) as f:
            nfr = int(min(f["meta/committed_frames"][0], f["actin/center"].shape[0]))
        stride = max(1, int(round(nfr / TARGET_FRAMES)))
        run = build_run(path, stride, label=label)
        run["key"] = "%s_%d" % (arm, seed)
        run["group"] = grp
        res = result_of(path)
        run["setupHtml"] = setup_html(prose, params, res, arm, seed)
        assert run["setupHtml"].isascii(), arm
        runs.append(run)
        summary.append((grp, arm, seed, label, res, len(params)))
        print("  %-12s seed %d  frames %d stride %d  -> %d drawn  tilt_x %.1f  heads %.1f  links %.2f"
              % (arm, seed, nfr, stride, run["frames"], res["tilt_x"], res["heads"], res["links"]))
    assert runs, "no runs built"
    for r in runs:
        assert r["layers"], r["key"]
    tpl = open(a.template, encoding="ascii").read()
    assert tpl.count("/*__RUNS__*/") == 1
    out = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":"))).replace("__TITLE__", a.title)
    assert out.isascii(), "output not ASCII-safe"
    open(a.out, "w", encoding="ascii").write(out)
    json.dump([dict(group=g, arm=ar, seed=s, label=l, **r, nflags=nf) for g, ar, s, l, r, nf in summary],
              open(a.out.replace(".html", ".summary.json"), "w"), indent=1)
    print("wrote %s  %.2f MB  %d runs" % (a.out, os.path.getsize(a.out) / 1e6, len(runs)))
