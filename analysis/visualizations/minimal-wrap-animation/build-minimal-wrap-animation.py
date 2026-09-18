#!/usr/bin/env python3
"""Interactive viewer for the MINIMAL WRAP-AROUND catch bond.

Follows VISUALIZATION_PIPELINE.md: single inline HTML, base64 float32 payload,
ASCII-safe, species-hue layer contract with artifact-owned colour fallbacks.

Rig: initial_structure=minimal_catch_balanced -- 4 actins, 2 myosins, periodic
x.  Actins (0,1) are the central Z-disc; (2,3) are its periodic image at the x
boundary, so the cell is one sarcomere repeat closing through the boundary.
Drawn in the LAB frame (not myosin-centred) because the wrap is the point.

Motor-head geometry is reconstructed exactly as compute_am_two_spring_head()
(interaction.cpp:555) and validated against state/am_head_strain_energy_current
before anything is written.  All four head states are drawn:
  off_srx(3) / detached(0)  -> anchor points on the myosin surface
  weak(1) / strong(2)       -> arm segments anchor -> actin attachment
"""
from __future__ import annotations
import argparse, base64, json, os, re, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np, h5py

BARE_FRAC = 0.10
BEAR_FRAC = 0.5 * (1.0 - BARE_FRAC)
DETACHED, WEAK, STRONG, OFF = 0, 1, 2, 3
NSCAL = 9


def encode(values, dtype=np.float32):
    return base64.b64encode(
        np.ascontiguousarray(values, dtype=dtype).tobytes()).decode("ascii")


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n > 1e-12, n, 1.0)


def soft_threshold(x, slack):
    return np.where(x > slack, x - slack, np.where(x < -slack, x + slack, 0.0))


def anchor_axial_table(hph, myosin_length, spokes=3):
    """head.anchor_axial exactly as _initialize_am_head_pools, sarcomere.cpp:2939.

    `spokes` is am_spokes_per_crown, read from the run's meta attrs (older files
    predate the flag and are 3 by definition).  Ceiling division, so spokes=3
    reproduces the previous hardcoded (hph + 2) // 3.
    """
    spokes = max(1, int(spokes))
    cc = max(1, (hph + spokes - 1) // spokes)
    out = np.empty(2 * hph)
    for half_index in (0, 1):
        half_sign = -1.0 if half_index == 0 else 1.0
        for h in range(hph):
            out[half_index * hph + h] = half_sign * (
                0.5 * BARE_FRAC * myosin_length
                + ((h // spokes + 0.5) / cc) * BEAR_FRAC * myosin_length)
    return out


def min_image(d, box, periodic):
    return d - np.where(periodic, box * np.round(d / box), 0.0)


# No clip_periodic_segments here.  Every segment this builder emits is already
# contiguous: filaments are centre +- half-length*direction, and ACTN2 links and
# head arms are min-imaged against their first endpoint.  None is ever drawn
# across the box, which is what the pipeline forbids.  Clipping by DROPPING
# boundary-crossing segments would delete actins 2 and 3 -- the periodic-image
# Z-disc that is the whole subject of this rig.


def band_order(x, box_x):
    """VISUALIZATION_PIPELINE.md band_order, verbatim."""
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(2j * np.pi * wrapped[..., None] * harmonics[None, None, :]
                   / box_x).mean(axis=1)
    return np.max(np.abs(phase), axis=1)


def read_material(handle, sl, nmyo):
    """Material frames live under either name depending on campaign vintage."""
    if "state/myosin_material_e1_current" in handle:
        return (handle["state/myosin_material_e1_current"][sl],
                handle["state/myosin_material_e2_current"][sl])
    e1 = np.stack([handle[f"myosin/material_e1_{c}"][sl, :, 0] for c in "xyz"], -1)
    e2 = np.stack([handle[f"myosin/material_e2_{c}"][sl, :, 0] for c in "xyz"], -1)
    return e1, e2


def build_run(path, stride=1, box_override=None, label=None):
    with h5py.File(path, "r", locking=False) as h:
        attrs = dict(h["meta"].attrs)
        frames = int(min(h["meta/committed_frames"][0],
                         h["actin/center"].shape[0],
                         h["state/am_head_state_current"].shape[0]))
        if "box_x" in attrs:
            box = np.array([attrs["box_x"], attrs["box_y"], attrs["box_z"]], float)
        elif "box" in h["meta"]:
            box = np.array(h["meta/box"][()], float).ravel()[:3]
        else:
            assert box_override is not None, f"{path}: no box in meta, pass --box"
            box = np.array(box_override, float)
        periodic = (np.array([attrs.get("periodic_x", 1), attrs.get("periodic_y", 1),
                              attrs.get("periodic_z", 1)], bool)
                    if "periodic_x" in attrs else np.array([True, True, True]))
        dt = float(attrs["dt"])
        La = float(attrs["actin_length"]); Lm = float(attrs["myosin_length"])
        Rm = float(attrs.get("myosin_radius", 0.05))
        hph = int(attrs["am_heads_per_half"])
        # Radial-channel mode.  Absent => an old file, which is legacy by
        # definition (surface anchor, 30 nm rest, unbounded linear spring).
        azimuth_on = bool(int(float(attrs.get("am_head_azimuth", 1))))
        # am_head_radial_force_max records the CONFIGURED cap, which the code
        # ignores in azimuth-ON mode (it passes 0 there).  Using the attribute
        # unconditionally applied a saturating energy to a linear legacy run and
        # broke the cross-check by 28.5%.
        radial_cap = (0.0 if azimuth_on
                      else float(attrs.get("am_head_radial_force_max", 0.0)))
        if azimuth_on:
            rho_rest = float(attrs["actin_radius"] + attrs["am_surface_gap"])
        else:
            rho_rest = float(attrs["am_preferred_centerline"])
        slack = float(attrs.get("am_tether_slack", 0.0))
        # Older campaigns (2026-07) predate these attrs, so the strain-energy
        # cross-check cannot run there.  Skip it explicitly and record the fact
        # in the artifact rather than silently trusting the reconstruction.
        can_validate = ("am_head_transverse_k" in attrs and "am_motor_k" in attrs)
        k_rho = float(attrs.get("am_head_transverse_k", 0.0))
        if not azimuth_on:
            # The aggregate pair tether is folded into each head at
            # am_tether_k / am_tether_saturation_heads.
            k_rho += float(attrs.get("am_tether_k", 0.0)) / max(
                1.0, float(attrs.get("am_tether_saturation_heads", 1)))
            Rm = 0.0          # anchor sits on the myosin axis
        k_ax_weak = float(attrs.get("am_motor_k", 0.0)) * 0.1

        sl = slice(0, frames, max(1, int(stride)))
        step = h["state/current_step"][sl, 0].astype(float)
        ac = h["actin/center"][sl]; ad = unit(h["actin/direction"][sl])
        mc = h["myosin/center"][sl]; md = unit(h["myosin/direction"][sl])
        e1r, e2r = read_material(h, sl, mc.shape[1])
        hs = h["state/am_head_state_current"][sl]
        ha = h["state/am_head_bound_actin_current"][sl].astype(int)
        haz = h["state/am_head_binding_azimuth_current"][sl]
        htg = h["state/am_head_target_axial_current"][sl]
        # stroke_progress and strain_energy are used ONLY by the validation
        # cross-check below, never by the drawn geometry, so a campaign that
        # predates them still renders -- it just renders unvalidated.
        have_val = ("state/am_head_stroke_progress_current" in h
                    and "state/am_head_strain_energy_current" in h)
        hst = (h["state/am_head_stroke_progress_current"][sl] if have_val
               else np.zeros_like(htg))
        hen = (h["state/am_head_strain_energy_current"][sl] if have_val
               else np.zeros_like(htg))
        # ACTN2 links, per the pipeline's local-link contract
        off = h["catch_bond/local_link_frame_offsets"][:]
        ident = h["catch_bond/local_link_identity"][:]
        coord = h["catch_bond/local_link_coordinates"][:]
        frame_idx = np.arange(0, frames, max(1, int(stride)))

    nA, nM = ac.shape[1], mc.shape[1]
    frames = len(step)
    spokes = int(float(attrs.get("am_spokes_per_crown", 3)))
    anchor_axial = anchor_axial_table(hph, Lm, spokes)
    assert hs.shape[2] == 2 * hph, f"{path}: head axis {hs.shape[2]} != 2*{hph}"

    e1 = unit(e1r - md * np.einsum("fmk,fmk->fm", e1r, md)[..., None])
    e2 = unit(np.cross(md, e1))
    e2 = np.where(np.einsum("fmk,fmk->fm", e2, e2r)[..., None] < 0.0, -e2, e2)

    # ---- filament segments (lab frame, PBC-clipped) ----------------------
    half_a, half_m = 0.5 * La, 0.5 * Lm
    actin_seg = np.concatenate([ac - half_a * ad, ac + half_a * ad], axis=2)
    myosin_seg = np.concatenate([mc - half_m * md, mc + half_m * md], axis=2)
    barbed = ac - half_a * ad          # barbed end is the MINUS tip

    # ---- heads -----------------------------------------------------------
    bound = (hs == WEAK) | (hs == STRONG)
    fi, mi, hi = np.nonzero(bound)
    radial = Rm * (np.cos(haz[fi, mi, hi])[:, None] * e1[fi, mi]
                   + np.sin(haz[fi, mi, hi])[:, None] * e2[fi, mi])
    m_att = mc[fi, mi] + anchor_axial[hi][:, None] * md[fi, mi] + radial
    aidx = ha[fi, mi, hi]
    d_att = min_image(ac[fi, aidx] - m_att, box, periodic)
    c_coord = np.clip(htg[fi, mi, hi], -half_a, half_a)
    conn = min_image(d_att + c_coord[:, None] * ad[fi, aidx], box, periodic)
    a_att = m_att + conn

    z = np.einsum("kj,kj->k", conn, ad[fi, aidx])
    rho = np.linalg.norm(conn - z[:, None] * ad[fi, aidx], axis=1)

    # ---- validate the reconstruction against the C++ energy --------------
    loaded = soft_threshold(rho - rho_rest, slack)
    if radial_cap > 0.0 and k_rho > 0.0:
        # saturating_radial: energy = f_max*scale*logcosh(ext/scale),
        # scale = f_max/k.  log cosh is evaluated in its stable form, the same
        # way stable_log_cosh does it in interaction.cpp.
        scale = radial_cap / k_rho
        xx = np.abs(loaded / scale)
        logcosh = np.where(xx > 20.0, xx - np.log(2.0),
                           np.log(np.cosh(np.minimum(xx, 20.0))))
        radial_e = radial_cap * scale * logcosh
    else:
        radial_e = 0.5 * k_rho * loaded ** 2
    is_weak = hs[fi, mi, hi] == WEAK
    mine = np.where(is_weak,
                    radial_e + 0.5 * k_ax_weak * z ** 2,
                    radial_e + 0.5 * k_ax_weak
                    * np.maximum(0.0, z - hst[fi, mi, hi]) ** 2)
    ref = hen[fi, mi, hi]
    if len(mine) and can_validate and have_val:
        err = np.abs(mine - ref) / np.maximum(1e-9, np.abs(ref))
        assert err.max() < 1e-6, (
            f"{path}: head reconstruction disagrees with "
            f"state/am_head_strain_energy_current, max rel err {err.max():.3e}")
        print(f"  {os.path.basename(path)}: head geometry VALIDATED against "
              f"am_head_strain_energy_current, max rel err {err.max():.2e}")
    elif len(mine):
        print(f"  {os.path.basename(path)}: WARNING head geometry NOT validated "
              f"(this campaign predates am_head_strain_energy_current / "
              f"am_head_transverse_k)")

    arm_state = hs[fi, mi, hi].astype(np.int32)
    arms = np.concatenate([m_att, a_att], axis=1)

    # unbound heads: anchor points on the myosin surface
    unb = (hs == DETACHED) | (hs == OFF)
    ufi, umi, uhi = np.nonzero(unb)
    uradial = Rm * (np.cos(haz[ufi, umi, uhi])[:, None] * e1[ufi, umi]
                    + np.sin(haz[ufi, umi, uhi])[:, None] * e2[ufi, umi])
    u_anchor = (mc[ufi, umi] + anchor_axial[uhi][:, None] * md[ufi, umi] + uradial)
    u_state = hs[ufi, umi, uhi].astype(np.int32)

    # ---- ACTN2 links -----------------------------------------------------
    aa_rows = []
    for fpos, fr in enumerate(frame_idx):
        row = off[fr]
        o, cnt = int(row[1]), int(row[2])
        if cnt <= 0:
            continue
        ids = ident[o:o + cnt]; cs = coord[o:o + cnt]
        ai = ids[:, 2].astype(int); aj = ids[:, 3].astype(int)
        stt = ids[:, 4].astype(int)
        p0 = ac[fpos, ai] + cs[:, 0][:, None] * ad[fpos, ai]
        p1raw = ac[fpos, aj] + cs[:, 1][:, None] * ad[fpos, aj]
        p1 = p0 + min_image(p1raw - p0, box, periodic)
        for k in range(cnt):
            aa_rows.append((fpos, stt[k], *p0[k], *p1[k]))
    aa = np.array(aa_rows, float) if aa_rows else np.zeros((0, 8))

    # ---- per-frame layer packing ----------------------------------------
    def pack(frame_of, seg, mask=None):
        idx = np.argsort(frame_of, kind="stable")
        fo = frame_of[idx]; sg = seg[idx]
        if mask is not None:
            keep = mask[idx]; fo, sg = fo[keep], sg[keep]
        out, offs = [], np.zeros(frames + 1, np.int32)
        counts = np.bincount(fo, minlength=frames)
        offs[1:] = np.cumsum(counts)
        return sg.reshape(-1), offs

    layers = {}
    for key, st in (("armWeak", WEAK), ("armStrong", STRONG)):
        m = arm_state == st
        v, o = pack(fi[m], arms[m])
        layers[key] = dict(values=encode(v), offsets=encode(o, np.int32), kind="segments")
    for key, st in (("anchorDet", DETACHED), ("anchorOff", OFF)):
        m = u_state == st
        v, o = pack(ufi[m], u_anchor[m])
        layers[key] = dict(values=encode(v), offsets=encode(o, np.int32), kind="points")
    for key, st in (("aa1", 1), ("aa2", 2)):
        if len(aa):
            m = aa[:, 1].astype(int) == st
            v, o = pack(aa[m, 0].astype(int), aa[m, 2:])
        else:
            v, o = np.zeros(0), np.zeros(frames + 1, np.int32)
        layers[key] = dict(values=encode(v), offsets=encode(o, np.int32), kind="segments")
    # filaments: fixed count per frame, so no offset table needed
    layers["actin"] = dict(values=encode(actin_seg.reshape(-1)), perFrame=nA, kind="segments")
    layers["myosin"] = dict(values=encode(myosin_seg.reshape(-1)), perFrame=nM, kind="segments")
    layers["barbed"] = dict(values=encode(barbed.reshape(-1)), perFrame=nA, kind="points")

    # ---- scalars ---------------------------------------------------------
    n_off = (hs == OFF).sum(axis=(1, 2)).astype(float)
    n_det = (hs == DETACHED).sum(axis=(1, 2)).astype(float)
    n_weak = (hs == WEAK).sum(axis=(1, 2)).astype(float)
    n_str = (hs == STRONG).sum(axis=(1, 2)).astype(float)
    n_aa1 = np.zeros(frames); n_aa2 = np.zeros(frames)
    if len(aa):
        for fpos, st in zip(aa[:, 0].astype(int), aa[:, 1].astype(int)):
            if st == 1: n_aa1[fpos] += 1
            elif st == 2: n_aa2[fpos] += 1
    m_order = band_order(mc[:, :, 0], box[0])

    # theta: per-frame mean over actins of the angle between the actin axis and
    # the axis of its NEAREST myosin (min-imaged centre distance).  The nearest
    # myosin, not myosin 0 and not the x axis: in this rig the thick filaments
    # are free, myosin 0 tumbles up to 44 deg off x and the two myosins differ
    # by up to 52 deg, so those references disagree with this one by up to
    # 35 deg.  Unsigned (|u_a . u_m|), so the range is [0, 90] deg and an
    # antiparallel actin reads 0 like a parallel one.
    d_cm = min_image(ac[:, :, None, :] - mc[:, None, :, :], box, periodic)
    nearest = np.linalg.norm(d_cm, axis=3).argmin(axis=2)            # (frames, nA)
    cos_am = np.abs(np.einsum('fak,fak->fa', ad, np.take_along_axis(
        md, nearest[:, :, None], axis=1)))
    theta = np.degrees(np.arccos(np.clip(cos_am, 0.0, 1.0))).mean(axis=1)

    scalars = np.stack([n_off, n_det, n_weak, n_str, n_aa1, n_aa2, m_order, step,
                        theta], axis=1)
    assert scalars.shape[1] == NSCAL

    def bounding_box(frame_mask):
        a = actin_seg[frame_mask]; m = myosin_seg[frame_mask]
        pts = np.concatenate([a.reshape(-1, 6)[:, :3], a.reshape(-1, 6)[:, 3:],
                              m.reshape(-1, 6)[:, :3], m.reshape(-1, 6)[:, 3:]], axis=0)
        lo = pts.min(axis=0); hi = pts.max(axis=0)
        pad = 0.06 * np.maximum(hi - lo, 0.1)
        return [(lo - pad).tolist(), (hi + pad).tolist()]

    bounds_all = bounding_box(np.ones(frames, bool))
    # The template's "frame the engaged episodes" control needs its own box.
    # For this rig engagement is rare, so fall back to the full-trajectory box
    # rather than handing the camera a degenerate one-frame extent.
    engaged = ((hs == WEAK) | (hs == STRONG)).any(axis=(1, 2))
    bounds_engaged = bounding_box(engaged) if engaged.sum() > 2 else bounds_all
    name = os.path.basename(path)
    return dict(
        key=name, label=label or name.replace(".h5", ""),
        frames=frames, nActin=nA, nMyosin=nM, headsPerHalf=hph,
        frameMs=float(step[1] - step[0]) * dt * 1e3 if frames > 1 else 0.0,
        box=box.tolist(), bounds=bounds_all, boundsEngaged=bounds_engaged,
        layers=layers, scalars=encode(scalars.reshape(-1)),
        nHeadTotal=int(2 * hph * nM), validated=bool(can_validate and have_val))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("runs", nargs="+")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--box", type=float, nargs=3, default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--template", required=True)
    p.add_argument("--title", default="Minimal wrap-around catch bond")
    a = p.parse_args()
    runs = [build_run(r, a.stride, a.box) for r in a.runs]
    tpl = open(a.template, encoding="ascii").read()
    payload = json.dumps(runs, separators=(",", ":"))
    assert "/*__RUNS__*/" in tpl, "template lacks /*__RUNS__*/"
    out = tpl.replace("/*__RUNS__*/", payload).replace("__TITLE__", a.title)
    assert out.isascii(), "output is not ASCII-safe"
    open(a.out, "w", encoding="ascii").write(out)
    print(f"wrote {a.out} ({os.path.getsize(a.out)/1e6:.2f} MB), "
          f"{len(runs)} run(s), frames={[r['frames'] for r in runs]}")
