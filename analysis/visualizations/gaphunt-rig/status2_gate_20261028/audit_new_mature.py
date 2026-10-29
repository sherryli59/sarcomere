"""Independent saved-pose cylinder audit, following post-resume 1->2 link IDs."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

CAMP = Path(os.environ["SCRATCH"]) / "sarcomere/gaphunt_20261023"
OUT = CAMP / "torque_out/status2_gate_20261028/maturation_audit"
R = .130


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def wrap(x, box):
    # No half-box ties occur in the measured configurations.
    return x - box * np.round(x / box)


def blocked_connectors(first, second, mc, mu, box, length):
    """Independently clip segments to each myosin slab and minimize radius.
    For this production geometry, the cylinder plus half connector length
    fits inside half every box dimension: only the nearest midpoint image
    can intersect. Assert that condition rather than silently assuming it.
    """
    d = wrap(second - first, box)
    reach = np.hypot(.5 * length, R) + .5*np.linalg.norm(d, axis=1)
    assert np.all(reach < .5*np.min(box)), "requires explicit additional images"
    midpoint = first + .5*d
    n = len(first)
    hit = np.zeros(n, bool)
    best_depth = np.full(n, -np.inf)
    best_radius = np.full(n, np.nan)
    best_m = np.full(n, -1, int)
    best_axial = np.full(n, np.nan)
    for m, (center, axis) in enumerate(zip(mc, mu)):
        relative_mid = wrap(midpoint - center, box)
        near = np.flatnonzero(np.sum(relative_mid**2, axis=1) <= reach**2)
        if not len(near): continue
        p = relative_mid[near] - .5*d[near]
        delta = d[near]
        z0 = np.sum(p * axis, axis=1)
        dz = np.sum(delta * axis, axis=1)
        parallel = dz == 0
        divisor = np.where(parallel, 1, dz)
        a = (-.5*length-z0)/divisor
        b = (.5*length-z0)/divisor
        lo = np.maximum(0, np.minimum(a,b))
        hi = np.minimum(1, np.maximum(a,b))
        lo[parallel] = 0
        hi[parallel] = 1
        eligible = (lo <= hi) & (~parallel | (np.abs(z0) <= .5*length))
        rp = p - z0[:,None]*axis
        rd = delta - dz[:,None]*axis
        vv = np.sum(rd**2, axis=1)
        t = np.clip(-np.sum(rp*rd, axis=1)/np.where(vv>0,vv,1), lo, hi)
        rho = np.linalg.norm(rp+t[:,None]*rd, axis=1)
        eligible &= rho <= R
        depth = .5*length - np.abs(z0+t*dz)
        improve = eligible & (depth > best_depth[near])
        idx = near[improve]
        hit[near[eligible]] = True
        best_depth[idx] = depth[improve]
        best_radius[idx] = rho[improve]
        best_m[idx] = m
        best_axial[idx] = (z0+t*dz)[improve]
    return hit, best_depth, best_radius, best_m, best_axial


def analyze(arm):
    path = CAMP / f"data/{arm}_seed20260805.h5"
    with h5py.File(path, "r", locking=False) as f:
        a = f["meta"].attrs
        box = np.array([a["box_x"],a["box_y"],a["box_z"]])
        assert all(a["periodic_"+k] for k in "xyz")
        dt = float(a["dt"])
        nframes = int(f["meta/committed_frames"][0])
        steps = f["state/current_step"][:nframes,0]
        ec = {s:i for i,s in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
        ev = f["catch_bond/local_link_events_v3"][:]
        mat = ev[(ev[:,ec["old_status"]]==1)&(ev[:,ec["new_status"]]==2)&(ev[:,ec["step"]]>200000)]
        mat_step = {int(row[ec["link_id"]]):int(row[ec["step"]]) for row in mat}
        assert len(mat_step)==len(mat), "multiple maturations need phase-specific tracking"
        born = set(ev[(ev[:,ec["old_status"]]==0)&(ev[:,ec["new_status"]]==1)&(ev[:,ec["step"]]>200000),ec["link_id"]].astype(int))
        off = f["catch_bond/local_link_frame_offsets"]
        ident = f["catch_bond/local_link_identity"]
        coords = f["catch_bond/local_link_coordinates"]
        row = off[400]
        source_rows = ident[int(row[1]):int(row[1]+row[2])]
        inherited = set(source_rows[source_rows[:,4]==2,1].astype(int))
        first_seen = {}
        first_clear_then_enter = {}
        first_blocked = []
        later_entries = []
        frames = []
        totals = {k:dict(samples=0,blocked=0,deep=0) for k in ("status1","inherited2","new2")}
        for frame in range(401,nframes):
            start,count = map(int,off[frame,1:3])
            ids = ident[start:start+count]
            sites = coords[start:start+count]
            ca = f["actin/center"][frame]
            ua = unit(f["actin/direction"][frame])
            mc = f["myosin/center"][frame]
            mu = unit(f["myosin/direction"][frame])
            p0 = ca[ids[:,2]] + sites[:,0,None]*ua[ids[:,2]]
            p1 = ca[ids[:,3]] + sites[:,1,None]*ua[ids[:,3]]
            hit,depth,rho,m,ax = blocked_connectors(p0,p1,mc,mu,box,float(a["myosin_length"]))
            mature = ids[:,4]==2
            new = mature & np.isin(ids[:,1],list(mat_step))
            old = mature & np.isin(ids[:,1],list(inherited))
            groups = {"status1":ids[:,4]==1,"inherited2":old,"new2":new}
            record = dict(frame=frame,step=int(steps[frame]))
            for name,q in groups.items():
                deep = q & hit & (depth>.050) & (rho<.120)
                record[name] = dict(n=int(q.sum()),blocked=int((q&hit).sum()),deep=int(deep.sum()))
                for key,value in (("samples",int(q.sum())),("blocked",int((q&hit).sum())),("deep",int(deep.sum()))):
                    totals[name][key] += value
            for k in np.flatnonzero(new):
                lid = int(ids[k,1])
                rec = dict(link_id=lid,actins=ids[k,2:4].astype(int).tolist(),
                           matured_step=mat_step[lid],born_after_resume=lid in born,
                           sample_step=int(steps[frame]),frame=frame,
                           delay_ms=(int(steps[frame])-mat_step[lid])*dt*1000,
                           blocked=bool(hit[k]),myosin=int(m[k]),
                           radial_nm=float(rho[k]*1000) if hit[k] else None,
                           axial_depth_nm=float(depth[k]*1000) if hit[k] else None)
                if lid not in first_seen:
                    first_seen[lid] = rec
                    if hit[k]: first_blocked.append(rec)
                elif hit[k] and not first_seen[lid]["blocked"] and lid not in first_clear_then_enter:
                    first_clear_then_enter[lid]=rec
                    later_entries.append(rec)
            frames.append(record)
            if frame%20==0: print(arm,frame,record,flush=True)
        first_blocked.sort(key=lambda r:(r["sample_step"],-r["axial_depth_nm"]))
        first_rows = list(first_seen.values())
        event_keys = {int(row[ec["link_id"]]):row for row in mat}
        examples = []
        for rec in first_blocked[:20]:
            row = event_keys[rec["link_id"]]
            examples.append(dict(rec,event_sites_um=[float(row[ec["s_i"]]),float(row[ec["s_j"]])]))
        deaths = ev[(ev[:,ec["step"]]>200000)&(ev[:,ec["new_status"]]==0)]
        new_deaths = deaths[np.isin(deaths[:,ec["link_id"]].astype(int),list(mat_step))]
        result = dict(
            arm=arm,input=str(path),radius_nm=130,window_steps=[200000,240000],
            post_resume_maturation_events=len(mat),observed_in_saved_frames=len(first_seen),
            died_before_first_saved_pose=len(mat)-len(first_seen),
            first_saved_pose_blocked=len(first_blocked),
            first_saved_pose_deep_and_radial_below_120nm=sum(
                r["axial_depth_nm"]>50 and r["radial_nm"]<120 for r in first_blocked),
            first_clear_then_later_entered=len(first_clear_then_enter),
            first_sample_delay_ms_percentiles=np.percentile([r["delay_ms"] for r in first_rows],[0,50,90,100]).tolist(),
            zero_delay_samples=[r for r in first_rows if r["delay_ms"]==0],
            occupancy_link_frames=totals,per_frame=frames,first_blocked_examples=examples,
            first_clear_then_enter_examples=later_entries[:20],
            new_mature_death_reasons={str(k):int((new_deaths[:,ec["detachment_reason"]]==k).sum()) for k in (1,2,3,4,5)},
            caveat="Stored poses follow integration and are up to 5 ms after the actual 1->2 decision; first-pose overlap is not proof of a maturation-time violation."
        )
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/(arm+".json")).write_text(json.dumps(result,indent=2)+"\n")
        print("AUDIT_SUMMARY",json.dumps({k:v for k,v in result.items() if k not in (
            "per_frame","first_blocked_examples","zero_delay_samples","first_clear_then_enter_examples")}),flush=True)
        print("EARLIEST_BLOCKED",json.dumps(examples[:4]),flush=True)
        return result


if __name__=="__main__":
    if not os.environ.get("SLURM_JOB_ID"): raise RuntimeError("Run this trajectory audit inside Slurm")
    for arm in ("r27zt0np","r27zt0npg130"): analyze(arm)
