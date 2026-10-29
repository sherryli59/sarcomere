"""Requested patch forces and recorded force budgets for the two entry clips.
The patch kernel is reconstructed from stored poses; recorded contact fields
refer to force evaluation before the final integrator movement within a step.
No selected-link ID is stored, so a matching gap/radius is evidence, not an
exact identity certificate.
"""
import json
import os
from pathlib import Path
import h5py
import numpy as np
from audit_new_mature import unit, wrap
from render_entry_movies import geometry

ROOT = Path(os.environ["SCRATCH"])/"sarcomere/gaphunt_20261023"
OUT = ROOT/"torque_out/status2_gate_20261028"
CASES = OUT/"entry_movies_v2/manifest.json"


def main():
    assert os.environ.get("SLURM_JOB_ID")
    cases=json.loads(CASES.read_text())["cases"]
    source=ROOT/"data/r27zt0npg130_seed20260805.h5"
    result=dict(input=str(source),cases=[],
        caveat="Kernel values reconstructed at saved pose; recorded forces evaluated just before that pose. Selected link ID is not recorded.")
    with h5py.File(source,"r",locking=False) as f:
        a=f["meta"].attrs
        box=np.array([a["box_x"],a["box_y"],a["box_z"]])
        lm=float(a["myosin_length"])
        flat=float(a["aa_zdisc_tip_transverse_flat"])
        cutoff=float(a["aa_zdisc_tip_transverse_cutoff"])
        onset=float(a["aa_zdisc_tip_onset"])
        stiffness=float(a["aa_zdisc_tip_stiffness"])
        fmax=float(a["aa_zdisc_tip_force_max"])
        result["parameters"]={k:float(a[k]) for k in (
            "aa_zdisc_tip_transverse_flat","aa_zdisc_tip_transverse_cutoff",
            "aa_zdisc_tip_force_max","max_actin_force","max_myosin_force",
            "aa_actn2_contact_radius","myosin_radius","aa_myosin_contact_shell_width")}
        for case in cases:
            lid,m,sign=case["lid"],case["m"],case["end_sign"]
            side="plus" if sign>0 else "minus"
            other="minus" if sign>0 else "plus"
            rows=[]
            for frame in range(case["start"],case["end"]+1):
                off=f["catch_bond/local_link_frame_offsets"][frame]
                start,count=map(int,off[1:3])
                ids=f["catch_bond/local_link_identity"][start:start+count]
                which=np.flatnonzero(ids[:,1]==lid)
                if not len(which):continue
                k=int(which[0]);i,j=map(int,ids[k,2:4])
                si,sj=f["catch_bond/local_link_coordinates"][start+k]
                ca=f["actin/center"][frame];ua=unit(f["actin/direction"][frame])
                cm=f["myosin/center"][frame,m];u=unit(f["myosin/direction"][frame,m])
                p0=ca[i]+si*ua[i];p1=p0+wrap(ca[j]+sj*ua[j]-p0,box)
                patch=.5*(p0+p1);outward=sign*u;tip=cm+.5*lm*outward
                rel=wrap(patch-tip,box)
                gap=float(rel@outward)
                transverse=rel-gap*outward
                rho=float(np.linalg.norm(transverse))
                eligible_half=bool(sign*float(wrap(patch-cm,box)@u)>=0)
                gate=0.
                if rho<cutoff:
                    q=np.clip((rho-flat)/(cutoff-flat),0,1)
                    gate=float(1-3*q*q+2*q*q*q)
                raw=float(gate*fmax*np.tanh(stiffness*(onset-gap)/fmax)) if eligible_half and gap<onset else 0.
                contact=geometry(p0,p1,cm,u,box,lm)
                sf=float(f["myosin/aa_zdisc_tip_force_"+side][frame,m,0])
                of=float(f["myosin/aa_zdisc_tip_force_"+other][frame,m,0])
                sg=float(f["myosin/aa_zdisc_tip_gap_"+side][frame,m,0])
                sr=float(f["myosin/aa_zdisc_tip_rho_"+side][frame,m,0])
                mf=f["myosin/force"][frame,m]
                record=dict(frame=frame,step=int(f["state/current_step"][frame,0]),
                    age_ms=(int(f["state/current_step"][frame,0])-case["mature_step"])*.01,
                    status=int(ids[k,4]),connector_excluded=contact["blocked"],
                    connector_min_radius_nm=None if contact["rho"] is None else 1000*contact["rho"],
                    midpoint_radius_nm=1000*rho,midpoint_gap_nm=1000*gap,
                    eligible_half=eligible_half,radial_gate=gate,
                    candidate_patch_force_pN=raw,
                    recorded_selected_tip_force_pN=sf,recorded_opposite_tip_force_pN=of,
                    two_tip_force_on_myosin_outward_pN=of-sf,
                    recorded_selected_gap_nm=1000*sg,recorded_selected_radius_nm=1000*sr,
                    gap_difference_from_selected_nm=1000*(gap-sg),
                    radius_difference_from_selected_nm=1000*(rho-sr),
                    myosin_net_force_mag_pN=float(np.linalg.norm(mf)),
                    myosin_net_force_outward_pN=float(mf@outward),actins=[])
                for idx in (i,j):
                    af=f["actin/force"][frame,idx]
                    cap=f["actin/cumulative_force_cap_events"]
                    freq=float((cap[frame,idx,0]-cap[frame-1,idx,0])/500) if frame>400 else None
                    record["actins"].append(dict(id=idx,
                        patch_force_total_pN=float(f["actin/aa_zdisc_tip_force"][frame,idx,0]),
                        force_mag_pN=float(np.linalg.norm(af)),
                        force_outward_pN=float(af@outward),cap_fraction_since_last_frame=freq,
                        barbed_endstop_force_pN=float(f["actin/aa_barbed_tip_endstop_force"][frame,idx,0])))
                rows.append(record)
            result["cases"].append(dict(case=case,rows=rows))
            for r in rows:
                if r["frame"] in (case["entry"]-2,case["entry"]-1,case["entry"],case["entry"]+1,case["entry"]+2,case["poster_frame"]):
                    print("FORCE_SAMPLE",case["case"],json.dumps(r),flush=True)
    (OUT/"entry_force_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    durable=Path(__file__).resolve().parent/"results/entry_force_audit.json"
    durable.write_text(json.dumps(result,indent=2)+"\n")
    print("ENTRY_FORCE_AUDIT_DONE",flush=True)


if __name__=="__main__": main()
