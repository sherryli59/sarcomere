"""Build an isolated logger-only replay, verify exact admission geometry in Python."""
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from audit_new_mature import blocked_connectors, unit


def main():
    assert os.environ.get("SLURM_JOB_ID")
    out = Path(sys.argv[1])
    build = Path(os.environ["SCRATCH"])/"sarcomere/status2_gate_20261028.983jQo"
    camp = Path(os.environ["SCRATCH"])/"sarcomere/gaphunt_20261023"
    compiler = Path(os.environ["GROUP_HOME"])/"sherry/mamba_root/envs/sarcomere/bin/g++"
    flags = (build/"CMakeFiles/my_library.dir/flags.make").read_text().splitlines()
    includes = shlex.split(next(s.split("=",1)[1] for s in flags if s.startswith("CXX_INCLUDES =")))
    options = shlex.split(next(s.split("=",1)[1] for s in flags if s.startswith("CXX_FLAGS =")))
    cmd = [str(compiler),*options,*includes,"-I"+str(HERE),"-c",str(out/"audit_actn2.cpp"),"-o",str(out/"audit_actn2.o")]
    print("COMPILE",shlex.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    link = shlex.split((build/"CMakeFiles/sarcomere.dir/link.txt").read_text())
    for k,token in enumerate(link):
        if token == "sarcomere": link[k] = str(out/"sarcomere_exact_audit")
        elif token.endswith(".o") or token == "libmy_library.a": link[k] = str(build/token)
        elif token.startswith("-Wl,--dependency-file="): link[k] = "-Wl,--dependency-file="+str(out/"link.d")
    link.insert(link.index("-o"),str(out/"audit_actn2.o"))
    subprocess.run(link,check=True,cwd=out)
    source = camp/"data/hph27m24_seed20260805.h5"
    target = out/"replay.h5"
    assert not target.exists(), "refusing to overwrite replay"
    shutil.copy2(source,target)
    with h5py.File(target,"r+",locking=False) as f:
        f["meta"].attrs["titin_k"] = 1e-6
    argv = (camp/"logs/r27zt0npg130_seed20260805.argv").read_text().splitlines()
    cmd = [str(out/"sarcomere_exact_audit")]
    for token in argv[1:]:
        key = token.split("=",1)[0]
        if key=="--filename": token="--filename="+str(target)
        elif key=="--nsteps": token="--nsteps=202000"
        cmd.append(token)
    env = dict(os.environ,ACTN2_GATE_AUDIT_PATH=str(out/"admission.jsonl"))
    (out/"argv.json").write_text(json.dumps(cmd,indent=2)+"\n")
    started = time.monotonic()
    with open(out/"simulation.log","w") as log:
        subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print("REPLAY_DONE",time.monotonic()-started,flush=True)
    accepted, rejected, violations = [], [], []
    examples = []
    for line in (out/"admission.jsonl").open():
        row = json.loads(line)
        if not row["admitted"]:
            rejected.append((row["step"],row["link_id"]))
            continue
        m = np.asarray(row["myosins"])
        hit,depth,rho,which,ax = blocked_connectors(
            np.array([row["first"]]),np.array([row["second"]]),
            m[:,:3],unit(m[:,3:]),np.array([6.,2.,2.]),1.5)
        accepted.append((row["step"],row["link_id"]))
        if hit[0]:
            violations.append(dict(step=row["step"],link_id=row["link_id"],
                myosin=int(which[0]),rho_nm=float(rho[0]*1000),depth_nm=float(depth[0]*1000)))
        if row["link_id"] in (164164,161458,162680,162686):
            examples.append(row)
    production = camp/"data/r27zt0npg130_seed20260805.h5"
    comparisons = []
    with h5py.File(target,"r",locking=False) as f, h5py.File(production,"r",locking=False) as g:
        for frame in (400,401,402,403,404):
            for species in ("actin","myosin"):
                for field in ("center","direction"):
                    x,y=f[species+"/"+field][frame],g[species+"/"+field][frame]
                    comparisons.append(dict(frame=frame,species=species,field=field,
                                            bitwise_equal=bool(np.array_equal(x,y)),
                                            max_abs_error=float(np.max(np.abs(x-y)))))
        ec={k:i for i,k in enumerate(g["meta"].attrs["aa_local_link_events_v3_columns"].decode().split(","))}
        def transitions(h):
            ev=h["catch_bond/local_link_events_v3"][:]
            q=(ev[:,ec["step"]]>200000)&(ev[:,ec["step"]]<=202000)&(ev[:,ec["old_status"]]==1)&(ev[:,ec["new_status"]]==2)
            return set((int(r[ec["step"]]),int(r[ec["link_id"]])) for r in ev[q])
        replay_events,production_events=transitions(f),transitions(g)
    result=dict(
        window_steps=[200000,202000],job_id=os.environ["SLURM_JOB_ID"],
        replay_threads=int(os.environ["OMP_NUM_THREADS"]),production_threads=16,
        accepted_proposals_independently_checked=len(accepted),
        rejected_proposals=len(rejected),
        accepted_unique_links=len(set(lid for _,lid in accepted)),
        exact_geometry_violations=violations,
        replay_maturation_events=len(replay_events),
        production_maturation_events=len(production_events),
        maturation_event_sets_equal=replay_events==production_events,
        all_actual_maturation_events_have_admission_record=production_events.issubset(set(accepted)),
        saved_pose_comparisons=comparisons,
        source_translation_unit_sha256=hashlib.sha256((out/"audit_actn2.cpp").read_bytes()).hexdigest(),
        output=str(out),example_link_ids=[r["link_id"] for r in examples],
    )
    (out/"exact_summary.json").write_text(json.dumps(result,indent=2)+"\n")
    (out/"example_exact_geometry.json").write_text(json.dumps(examples,indent=2)+"\n")
    print("EXACT_AUDIT",json.dumps(result),flush=True)
    assert not violations, "the gate admitted excluded geometry"
    assert replay_events==production_events, "replay changed the maturation decisions"
    assert all(c["max_abs_error"] < 1e-9 for c in comparisons), "replay diverged from production by more than 1e-9 um"


if __name__=="__main__": main()
