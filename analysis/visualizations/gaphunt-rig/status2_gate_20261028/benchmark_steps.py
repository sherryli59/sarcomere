"""Short matched continuations for timing; production inputs remain untouched."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import h5py
import numpy as np
from connector_metrics import analyze as geometry_overshoot


def digest(path):
    return hashlib.file_digest(open(path, "rb"), "sha256").hexdigest()


def main():
    assert os.environ.get("SLURM_JOB_ID")
    out = Path(sys.argv[1])
    assert not (out/"step_benchmark.json").exists(), "refusing to overwrite timing results"
    camp = Path(os.environ["SCRATCH"])/"sarcomere/gaphunt_20261023"
    baseline = Path(os.environ["SCRATCH"])/"sarcomere/connector_force_20261028/sarcomere"
    assert digest(baseline) == "c23ae96b5e6024983ee333b6613dced1f34143a669ec51a1e3f28cf39596685d"
    source = camp/"data/hph27m24_seed20260805.h5"
    template = (camp/"logs/r27zt0npg130_seed20260805.argv").read_text().splitlines()[1:]
    steps = 200 if "--quick" in sys.argv[2:] else 500
    end_step = 200000+steps
    result = dict(job=os.environ["SLURM_JOB_ID"], host=os.uname().nodename,
                  threads=int(os.environ["OMP_NUM_THREADS"]), steps=steps,
                  baseline_sha256=digest(baseline), optimized_sha256=digest(out/"sarcomere"),
                  runs=[], caveat="Short-run timings on one allocated node; includes startup and checkpoint I/O.")
    saved = {}
    sequence = ("baseline", "optimized") if "--quick" in sys.argv[2:] else ("baseline", "optimized", "optimized", "baseline")
    for number, name in enumerate(sequence):
        target = out/f"timing_{number}_{name}.h5"
        assert not target.exists()
        subprocess.run(["cp", "--reflink=auto", str(source), str(target)], check=True)
        with h5py.File(target, "r+", locking=False) as f:
            f["meta"].attrs["titin_k"] = 1e-6
            f["meta"].attrs["aa_zdisc_tip_hard_projection"] = np.int32(0)
        binary = baseline if name == "baseline" else out/"sarcomere"
        overrides = {"--filename":str(target), "--nsteps":str(end_step),
                     "--aa_zdisc_tip_hard_projection":"false", "--save_every":str(steps)}
        cmd = [str(binary)] + [arg.split("=",1)[0]+"="+overrides[arg.split("=",1)[0]]
               if arg.split("=",1)[0] in overrides else arg for arg in template]
        (out/f"timing_{number}_{name}.argv.json").write_text(json.dumps(cmd,indent=2)+"\n")
        log = out/f"timing_{number}_{name}.log"
        start = time.monotonic()
        with log.open("w") as stream:
            subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=True)
        elapsed = time.monotonic()-start
        with h5py.File(target, "r", locking=False) as f:
            frames = int(f["meta/committed_frames"][0])
            assert int(f["state/current_step"][frames-1,0]) == end_step
            state = {species+"/"+field:f[species+"/"+field][frames-1]
                     for species in ("actin","myosin") for field in ("center","direction")}
            ev = f["catch_bond/local_link_events_v3"][:]
            cols = f["meta"].attrs["aa_local_link_events_v3_columns"]
            if isinstance(cols,bytes): cols=cols.decode()
            ec = {k:i for i,k in enumerate(cols.split(","))}
            ev = ev[ev[:,ec["step"]]>200000]
            death = ev[ev[:,ec["new_status"]]==0]
            counts = {str(k):int(np.sum(death[:,ec["detachment_reason"]]==k)) for k in range(1,6)}
            if name == "optimized":
                assert int(f["resume"].attrs["aa_mature_connector_detachment_reason"]) == 5
                assert counts["5"] > 0
            normalized = ev.copy()
            normalized[normalized[:,ec["detachment_reason"]]==5,ec["detachment_reason"]] = 3
        saved[number] = (state, normalized)
        profile = [line for line in log.read_text().splitlines()
                   if line.startswith(("runtime_profile cumulative:","aa_phase_profile:"))]
        result["runs"].append(dict(index=number,kind=name,wall_s=elapsed,
                wall_ms_per_step=elapsed*1000/steps,death_counts=counts,profile=profile))
        geometry = geometry_overshoot(target,200000,1)
        assert geometry["death_counts"]["all"]["5"] == counts["5"]
        result["runs"][-1]["geometry_overshoot"] = geometry["summary"]
        target.unlink()  # timing copies are disposable; retain logs and measurements
        (out/"step_benchmark.json").write_text(json.dumps(result,indent=2)+"\n")
        print("STEP_BENCHMARK",name,number,elapsed,counts,flush=True)
    result["comparisons"] = []
    comparisons = ((0,1),) if len(sequence)==2 else ((0,1),(1,2),(0,3),(2,3))
    for left,right in comparisons:
        x,e = saved[left]; y,f = saved[right]
        result["comparisons"].append(dict(left=left,right=right,
            state_bitwise_equal=all(np.array_equal(x[k],y[k]) for k in x),
            max_pose_error=max(float(np.max(np.abs(x[k]-y[k]))) for k in x),
            events_equal_after_reason_remap=bool(np.array_equal(e,f,equal_nan=True))))
    means = {name:float(np.mean([r["wall_s"] for r in result["runs"] if r["kind"]==name]))
             for name in ("baseline","optimized")}
    result["mean_wall_s"] = means
    result["whole_run_speedup"] = means["baseline"]/means["optimized"]
    (out/"step_benchmark.json").write_text(json.dumps(result,indent=2)+"\n")
    print("STEP_BENCHMARK_DONE",result["whole_run_speedup"],result["comparisons"],flush=True)


if __name__ == "__main__": main()
