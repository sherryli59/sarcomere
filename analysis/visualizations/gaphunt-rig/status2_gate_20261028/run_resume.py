"""Stage and run the matched status-2 exclusion resume, inside Slurm only."""
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import h5py
import numpy as np

EXPECTED_BINARY = "dfa63b085df368246e93f8c27f121f4884b6f1e08fa0fc3f77379b65bb84c47b"
ARM = "r27zt0npg130"
SEED = 20260805


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this production preparation and simulation in Slurm")
    camp = Path(os.environ["SCRATCH"]) / "sarcomere/gaphunt_20261023"
    validated = Path(os.environ["SCRATCH"]) / "sarcomere/status2_gate_20261028.983jQo/sarcomere"
    binary = camp / "bin/sarcomere_status2_gate_20261028"
    result_dir = camp / "torque_out/status2_gate_20261028"
    result_dir.mkdir(parents=True, exist_ok=True)
    if binary.exists():
        assert digest(binary) == EXPECTED_BINARY, "pinned binary differs"
    else:
        assert digest(validated) == EXPECTED_BINARY, "validated binary differs"
        shutil.copy2(validated, binary)
    original_argv = (camp / f"logs/hph27m24_seed{SEED}.argv").read_text().splitlines()
    source = camp / f"data/hph27m24_seed{SEED}.h5"
    outputs = [camp / f"data/{ARM}{suffix}_seed{SEED}.h5" for suffix in ("_smoke", "")]
    for target in outputs:
        if target.exists():
            raise FileExistsError(f"refusing duplicate or overwrite: {target}")
    source_hash = digest(source)
    with h5py.File(source, "r", locking=False) as f:
        n = int(f["meta/committed_frames"][0])
        assert int(f["state/current_step"][n - 1, 0]) == 200000
        assert float(f["meta"].attrs["dt"]) == 1e-5
        assert float(f["meta"].attrs["k_mm"]) == 0.0
        assert float(f["meta"].attrs["box_x"]) == 6.0
        assert int(f["meta"].attrs["am_heads_per_half"]) == 27
        offsets = f["catch_bond/local_link_frame_offsets"][400]
        start, count = int(offsets[1]), int(offsets[2])
        source_mature = int((f["catch_bond/local_link_identity"][start:start + count, 4] == 2).sum())
        initial_pose = f["actin/center"][400]
    manifest = dict(
        arm=ARM, seed=SEED, source=str(source), source_sha256=source_hash,
        binary=str(binary), binary_sha256=EXPECTED_BINARY,
        source_argv=str(camp / f"logs/hph27m24_seed{SEED}.argv"),
        window_steps=[200000, 240000], window_s=[2.0, 2.4],
        projection=False, lifetime_coeff=0.4, titin_k=1e-6, k_mm=0.0,
        status2_exclusion_radius_um=0.130,
        status2_exclusion_geometry="full_connector_finite_cylinder_axial_overlap_v1",
        inherited_mature_links=source_mature, slurm_job_id=os.environ["SLURM_JOB_ID"],
        host=os.uname().nodename,
        cpu_info=subprocess.check_output(["lscpu"], text=True),
        original_argv=original_argv,
        preparation_script_sha256=digest(Path(__file__)),
        runs=[],
    )
    manifest_path = result_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    for target, end_step in zip(outputs, (200500, 240000)):
        shutil.copy2(source, target)
        with h5py.File(target, "r+", locking=False) as f:
            a = f["meta"].attrs
            a["titin_k"] = 1e-6
            a["aa_zdisc_tip_contact_requires_link"] = np.int32(1)
            a["aa_zdisc_tip_hard_projection"] = np.int32(0)
            a["aa_status2_connector_exclusion_radius"] = 0.130
            a["aa_status2_connector_exclusion_geometry"] = np.bytes_(
                "full_connector_finite_cylinder_axial_overlap_v1")
            a["ab_binary_sha256"] = np.bytes_(EXPECTED_BINARY)
            a["ab_source_sha256"] = np.bytes_(source_hash)
            a["resume_ab_note"] = np.bytes_(
                "Frames 0-400 are copied from hph27m24 (legacy maturation, titin_k=100, "
                "projection on). At step 200000: titin_k=1e-6, k_mm=0, "
                "requires_link=true, projection=false, lifetime_coeff=0.4; new status 1->2 "
                "requires full connector outside every 130-nm finite myosin cylinder. "
                "Inherited status-2 links keep legacy retention. See run_manifest.json.")
        overrides = {
            "--filename": str(target), "--nsteps": str(end_step), "--resume": "true",
            "--titin_k": "1e-6", "--aa_zdisc_tip_contact_requires_link": "true",
            "--aa_zdisc_tip_hard_projection": "false", "--lifetime_coeff": "0.4",
        }
        cmd = [str(binary)]
        remaining = dict(overrides)
        for arg in original_argv[1:]:
            key = arg.split("=", 1)[0]
            if key in overrides:
                cmd.append(key + "=" + overrides[key])
                remaining.pop(key, None)
            else:
                cmd.append(arg)
        cmd.extend(k + "=" + v for k, v in remaining.items())
        argv_path = camp / ("logs/" + target.stem + ".argv")
        argv_path.write_text("\n".join(cmd) + "\n")
        print("".join(difflib.unified_diff(original_argv, cmd, lineterm="\n")), flush=True)
        started = time.monotonic()
        print(f"RUN {target.stem} to {end_step} host={os.uname().nodename}", flush=True)
        # Each stage has its own log so runtime profiles are unambiguous.
        log_path = camp / ("logs/" + target.stem + ".simulation.out")
        with open(log_path, "x") as log:
            completed = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        wall = time.monotonic() - started
        entry = dict(path=str(target), argv=str(argv_path), log=str(log_path),
                     end_step=end_step, wall_s=wall, exit_code=completed.returncode)
        manifest["runs"].append(entry)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"exit={completed.returncode} {target.stem} wall_s={wall:.3f}", flush=True)
        completed.check_returncode()
        with h5py.File(target, "r", locking=False) as f:
            n = int(f["meta/committed_frames"][0])
            assert int(f["state/current_step"][n - 1, 0]) == end_step
            assert np.array_equal(f["actin/center"][400], initial_pose)
            assert f["resume"].attrs["aa_status2_connector_exclusion_radius"] == 0.130
        if end_step == 200500:
            print("SMOKE_OK: load, advance, save and checkpoint provenance verified", flush=True)
    assert digest(source) == source_hash, "source changed during run"
    manifest["complete"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
