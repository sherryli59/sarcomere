"""Resume a full checkpoint and its slim cut side by side; require identical output.

    python validate_slim_checkpoint.py FULL.h5 SLIM.h5 ARGV_TEMPLATE BINARY WORKDIR [--steps N]

Both copies get the same meta edits as the axial-stop A/B runs, then run
--steps steps (saving every --save-every) on one OpenMP thread, which makes
the trajectory deterministic. Every appended frame, ACTN2 link row, event-log
row and the final /resume snapshot must match bit for bit (NaN == NaN).
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np

from slim_checkpoint import classify, same


def prepare(source, target):
    subprocess.run(["cp", "--reflink=auto", str(source), str(target)], check=True)
    with h5py.File(target, "r+", locking=False) as f:
        a = f["meta"].attrs
        a["titin_k"] = 1e-6
        a["aa_zdisc_tip_contact_requires_link"] = np.int32(1)
        a["aa_zdisc_tip_hard_projection"] = np.int32(0)


def command(binary, template, target, end_step, save_every):
    overrides = {"--filename": str(target), "--nsteps": str(end_step),
                 "--save_every": str(save_every)}
    args = Path(template).read_text().splitlines()[1:]
    keys = [a.split("=", 1)[0] for a in args]
    assert all(keys.count(k) == 1 for k in overrides) and "--resume=true" in args
    return [str(binary)] + [k + "=" + overrides[k] if k in overrides else a for k, a in zip(keys, args)]


def run(cmd, log):
    env = dict(os.environ, OMP_NUM_THREADS="1", HDF5_USE_FILE_LOCKING="FALSE")
    started = time.monotonic()
    with open(log, "x") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, env=env, check=True)
    return time.monotonic() - started


def compare(full, slim, source_frames, source_rows):
    """source_rows: original row count of each log/ragged table in the full file."""
    problems, compared = [], 0
    with h5py.File(full, "r", locking=False) as f, h5py.File(slim, "r", locking=False) as s:
        nf, ns = int(f["meta/committed_frames"][0]), int(s["meta/committed_frames"][0])
        new_frames = nf - source_frames
        assert new_frames > 0 and ns == 1 + new_frames, (nf, ns)
        fo = f["catch_bond/local_link_frame_offsets"][source_frames:nf]
        so = s["catch_bond/local_link_frame_offsets"][1:ns]
        fv = f["catch_bond/local_link_mechanics_v2_frame_offsets"][source_frames:nf]
        sv = s["catch_bond/local_link_mechanics_v2_frame_offsets"][1:ns]

        def visit(name, obj):
            nonlocal compared
            if isinstance(obj, h5py.Group):
                return
            kind = classify(name, obj, nf)
            other = s[name]
            if kind == "copy":
                if name == "meta/committed_frames":
                    return
                a, b = obj[()], other[()]
            elif kind == "frame":
                a, b = obj[source_frames:nf], other[1:ns]
                if name.endswith("local_link_frame_offsets"):
                    a, b = a[:, [0, 2]], b[:, [0, 2]]          # step, count
                elif name.endswith("mechanics_v2_frame_offsets"):
                    a, b = a[:, [0, 3]], b[:, [0, 3]]
            elif kind in ("ragged", "ragged_v2"):
                col = 1 if kind == "ragged" else 2
                fo_, so_ = (fo, so) if kind == "ragged" else (fv, sv)
                cnt = 2 if kind == "ragged" else 3
                a = np.concatenate([obj[r[col]:r[col] + r[cnt]] for r in fo_])
                b = np.concatenate([other[r[col]:r[col] + r[cnt]] for r in so_])
            else:
                a, b = obj[source_rows[name]:], other[()]
            compared += 1
            if not same(a, b):
                problems.append(f"{name} ({kind}): differs, shapes {np.shape(a)} vs {np.shape(b)}")

        f.visititems(visit)
        for key, value in f["meta"].attrs.items():
            if not same(s["meta"].attrs[key], value):
                problems.append(f"meta attr {key} differs")
        steps = f["state/current_step"][source_frames:nf, 0].tolist()
    return problems, compared, steps


def main():
    p = argparse.ArgumentParser()
    p.add_argument("full"); p.add_argument("slim"); p.add_argument("template")
    p.add_argument("binary"); p.add_argument("workdir")
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--save-every", type=int, default=100)
    a = p.parse_args()
    work = Path(a.workdir)
    work.mkdir(parents=True, exist_ok=False)
    full, slim = work / "full.h5", work / "slim.h5"
    prepare(a.full, full)
    prepare(a.slim, slim)
    with h5py.File(full, "r", locking=False) as f:
        source_frames = int(f["meta/committed_frames"][0])
        start = int(f["state/current_step"][source_frames - 1, 0])
        rows = {}
        f.visititems(lambda n, o: rows.__setitem__(n, o.shape[0]) if isinstance(o, h5py.Dataset) else None)
    end = start + a.steps
    with ThreadPoolExecutor(2) as pool:
        walls = list(pool.map(lambda pair: run(command(a.binary, a.template, pair[0], end, a.save_every), pair[1]),
                              [(full, work / "full.log"), (slim, work / "slim.log")]))
    problems, compared, steps = compare(full, slim, source_frames, rows)
    print(f"resumed {start} -> {end}, appended steps {steps}; wall full {walls[0]:.0f} s, "
          f"slim {walls[1]:.0f} s; compared {compared} datasets")
    if problems:
        raise SystemExit("VALIDATION FAILED\n" + "\n".join(problems))
    print("VALIDATION OK: slim resume is bit-identical to the full resume")


if __name__ == "__main__":
    main()
