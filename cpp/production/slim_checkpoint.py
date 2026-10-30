"""Cut a trajectory HDF5 down to its last committed frame, keeping it resumable.

    python slim_checkpoint.py SOURCE.h5 SLIM.h5

Resume (Sarcomere::load_state) reads the last committed frame of the
per-frame datasets and the /resume snapshot; it never reads the append-only
event logs (pending events live in /resume). The slim file therefore keeps:

- per-frame datasets (actin/, myosin/, state/, actin_myo bonds/distance,
  the two catch_bond frame-offset tables): the last frame only;
- frame-ragged ACTN2 tables (local_link_identity and its row-aligned
  siblings): the last frame's rows, with offset start columns rebased to 0;
- event logs: zero rows;
- /meta (attributes intact, committed_frames = 1) and /resume, copied as is.

Every dataset keeps its dtype, chunking, maxshape and fill value, so the
simulator can append to the slim file exactly as to the original. Unknown
datasets are an error, not a guess. In the resumed output the start frame is
frame 0, not frame 400: analysis code must locate it by step.
"""
import argparse
import os
from pathlib import Path

import h5py
import numpy as np

PER_FRAME_GROUPS = ("actin", "myosin", "state")
PER_FRAME = {"actin_myo/bonds", "actin_myo/distance",
             "catch_bond/local_link_frame_offsets",
             "catch_bond/local_link_mechanics_v2_frame_offsets"}
# Row tables sliced by local_link_frame_offsets (columns step, start, count).
RAGGED = {"catch_bond/local_link_identity", "catch_bond/local_link_coordinates",
          "catch_bond/local_link_mechanics", "catch_bond/local_link_arm_state",
          "catch_bond/local_link_motor_support"}
# Sliced by local_link_mechanics_v2_frame_offsets (step, identity_start, v2_start, count).
RAGGED_V2 = {"catch_bond/local_link_mechanics_v2"}
LOG_GROUPS = ("actin_myo", "catch_bond", "titin")
COPIED_GROUPS = ("meta", "resume")


def classify(name, ds, frames):
    top = name.split("/")[0]
    if top in COPIED_GROUPS:
        return "copy"
    if top in PER_FRAME_GROUPS or name in PER_FRAME:
        if ds.shape[0] != frames or ds.maxshape[0] is not None:
            raise ValueError(f"{name}: expected {frames} appendable frames, got {ds.shape}")
        return "frame"
    if name in RAGGED:
        return "ragged"
    if name in RAGGED_V2:
        return "ragged_v2"
    if top in LOG_GROUPS and ds.maxshape[0] is None:
        return "log"
    raise ValueError(f"{name}: unknown dataset {ds.shape}; extend the classification first")


def create_like(dst, name, ds, data):
    out = dst.create_dataset(
        name, shape=data.shape, maxshape=ds.maxshape, chunks=ds.chunks, dtype=ds.dtype,
        compression=ds.compression, compression_opts=ds.compression_opts,
        shuffle=ds.shuffle, fletcher32=ds.fletcher32, fillvalue=ds.fillvalue)
    if data.size:
        out[...] = data
    out.attrs.update(ds.attrs)


def slim(source, target):
    target = Path(target)
    if target.exists():
        raise FileExistsError(target)
    tmp = target.with_name(target.name + f".tmp{os.getpid()}")
    counts = {}
    with h5py.File(source, "r", locking=False) as src, h5py.File(tmp, "x") as dst:
        frames = int(src["meta/committed_frames"][0])
        last = frames - 1
        step = int(src["state/current_step"][last, 0])
        assert int(src["state/state_commit_count"].shape[0]) == frames, "state and geometry frames differ"
        assert int(src["resume/current_step"][0]) == step, "/resume is not at the last committed frame"
        offsets = src["catch_bond/local_link_frame_offsets"][last]
        assert int(offsets[0]) == step
        start, count = int(offsets[1]), int(offsets[2])
        v2 = src["catch_bond/local_link_mechanics_v2_frame_offsets"][last]
        assert int(v2[0]) == step and int(v2[1]) == start and int(v2[3]) == count
        v2_start = int(v2[2])
        dst.attrs.update(src.attrs)

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                dst.require_group(name).attrs.update(obj.attrs)
                return
            kind = classify(name, obj, frames)
            counts[kind] = counts.get(kind, 0) + 1
            if kind == "copy":
                src.copy(obj, dst.require_group(os.path.dirname(name)), name=os.path.basename(name))
            elif kind == "frame":
                data = obj[last:frames]
                if name == "catch_bond/local_link_frame_offsets":
                    data = data.copy(); data[0, 1] = 0
                elif name == "catch_bond/local_link_mechanics_v2_frame_offsets":
                    data = data.copy(); data[0, 1] = 0; data[0, 2] = 0
                create_like(dst, name, obj, data)
            elif kind == "ragged":
                assert obj.shape[0] >= start + count, name
                create_like(dst, name, obj, obj[start:start + count])
            elif kind == "ragged_v2":
                create_like(dst, name, obj, obj[v2_start:v2_start + count])
            else:
                create_like(dst, name, obj, np.empty((0,) + obj.shape[1:], obj.dtype))

        src.visititems(visit)
        dst["meta/committed_frames"][0] = 1
        meta = dst["meta"].attrs
        meta["slim_checkpoint_source"] = np.bytes_(str(Path(source).resolve()))
        meta["slim_checkpoint_source_frame"] = np.int32(last)
        meta["slim_checkpoint_step"] = np.int32(step)
    os.replace(tmp, target)
    return dict(step=step, source_frames=frames, links=count, datasets=counts,
                source_bytes=os.path.getsize(source), slim_bytes=os.path.getsize(target))


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    nan_ok = a.dtype.kind in "fc" and b.dtype.kind in "fc"
    return a.shape == b.shape and np.array_equal(a, b, equal_nan=nan_ok)


def check(source, target):
    """Every kept value equals the source's last frame; /resume and /meta unchanged."""
    with h5py.File(source, "r", locking=False) as src, h5py.File(target, "r", locking=False) as dst:
        frames = int(src["meta/committed_frames"][0])
        last = frames - 1
        o = src["catch_bond/local_link_frame_offsets"][last]
        v2 = src["catch_bond/local_link_mechanics_v2_frame_offsets"][last]
        problems = []

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                return
            kind = classify(name, obj, frames)
            got = dst[name]
            if (got.dtype, got.maxshape, got.chunks) != (obj.dtype, obj.maxshape, obj.chunks):
                problems.append(f"{name}: layout differs")
            if kind == "copy":
                want = obj[()]
                if name == "meta/committed_frames":
                    want = np.array([1], obj.dtype)
            elif kind == "frame":
                want = obj[last:frames].copy()
                if name.endswith("local_link_frame_offsets"):
                    want[0, 1] = 0
                elif name.endswith("mechanics_v2_frame_offsets"):
                    want[0, 1:3] = 0
            elif kind == "ragged":
                want = obj[int(o[1]):int(o[1]) + int(o[2])]
            elif kind == "ragged_v2":
                want = obj[int(v2[2]):int(v2[2]) + int(v2[3])]
            else:
                want = np.empty((0,) + obj.shape[1:], obj.dtype)
            if not same(got[()], want):
                problems.append(f"{name}: values differ")
            if set(got.attrs) != set(obj.attrs):
                problems.append(f"{name}: attributes differ")

        src.visititems(visit)
        for key, value in src["meta"].attrs.items():
            if not same(dst["meta"].attrs[key], value):
                problems.append(f"meta attr {key} differs")
        return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("source")
    parser.add_argument("target")
    parser.add_argument("--no-check", action="store_true")
    args = parser.parse_args()
    report = slim(args.source, args.target)
    print(f"step {report['step']}: kept frame {report['source_frames'] - 1} of "
          f"{report['source_frames']}, {report['links']} ACTN2 links; "
          f"{report['source_bytes'] / 1e9:.2f} GB -> {report['slim_bytes'] / 1e6:.1f} MB; "
          f"datasets {report['datasets']}", flush=True)
    if not args.no_check:
        problems = check(args.source, args.target)
        if problems:
            raise SystemExit("CHECK FAILED\n" + "\n".join(problems))
        print("CHECK OK: slim file matches the source's last frame", flush=True)


if __name__ == "__main__":
    main()
