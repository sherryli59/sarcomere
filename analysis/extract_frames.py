#!/usr/bin/env python3
"""Write a small HDF5 holding only the frames and datasets visualize_traj.py reads.

visualize_traj.py calls hdf5_to_dict(), which materialises the WHOLE file as
numpy arrays before rendering a single frame, and then hands that dict to every
joblib worker. On a 14.6 GB trajectory that is an immediate OOM. It only ever
reads twelve frame-indexed datasets, so a subset file is a drop-in.
"""
import argparse, os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py
import numpy as np

WANTED = ["/actin/center", "/actin/direction", "/actin/cb_status", "/actin/f_load",
          "/actin/force", "/actin/velocity",
          "/myosin/center", "/myosin/direction", "/myosin/force",
          "/myosin/velocity", "/myosin/bonds", "/actin_myo/bonds"]

ap = argparse.ArgumentParser()
ap.add_argument("--source", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--frames", type=int, default=1, help="how many, evenly spaced")
ap.add_argument("--first", type=int, default=None, help="single frame index")
a = ap.parse_args()

with h5py.File(a.source, "r", locking=False) as src:
    n_committed = int(src["meta/committed_frames"][0]) if "meta/committed_frames" in src else None
    lengths = [src[k].shape[0] for k in WANTED if k in src]
    usable = min([n for n in ([n_committed] if n_committed else []) + lengths])
    picks = ([a.first] if a.first is not None
             else list(np.unique(np.linspace(0, usable - 1, a.frames).round().astype(int))))
    print(f"source has {usable} usable frames; extracting {len(picks)}: "
          f"{picks[:4]}{' ...' if len(picks) > 4 else ''}")
    with h5py.File(a.out, "w") as dst:
        for key in WANTED:
            if key not in src:
                print(f"  {key:26s} absent, skipped")
                continue
            arr = np.stack([np.asarray(src[key][int(f)]) for f in picks])
            dst.create_dataset(key.lstrip("/"), data=arr, compression="gzip",
                               compression_opts=1)
            print(f"  {key:26s} {arr.shape}  {arr.nbytes/1048576:7.1f} MB")
        # ACTN2 tables: concatenate the per-frame slices and rebuild offsets.
        op = "catch_bond/local_link_frame_offsets"
        ip = "catch_bond/local_link_identity"
        cp = "catch_bond/local_link_coordinates"
        if all(k in src for k in (op, ip, cp)):
            rows, ident, coord, cursor = [], [], [], 0
            for f in picks:
                step, offset, count = np.asarray(src[op][int(f)], dtype=int)
                ident.append(np.asarray(src[ip][offset:offset + count]))
                coord.append(np.asarray(src[cp][offset:offset + count]))
                rows.append([step, cursor, count])
                cursor += count
            dst.create_dataset(op, data=np.asarray(rows, dtype=np.int64))
            dst.create_dataset(ip, data=(np.concatenate(ident) if ident else
                                         np.zeros((0, 5), np.int64)),
                               compression="gzip", compression_opts=1)
            dst.create_dataset(cp, data=(np.concatenate(coord) if coord else
                                         np.zeros((0, 2), float)),
                               compression="gzip", compression_opts=1)
            print(f"  {'catch_bond (ACTN2)':26s} {cursor} links across "
                  f"{len(picks)} frames")
        else:
            print(f"  {'catch_bond (ACTN2)':26s} absent, skipped")

        meta = dst.create_group("meta")
        meta.create_dataset("committed_frames", data=np.array([len(picks)]))
        if "meta" in src:
            for name, value in src["meta"].attrs.items():
                try:
                    meta.attrs[name] = value
                except Exception:
                    pass
        # keep the true step numbers so captions stay honest
        if "state/current_step" in src:
            dst.create_dataset("state/current_step",
                               data=np.stack([np.asarray(src["state/current_step"][int(f)])
                                              for f in picks]))
            dst.create_dataset("meta/source_frame", data=np.array(picks))
print(f"wrote {a.out}  {os.path.getsize(a.out)/1048576:.1f} MB")
