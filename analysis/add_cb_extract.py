#!/usr/bin/env python3
"""Carry the ACTN2 tables into the extracted subset file.

local_link_identity and local_link_coordinates are one long ragged array
indexed by local_link_frame_offsets rows of (step, offset, count). Copying a
subset of frames means concatenating each frame's slice and rebuilding the
offsets cumulatively -- copying the offset rows verbatim would point into the
original array and read the wrong bridges.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "extract_frames.py")
s = p.read_text()
if "local_link" in s:
    sys.exit("already added")

old = """        meta = dst.create_group("meta")"""
new = """        # ACTN2 tables: concatenate the per-frame slices and rebuild offsets.
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

        meta = dst.create_group("meta")"""
if old not in s:
    sys.exit("anchor not found")
p.write_text(s.replace(old, new, 1))
print("extract_frames.py now carries the ACTN2 tables")
