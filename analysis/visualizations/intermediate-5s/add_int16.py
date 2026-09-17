#!/usr/bin/env python3
"""Halve the payload with int16 coordinates.

Coordinates live in a 6 x 2 x 2 um box, so a 16-bit fixed-point mapping of
[-L/2, L/2] per axis resolves 6000/65535 = 0.09 nm on the long axis -- three
orders of magnitude finer than the 15 nm filament radius, i.e. far below one
screen pixel. float32 spends 4 bytes to carry precision the display cannot
show; int16 spends 2 and buys twice the frames.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/intermediate-5s-animation-template.html")
s = p.read_text()
if "decodeSegments" in s:
    sys.exit("already added")

old = """      decodedFrames.set(run.key,run.frames.map(frame => Object.fromEntries(
        Object.entries(frame.layers).map(([key,value]) => [key,decode(value,Float32Array)])
      )));"""
new = """      // Coordinates are int16 fixed-point over [-L/2, L/2] per axis; see the
      // builder. Expand to Float32Array once at load so the draw path is
      // unchanged.
      const halfBox = run.box.map(v => v / 2);
      const decodeSegments = value => {
        const q = decode(value,Int16Array);
        const out = new Float32Array(q.length);
        for (let i = 0; i < q.length; i++) out[i] = q[i] * halfBox[i % 3] / 32767;
        return out;
      };
      decodedFrames.set(run.key,run.frames.map(frame => Object.fromEntries(
        Object.entries(frame.layers).map(([key,value]) => [key,decodeSegments(value)])
      )));"""
if old not in s:
    sys.exit("decode block not found")
p.write_text(s.replace(old, new, 1))
print("template now decodes int16 coordinates")
