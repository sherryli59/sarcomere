#!/usr/bin/env python3
"""Collapse the status split into three species for the legend.

Showing three legend entries while six colours are on screen would leave the
pale blue and pale pink unexplained, so the merge applies to the DRAWING as
well: all actin statuses share one blue, both ACTN2 statuses share one pink.
Every segment is still drawn -- only the styling is collapsed.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/html_to_mp4.py")
s = p.read_text()
if "merge_species" in s:
    sys.exit("already patched")

s = s.replace('    ap.add_argument("--all-layers", action="store_true",',
              '    ap.add_argument("--merge-species", action="store_true",\n'
              '                    help="one style and one legend entry per species")\n'
              '    ap.add_argument("--all-layers", action="store_true",')

old = """    draw = [l for l in layers if a.all_layers or l["on"]]"""
new = """    draw = [l for l in layers if a.all_layers or l["on"]]
    if a.merge_species:
        SPECIES = [("actin",  ["actin0", "actin1", "actin2"], "#08306b", 1.7, 0.80),
                   ("myosin", ["myosin"],                     "#f5a45b", 6.6, 0.92),
                   ("ACTN2",  ["aa1", "aa2"],                 "#e87ba4", 1.35, 0.85)]
        present = {l["key"] for l in draw}
        draw = [{"key": keys, "label": label, "color": c, "width": w, "alpha": al,
                 "on": True}
                for label, keys, c, w, al in SPECIES
                if present.intersection(keys)]
        for d in draw:
            d["keys"] = d["key"] if isinstance(d["key"], list) else [d["key"]]
    else:
        for d in draw:
            d["keys"] = [d["key"]]"""
if old not in s:
    sys.exit("draw list anchor not found")
s = s.replace(old, new, 1)

old_loop = """            payload = fr["layers"].get(layer["key"], "")
            if not payload:
                continue
            q = np.frombuffer(base64.b64decode(payload), np.int16).reshape(-1, 2, 3)"""
new_loop = """            chunks = [np.frombuffer(base64.b64decode(fr["layers"][k]), np.int16)
                      for k in layer["keys"]
                      if fr["layers"].get(k)]
            if not chunks:
                continue
            q = np.concatenate(chunks).reshape(-1, 2, 3)"""
if old_loop not in s:
    sys.exit("draw loop anchor not found")
s = s.replace(old_loop, new_loop, 1)
p.write_text(s)
print("--merge-species added")
