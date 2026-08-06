#!/usr/bin/env python3
"""Embed an explicit-titin diagnostic trajectory in the established viewer."""

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def number(value):
    parsed = float(value)
    return None if not math.isfinite(parsed) else round(parsed, 6)


def grouped(path):
    result = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            result[int(row["step"])].append(row)
    return result


def finite_or_none(value):
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def main():
    if len(sys.argv) != 5:
        raise SystemExit("usage: render_thermo_titin_visualization.py PREFIX TEMPLATE OUTPUT TITLE")
    prefix = Path(sys.argv[1])
    template_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])
    title = sys.argv[4]

    with (Path(str(prefix) + ".csv")).open(newline="") as handle:
        main_rows = list(csv.DictReader(handle))
    heads = grouped(Path(str(prefix) + "_heads.csv"))
    links = grouped(Path(str(prefix) + "_links.csv"))
    titins = grouped(Path(str(prefix) + "_titin.csv"))

    # A killed diagnostic can leave the last set of auxiliary CSV rows only
    # partially written.  Render only frames with the full physical inventory.
    complete_rows = [
        row for row in main_rows
        if len(heads[int(row["step"])]) == 140
        and len(titins[int(row["step"])]) == 12
    ]
    if not complete_rows:
        raise RuntimeError("no complete trajectory frames found")

    frames = []
    for row in complete_rows:
        step = int(row["step"])
        actins = [[number(row[f"a{i}_{key}"]) for key in
                   ("x", "y", "z", "ux", "uy", "uz")]
                  for i in range(12)]
        myosins = [[number(row[f"m{i}_{key}"]) for key in
                    ("x", "y", "z", "ux", "uy", "uz")]
                   for i in range(2)]
        head_data = []
        for item in heads[step]:
            head_data.append([
                int(item["myosin"]), int(item["half"]), int(item["head_id"]),
                int(item["state"]), int(item["bound_actin"]),
                number(item["anchor_x"]), number(item["anchor_y"]),
                number(item["anchor_z"]), number(item["binding_x"]),
                number(item["binding_y"]), number(item["binding_z"]),
                number(item["active_force"]),
            ])
        link_data = []
        for item in links[step]:
            link_data.append([
                int(item["link_id"]), int(item["actin_i"]),
                int(item["actin_j"]), int(item["status"]),
                number(item["first_x"]), number(item["first_y"]),
                number(item["first_z"]), number(item["second_x"]),
                number(item["second_y"]), number(item["second_z"]),
                number(item["length"]),
            ])
        titin_data = []
        for item in titins[step]:
            titin_data.append([
                int(item["titin_id"]), int(item["myosin"]), int(item["half"]),
                int(item["track"]), int(item["bound"]),
                int(item["bound_actin"]), number(item["distance"]),
                number(item["extension"]), number(item["force"]),
                number(item["thick_x"]), number(item["thick_y"]),
                number(item["thick_z"]), number(item["actin_x"]),
                number(item["actin_y"]), number(item["actin_z"]),
                number(item["attachment_force"]),
                number(item["attachment_energy"]),
            ])
        frames.append({
            "step": step,
            "t": round(float(row["time"]), 6),
            "L": number(row["Lx"]),
            "aa": number(finite_or_none(row["aa_linked_overlap"])
                         if finite_or_none(row["aa_linked_overlap"]) is not None
                         else row["aa_overlap"]),
            "am": number(row["am_overlap"]),
            "weak": int(row["weak_heads"]),
            "strong": int(row["strong_heads"]),
            "tb": int(row["titin_bound"]),
            "tf": number(row["titin_total_force"]),
            "zcv": number(row["z_radius_cv"]),
            "a": actins,
            "m": myosins,
            "h": head_data,
            "l": link_data,
            "ti": titin_data,
        })

    duration_ms = round(1000 * frames[-1]["t"])
    data = {"duration": duration_ms, "La": 1.0, "Lm": 1.5, "frames": frames}
    fragment = template_path.read_text()
    root_id = re.sub(r"[^a-z0-9-]+", "-", output_path.stem.lower()).strip("-")
    fragment = fragment.replace("tuned-actn2-trajectory", root_id)
    fragment = fragment.replace(
        "Tuned ACTN2 and explicit titin trajectory", title)
    fragment = fragment.replace(
        "A true-scale axial projection and enlarged end-on views show the actual 200 millisecond full-diffusion trajectory with two myosins, twelve actins, ACTN2, titin and myosin heads.",
        f"A true-scale axial projection and enlarged end-on views show the {duration_ms} millisecond full-diffusion annealing trajectory from high actin-actin and low actin-myosin overlap.")
    fragment = fragment.replace('max="200"', f'max="{duration_ms}"')
    fragment = fragment.replace(
        '<span>titin <b data-value="titin">—</b></span>',
        '<span>titin bound <b data-value="titin">—</b></span>\n'
        '    <span>titin force <b data-value="titin-force">—</b></span>\n'
        '    <span>max attachment <b data-value="attach-force">—</b></span>\n'
        '    <span>ACTN2 pairs <b data-value="pair-count">—</b></span>')
    fragment = fragment.replace(
        'root.querySelector(\'[data-value="titin"]\').textContent=Math.round(f.tb)+"/12";',
        'root.querySelector(\'[data-value="titin"]\').textContent=Math.round(f.tb)+"/12";'
        'root.querySelector(\'[data-value="titin-force"]\').textContent=f.tf.toFixed(2)+" pN";'
        'var af=f.ti.reduce(function(v,t){return Math.max(v,t[15]||0);},0);'
        'root.querySelector(\'[data-value="attach-force"]\').textContent=af.toFixed(2)+" pN";'
        'var ps={};f.l.forEach(function(l){if(l[3]>0){var a=Math.min(l[1],l[2]),b=Math.max(l[1],l[2]);ps[a+":"+b]=1;}});'
        'root.querySelector(\'[data-value="pair-count"]\').textContent=Object.keys(ps).length+"/12";')
    fragment = fragment.replace(
        'root.querySelector(\'[data-value="aa"]\').textContent=isFinite(aa)?Math.round(aa*1000)+" nm":"unbound";',
        'root.querySelector(\'[data-value="aa"]\').textContent=Math.round(f.aa*1000)+" nm";')
    fragment = fragment.replace('var axialScale=688/4.7;',
                                'var axialScale=688/DATA.frames[0].L;')
    replacement = "var DATA=" + json.dumps(data, separators=(",", ":")) + ", root="
    fragment, count = re.subn(r"var DATA=.*?, root=", replacement, fragment, count=1)
    if count != 1:
        raise RuntimeError("could not replace embedded trajectory data")
    output_path.write_text(fragment)


if __name__ == "__main__":
    main()
