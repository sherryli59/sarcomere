#!/usr/bin/env python3
"""Build the deliberately sparse final-snapshot ACTN2 network explorer."""

from __future__ import annotations

import json
from pathlib import Path

from build_packet import RUNS, extract_frames


ROOT = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere")
OUT = ROOT / "analysis/collaborator_model_packet_20260805"
FRAGMENT = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/final-snapshots-3d.html"
)


def build_data() -> dict:
    runs = []
    for label, path in RUNS:
        extracted = extract_frames(path, label)
        runs.append(
            {
                "label": label,
                "box": extracted["box"],
                "actin_length": extracted["actin_length"],
                "myosin_length": extracted["myosin_length"],
                "frame": extracted["frames"][-1],
            }
        )
    return {"runs": runs}


def fragment(data: dict) -> str:
    payload = json.dumps(data, separators=(",", ":"))
    return f'''<div id="final-snapshots-3d-root">
  <div class="viz-controls" aria-label="Final snapshot controls">
    <label class="form-label" for="fs3d-seed">Final 200,000-step snapshot
      <select class="form-select" id="fs3d-seed">
        <option value="0">Seed 2080201</option>
        <option value="1">Seed 2080202</option>
      </select>
    </label>
    <label class="form-check form-switch">
      <input class="form-check-input" id="fs3d-hide-status1" type="checkbox">
      <span class="form-check-label">Hide status-1 actins</span>
    </label>
  </div>
  <div class="viz-row text-small" id="fs3d-legend" aria-label="Color and line legend">
    <span><i class="fs3d-line fs3d-actin-dark"></i>actin polarity +x</span>
    <span><i class="fs3d-line fs3d-actin-light"></i>actin polarity −x</span>
    <span><i class="fs3d-line fs3d-status1"></i>status 1</span>
    <span><i class="fs3d-line fs3d-myosin"></i>attached myosin</span>
  </div>
  <div id="fs3d-canvas-wrap">
    <canvas id="fs3d-canvas" role="img" aria-label="Interactive three-dimensional final snapshot containing ACTN2-linked actins and myosins attached to those actins"></canvas>
  </div>
  <div class="viz-row text-small" id="fs3d-status" aria-live="polite"></div>
</div>
<style>
  #final-snapshots-3d-root {{ width: 100%; color: var(--foreground); }}
  #fs3d-canvas-wrap {{ position: relative; width: 100%; height: 500px; touch-action: none; }}
  #fs3d-canvas {{ display: block; width: 100%; height: 100%; }}
  #fs3d-legend {{ gap: 1rem; color: var(--muted-foreground); }}
  #fs3d-status {{ color: var(--muted-foreground); }}
  .fs3d-line {{ display: inline-block; width: 1.6rem; margin-right: 0.35rem; vertical-align: 0.2em; border-top: 3px solid currentColor; }}
  .fs3d-actin-dark {{ color: color-mix(in srgb, var(--viz-series-1) 82%, var(--foreground)); }}
  .fs3d-actin-light {{ color: color-mix(in srgb, var(--viz-series-1) 48%, var(--background)); }}
  .fs3d-status1 {{ color: var(--viz-series-1); border-top-style: dashed; }}
  .fs3d-myosin {{ color: var(--viz-series-4); border-top-width: 5px; }}
  @media (max-width: 520px) {{ #fs3d-canvas-wrap {{ height: 380px; }} }}
</style>
<script>
(() => {{
  const root = document.getElementById('final-snapshots-3d-root');
  const canvas = document.getElementById('fs3d-canvas');
  const wrapElement = document.getElementById('fs3d-canvas-wrap');
  const seed = document.getElementById('fs3d-seed');
  const hideStatus1 = document.getElementById('fs3d-hide-status1');
  const status = document.getElementById('fs3d-status');
  const context = canvas.getContext('2d');
  const DATA = {payload};
  let rotationX = -0.34;
  let rotationZ = -0.12;
  let zoom = 1;
  let dragging = false;
  let lastX = 0;
  let lastY = 0;

  function resolvedColor(expression, fallback) {{
    const probe = document.createElement('span');
    probe.style.color = expression;
    probe.style.display = 'none';
    root.appendChild(probe);
    const value = getComputedStyle(probe).color;
    probe.remove();
    return value || fallback;
  }}

  function palette() {{
    return {{
      actinDark: resolvedColor('color-mix(in srgb, var(--viz-series-1) 82%, var(--foreground))', 'rgb(22,72,145)'),
      actinLight: resolvedColor('color-mix(in srgb, var(--viz-series-1) 48%, var(--background))', 'rgb(105,178,225)'),
      myosin: resolvedColor('var(--viz-series-4)', 'rgb(222,145,45)'),
      structure: resolvedColor('var(--muted-foreground)', 'rgb(110,120,128)')
    }};
  }}

  function rotate(point) {{
    const cz = Math.cos(rotationZ), sz = Math.sin(rotationZ);
    const cx = Math.cos(rotationX), sx = Math.sin(rotationX);
    const x1 = cz * point[0] - sz * point[1];
    const y1 = sz * point[0] + cz * point[1];
    return [x1, cx * y1 - sx * point[2], sx * y1 + cx * point[2]];
  }}

  function projector(run) {{
    const width = canvas.clientWidth;
    const height = canvas.clientHeight;
    const scale = Math.min(width / (run.box[0] * 1.12), height / (run.box[2] * 2.6)) * zoom;
    return point => {{
      const q = rotate(point);
      const perspective = 1 / Math.max(0.62, 1 + 0.045 * q[1]);
      return [width * 0.5 + scale * q[0] * perspective,
              height * 0.5 - scale * q[2] * perspective,
              q[1]];
    }};
  }}

  function endpoint(row, halfLength, sign) {{
    return [row[0] + sign * halfLength * row[3],
            row[1] + sign * halfLength * row[4],
            row[2] + sign * halfLength * row[5]];
  }}

  function boxSegments(box) {{
    const half = box.map(value => value * 0.5);
    const corners = [];
    for (const x of [-half[0], half[0]]) for (const y of [-half[1], half[1]]) for (const z of [-half[2], half[2]]) corners.push([x,y,z]);
    const segments = [];
    for (let i = 0; i < corners.length; i += 1) for (let j = i + 1; j < corners.length; j += 1) {{
      if ([0,1,2].filter(axis => corners[i][axis] !== corners[j][axis]).length === 1) segments.push([corners[i], corners[j]]);
    }}
    return segments;
  }}

  function command(project, first, second, stroke, width, alpha, dashed) {{
    const a = project(first), b = project(second);
    return {{a, b, depth: 0.5 * (a[2] + b[2]), stroke, width, alpha, dashed}};
  }}

  function draw() {{
    const run = DATA.runs[Number(seed.value)];
    const frame = run.frame;
    const hide = hideStatus1.checked;
    const colors = palette();
    const project = projector(run);
    const commands = [];
    const actins = frame.actins.filter(row => !hide || row[6] === 2);
    const myosins = frame.myosins.filter(row => !hide || row[6] === 2);

    for (const edge of boxSegments(run.box)) commands.push(command(project, edge[0], edge[1], colors.structure, 1, 0.24, false));
    for (const row of myosins) commands.push(command(
      project,
      endpoint(row, run.myosin_length * 0.5, -1),
      endpoint(row, run.myosin_length * 0.5, 1),
      colors.myosin, 4.5, 0.72, false
    ));
    for (const row of actins) commands.push(command(
      project,
      endpoint(row, run.actin_length * 0.5, -1),
      endpoint(row, run.actin_length * 0.5, 1),
      row[3] >= 0 ? colors.actinDark : colors.actinLight,
      row[6] === 2 ? 4.2 : 3.1,
      row[6] === 2 ? 0.98 : 0.82,
      row[6] === 1
    ));

    commands.sort((a, b) => b.depth - a.depth);
    context.clearRect(0, 0, canvas.clientWidth, canvas.clientHeight);
    context.lineCap = 'round';
    for (const item of commands) {{
      context.globalAlpha = item.alpha;
      context.strokeStyle = item.stroke;
      context.lineWidth = item.width;
      context.setLineDash(item.dashed ? [7, 5] : []);
      context.beginPath();
      context.moveTo(item.a[0], item.a[1]);
      context.lineTo(item.b[0], item.b[1]);
      context.stroke();
    }}
    context.setLineDash([]);
    context.globalAlpha = 1;
    const visibleStatus1 = hide ? 0 : frame.status1;
    status.textContent = `${{run.label}} · final ${{frame.time_ms.toFixed(1)}} ms snapshot · ${{visibleStatus1}} status-1 actins shown · ${{frame.status2}} status-2 actins shown · drag to rotate · wheel to zoom`;
  }}

  function resize() {{
    const ratio = Math.min(globalThis.devicePixelRatio || 1, 2);
    const width = Math.max(1, wrapElement.clientWidth);
    const height = Math.max(1, wrapElement.clientHeight);
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    canvas.style.width = `${{width}}px`;
    canvas.style.height = `${{height}}px`;
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    draw();
  }}

  canvas.addEventListener('pointerdown', event => {{
    dragging = true; lastX = event.clientX; lastY = event.clientY;
    canvas.setPointerCapture(event.pointerId);
  }});
  canvas.addEventListener('pointermove', event => {{
    if (!dragging) return;
    rotationZ += (event.clientX - lastX) * 0.007;
    rotationX += (event.clientY - lastY) * 0.007;
    lastX = event.clientX; lastY = event.clientY; draw();
  }});
  canvas.addEventListener('pointerup', () => {{ dragging = false; }});
  canvas.addEventListener('wheel', event => {{
    event.preventDefault();
    zoom = Math.min(3, Math.max(0.45, zoom * (event.deltaY > 0 ? 0.9 : 1.1)));
    draw();
  }}, {{passive: false}});
  seed.addEventListener('change', draw);
  hideStatus1.addEventListener('change', draw);
  new ResizeObserver(resize).observe(wrapElement);
  resize();
}})();
</script>
'''


def main() -> None:
    data = build_data()
    text = fragment(data)
    if len(text.encode("utf-8")) >= 2_000_000:
        raise RuntimeError("fragment exceeds 2 MB")
    FRAGMENT.parent.mkdir(parents=True, exist_ok=True)
    FRAGMENT.write_text(text)
    print(FRAGMENT)
    print(len(text.encode("utf-8")))


if __name__ == "__main__":
    main()
