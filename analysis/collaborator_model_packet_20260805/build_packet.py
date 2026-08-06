#!/usr/bin/env python3
"""Build the collaborator-facing model brief and successful-run 3-D explorer.

Only the two frozen 100 ms PRODALIGN120 trajectories are used for production
scale results.  In particular, no later collapsed-ACTN2-retention trajectory is
read by this script.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from textwrap import wrap

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_collaborator_mpl")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere")
OUT = ROOT / "analysis/collaborator_model_packet_20260805"
VIS_FRAGMENT = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/successful-production-3d.html"
)

RUNS = [
    (
        "seed 2080201",
        Path("/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/")
        / "PRODALIGN120_A1280_M640_h35_100ms_seed2080201.h5",
    ),
    (
        "seed 2080202",
        Path("/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/")
        / "PRODALIGN120_A1280_M640_h35_100ms_seed2080202.h5",
    ),
]

ANIMATIONS = [
    (
        "minimal-free-actin-walkoff.mp4",
        ROOT / "analysis/free_actin_head_tether_20260805/free-actin-myosin-walk-off.mp4",
        "Free actin on one myosin: weak/strong heads, material binding sites, connector strain, transverse tether, force, torque, and walk-off.",
    ),
    (
        "minimal-one-sided-crosslink.mp4",
        ROOT / "cpp/diagnostics/self_loaded_balance_20260804/animations/one_sided_300ms.mp4",
        "One motor-bound actin crosslinked to an unsupported actin: the pair co-moves while ACTN2 remains turnover-capable.",
    ),
    (
        "minimal-force-balanced.mp4",
        ROOT / "cpp/diagnostics/self_loaded_balance_20260804/axial_inline/perfect_inline_1000ms.mp4",
        "Two myosins and four actins from the registered start: endogenous opposing motor load and a stalled A-A register.",
    ),
    (
        "minimal-annealing-start.mp4",
        ROOT / "cpp/diagnostics/self_loaded_balance_20260804/axial_inline_annealing/annealing_1000ms.mp4",
        "The same balanced cell started with excess A-A overlap and reduced A-M overlap: sliding toward the preferred A-A register before stalling.",
    ),
]


def rounded(values: np.ndarray, digits: int = 3) -> list:
    return np.round(np.asarray(values, dtype=float), digits).tolist()


def pbc_link_points(
    centers: np.ndarray,
    directions: np.ndarray,
    i: int,
    j: int,
    s_i: float,
    s_j: float,
    box: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    first = centers[i] + s_i * directions[i]
    second_raw = centers[j] + s_j * directions[j]
    delta = second_raw - first
    delta -= box * np.round(delta / box)
    return first, first + delta


def extract_frames(path: Path, label: str) -> dict:
    with h5py.File(path, "r", locking=False) as h:
        nframes = int(h["meta/committed_frames"][0])
        sample = np.unique(np.linspace(0, nframes - 1, 10).round().astype(int))
        box = np.array(
            [h["meta"].attrs["box_x"], h["meta"].attrs["box_y"], h["meta"].attrs["box_z"]],
            dtype=float,
        )
        actin_length = float(h["meta"].attrs["actin_length"])
        myosin_length = float(h["meta"].attrs["myosin_length"])
        dt = float(h["meta"].attrs["dt"])
        offsets = np.asarray(h["catch_bond/local_link_frame_offsets"], dtype=int)
        identities = h["catch_bond/local_link_identity"]
        coordinates = h["catch_bond/local_link_coordinates"]
        frames = []

        for frame_index in sample:
            ac = np.asarray(h["actin/center"][frame_index], dtype=float)
            ad = np.asarray(h["actin/direction"][frame_index], dtype=float)
            mc = np.asarray(h["myosin/center"][frame_index], dtype=float)
            md = np.asarray(h["myosin/direction"][frame_index], dtype=float)
            cb = np.asarray(h["actin/cb_status"][frame_index, :, 0], dtype=int)
            head_state = np.asarray(h["state/am_head_state_current"][frame_index], dtype=int)
            head_actin = np.asarray(h["state/am_head_bound_actin_current"][frame_index], dtype=int)

            selected_actins = np.flatnonzero(cb > 0)
            actin_rows = []
            for a in selected_actins:
                actin_rows.append(
                    rounded(np.r_[ac[a], ad[a], cb[a]], 3)
                )

            selected_set = set(int(a) for a in selected_actins)
            myosin_rows = []
            for m in range(len(mc)):
                bound = head_actin[m][head_state[m] > 0]
                target_statuses = [int(cb[a]) for a in bound if int(a) in selected_set]
                if not target_statuses:
                    continue
                myosin_rows.append(
                    rounded(np.r_[mc[m], md[m], max(target_statuses)], 3)
                )

            _, start, count = offsets[frame_index]
            ident = np.asarray(identities[start : start + count], dtype=int)
            coords = np.asarray(coordinates[start : start + count], dtype=float)
            link_rows = []
            for row, coord in zip(ident, coords):
                i, j, status = int(row[2]), int(row[3]), int(row[4])
                first, second = pbc_link_points(
                    ac, ad, i, j, float(coord[0]), float(coord[1]), box
                )
                link_rows.append(rounded(np.r_[first, second, status], 3))

            step = int(h["state/current_step"][frame_index, 0])
            frames.append(
                {
                    "time_ms": round(1000.0 * step * dt, 2),
                    "actins": actin_rows,
                    "myosins": myosin_rows,
                    "links": link_rows,
                    "status1": int(np.count_nonzero(cb == 1)),
                    "status2": int(np.count_nonzero(cb == 2)),
                }
            )

        return {
            "label": label,
            "box": rounded(box, 3),
            "actin_length": actin_length,
            "myosin_length": myosin_length,
            "frames": frames,
        }


def make_fragment(datasets: list[dict]) -> str:
    data = json.dumps({"runs": datasets}, separators=(",", ":"))
    return f'''<div id="successful-production-3d-root">
  <div class="viz-controls" aria-label="Three-dimensional trajectory controls">
    <label class="form-label" for="sp3d-seed">Successful run
      <select class="form-select" id="sp3d-seed">
        <option value="0">Seed 2080201</option>
        <option value="1">Seed 2080202</option>
      </select>
    </label>
    <label class="form-label" for="sp3d-time"><span id="sp3d-time-label">Time</span>
      <input class="form-range" id="sp3d-time" type="range" min="0" max="9" value="9" step="1">
    </label>
    <button class="btn" type="button" id="sp3d-play" aria-pressed="false">Play</button>
    <label class="form-check">
      <input class="form-check-input" id="sp3d-mature" type="checkbox">
      <span class="form-check-label">Mature network only</span>
    </label>
    <label class="form-check">
      <input class="form-check-input" id="sp3d-links" type="checkbox" checked>
      <span class="form-check-label">ACTN2 connectors</span>
    </label>
    <label class="form-check">
      <input class="form-check-input" id="sp3d-barbed" type="checkbox" checked>
      <span class="form-check-label">Barbed ends</span>
    </label>
  </div>
  <div id="sp3d-canvas-wrap">
    <canvas id="sp3d-canvas" role="img" aria-label="Interactive three-dimensional view of ACTN2-bound actins and their attached myosins in the successful 1280-actin production run"></canvas>
  </div>
  <div class="viz-row text-small" id="sp3d-status" aria-live="polite"></div>
</div>
<style>
  #successful-production-3d-root {{ width: 100%; color: var(--foreground); }}
  #sp3d-canvas-wrap {{ position: relative; width: 100%; height: 520px; touch-action: none; }}
  #sp3d-canvas {{ display: block; width: 100%; height: 100%; }}
  #sp3d-time {{ min-width: 210px; }}
  #sp3d-status {{ justify-content: space-between; color: var(--muted-foreground); }}
  @media (max-width: 520px) {{
    #sp3d-canvas-wrap {{ height: 390px; }}
    #sp3d-time {{ min-width: 150px; }}
  }}
</style>
<script type="module">
  import * as THREE from 'https://unpkg.com/three@0.166.1/build/three.module.js';

  const root = document.getElementById('successful-production-3d-root');
  const canvas = document.getElementById('sp3d-canvas');
  const wrap = document.getElementById('sp3d-canvas-wrap');
  const seed = document.getElementById('sp3d-seed');
  const slider = document.getElementById('sp3d-time');
  const play = document.getElementById('sp3d-play');
  const mature = document.getElementById('sp3d-mature');
  const linksToggle = document.getElementById('sp3d-links');
  const barbedToggle = document.getElementById('sp3d-barbed');
  const status = document.getElementById('sp3d-status');
  const timeLabel = document.getElementById('sp3d-time-label');
  const DATA = {data};

  function cssColor(variable) {{
    const probe = document.createElement('span');
    probe.style.color = `var(${{variable}})`;
    probe.style.display = 'none';
    root.appendChild(probe);
    const value = getComputedStyle(probe).color;
    probe.remove();
    return new THREE.Color(value || 'rgb(100,100,100)');
  }}

  const palette = {{
    status1: cssColor('--viz-series-2'),
    status2: cssColor('--viz-series-3'),
    myosin1: cssColor('--viz-series-1'),
    myosin2: cssColor('--viz-series-4'),
    link1: cssColor('--viz-series-2'),
    link2: cssColor('--viz-series-3'),
    barbed: cssColor('--viz-series-5'),
    structure: cssColor('--muted-foreground')
  }};

  const renderer = new THREE.WebGLRenderer({{canvas, antialias: true, alpha: true}});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(38, 1, 0.01, 50);
  camera.position.set(0, -7.1, 4.8);
  const world = new THREE.Group();
  world.rotation.x = -0.12;
  world.rotation.z = -0.06;
  scene.add(world);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x777777, 2.2));
  const key = new THREE.DirectionalLight(0xffffff, 2.8);
  key.position.set(2, -4, 7);
  scene.add(key);

  let dynamicObjects = [];
  let boxObject = null;
  let timer = null;
  const yAxis = new THREE.Vector3(0, 1, 0);
  const tmpMatrix = new THREE.Matrix4();
  const tmpPosition = new THREE.Vector3();
  const tmpDirection = new THREE.Vector3();
  const tmpQuaternion = new THREE.Quaternion();
  const tmpScale = new THREE.Vector3();

  function disposeDynamic() {{
    for (const obj of dynamicObjects) {{
      world.remove(obj);
      if (obj.geometry) obj.geometry.dispose();
      if (obj.material) obj.material.dispose();
    }}
    dynamicObjects = [];
  }}

  function addCylinders(rows, length, radius, color) {{
    if (!rows.length) return;
    const geometry = new THREE.CylinderGeometry(radius, radius, 1, 8, 1, false);
    const material = new THREE.MeshStandardMaterial({{color, roughness: 0.52, metalness: 0.02}});
    const mesh = new THREE.InstancedMesh(geometry, material, rows.length);
    rows.forEach((row, index) => {{
      tmpPosition.set(row[0], row[1], row[2]);
      tmpDirection.set(row[3], row[4], row[5]).normalize();
      tmpQuaternion.setFromUnitVectors(yAxis, tmpDirection);
      tmpScale.set(1, length, 1);
      tmpMatrix.compose(tmpPosition, tmpQuaternion, tmpScale);
      mesh.setMatrixAt(index, tmpMatrix);
    }});
    mesh.instanceMatrix.needsUpdate = true;
    world.add(mesh);
    dynamicObjects.push(mesh);
  }}

  function addBarbedEnds(rows, length) {{
    if (!rows.length || !barbedToggle.checked) return;
    const geometry = new THREE.SphereGeometry(0.042, 10, 7);
    const material = new THREE.MeshStandardMaterial({{color: palette.barbed, roughness: 0.5}});
    const mesh = new THREE.InstancedMesh(geometry, material, rows.length);
    rows.forEach((row, index) => {{
      tmpPosition.set(
        row[0] - 0.5 * length * row[3],
        row[1] - 0.5 * length * row[4],
        row[2] - 0.5 * length * row[5]
      );
      tmpMatrix.compose(tmpPosition, new THREE.Quaternion(), new THREE.Vector3(1, 1, 1));
      mesh.setMatrixAt(index, tmpMatrix);
    }});
    mesh.instanceMatrix.needsUpdate = true;
    world.add(mesh);
    dynamicObjects.push(mesh);
  }}

  function addLinks(rows, wantedStatus, color) {{
    if (!linksToggle.checked) return;
    const positions = [];
    for (const row of rows) {{
      if (row[6] !== wantedStatus) continue;
      positions.push(...row.slice(0, 6));
    }}
    if (!positions.length) return;
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    const material = new THREE.LineBasicMaterial({{color, transparent: true, opacity: 0.95}});
    const lines = new THREE.LineSegments(geometry, material);
    world.add(lines);
    dynamicObjects.push(lines);
  }}

  function drawBox(box) {{
    if (boxObject) {{ world.remove(boxObject); boxObject.geometry.dispose(); boxObject.material.dispose(); }}
    const geometry = new THREE.EdgesGeometry(new THREE.BoxGeometry(box[0], box[1], box[2]));
    const material = new THREE.LineBasicMaterial({{color: palette.structure, transparent: true, opacity: 0.32}});
    boxObject = new THREE.LineSegments(geometry, material);
    world.add(boxObject);
  }}

  function update() {{
    disposeDynamic();
    const run = DATA.runs[Number(seed.value)];
    slider.max = String(run.frames.length - 1);
    const frame = run.frames[Math.min(Number(slider.value), run.frames.length - 1)];
    const onlyMature = mature.checked;
    const actin1 = onlyMature ? [] : frame.actins.filter(row => row[6] === 1);
    const actin2 = frame.actins.filter(row => row[6] === 2);
    const myosin1 = onlyMature ? [] : frame.myosins.filter(row => row[6] === 1);
    const myosin2 = frame.myosins.filter(row => row[6] === 2);
    addCylinders(actin1, run.actin_length, 0.020, palette.status1);
    addCylinders(actin2, run.actin_length, 0.024, palette.status2);
    addCylinders(myosin1, run.myosin_length, 0.030, palette.myosin1);
    addCylinders(myosin2, run.myosin_length, 0.036, palette.myosin2);
    addBarbedEnds(actin1.concat(actin2), run.actin_length);
    if (!onlyMature) addLinks(frame.links, 1, palette.link1);
    addLinks(frame.links, 2, palette.link2);
    drawBox(run.box);
    timeLabel.textContent = `Time: ${{frame.time_ms.toFixed(1)}} ms`;
    status.textContent = `${{run.label}} · ACTN2-bound actins: status 1 = ${{frame.status1}}, status 2 = ${{frame.status2}} · drag to rotate · wheel to zoom`;
  }}

  function resize() {{
    const width = Math.max(1, wrap.clientWidth);
    const height = Math.max(1, wrap.clientHeight);
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
  }}

  let dragging = false;
  let lastX = 0;
  let lastY = 0;
  canvas.addEventListener('pointerdown', event => {{
    dragging = true; lastX = event.clientX; lastY = event.clientY;
    canvas.setPointerCapture(event.pointerId);
  }});
  canvas.addEventListener('pointermove', event => {{
    if (!dragging) return;
    const dx = event.clientX - lastX;
    const dy = event.clientY - lastY;
    world.rotation.z += dx * 0.006;
    world.rotation.x += dy * 0.006;
    lastX = event.clientX; lastY = event.clientY;
  }});
  canvas.addEventListener('pointerup', () => {{ dragging = false; }});
  canvas.addEventListener('wheel', event => {{
    event.preventDefault();
    camera.position.multiplyScalar(event.deltaY > 0 ? 1.08 : 0.92);
    const radius = camera.position.length();
    if (radius < 4.0) camera.position.setLength(4.0);
    if (radius > 18.0) camera.position.setLength(18.0);
  }}, {{passive: false}});

  function stop() {{
    if (timer) clearInterval(timer);
    timer = null;
    play.textContent = 'Play';
    play.setAttribute('aria-pressed', 'false');
  }}
  play.addEventListener('click', () => {{
    if (timer) {{ stop(); return; }}
    play.textContent = 'Pause';
    play.setAttribute('aria-pressed', 'true');
    timer = setInterval(() => {{
      const next = (Number(slider.value) + 1) % (Number(slider.max) + 1);
      slider.value = String(next);
      update();
    }}, 650);
  }});
  seed.addEventListener('change', () => {{ slider.value = slider.max; update(); }});
  slider.addEventListener('input', update);
  mature.addEventListener('change', update);
  linksToggle.addEventListener('change', update);
  barbedToggle.addEventListener('change', update);
  new ResizeObserver(resize).observe(wrap);
  resize();
  update();
  function animate() {{ renderer.render(scene, camera); requestAnimationFrame(animate); }}
  animate();
</script>
'''


def make_offline_fragment(datasets: list[dict]) -> str:
    """Dependency-free canvas renderer that works from a downloaded file."""
    data = json.dumps({"runs": datasets}, separators=(",", ":"))
    return f'''<div id="successful-production-3d-offline-root">
  <div class="viz-controls" aria-label="Three-dimensional trajectory controls">
    <label class="form-label" for="sp3do-seed">Successful run
      <select class="form-select" id="sp3do-seed">
        <option value="0">Seed 2080201</option>
        <option value="1">Seed 2080202</option>
      </select>
    </label>
    <label class="form-label" for="sp3do-time"><span id="sp3do-time-label">Time</span>
      <input class="form-range" id="sp3do-time" type="range" min="0" max="9" value="9" step="1">
    </label>
    <button class="btn" type="button" id="sp3do-play" aria-pressed="false">Play</button>
    <label class="form-check">
      <input class="form-check-input" id="sp3do-mature" type="checkbox">
      <span class="form-check-label">Mature network only</span>
    </label>
    <label class="form-check">
      <input class="form-check-input" id="sp3do-links" type="checkbox" checked>
      <span class="form-check-label">ACTN2 connectors</span>
    </label>
    <label class="form-check">
      <input class="form-check-input" id="sp3do-barbed" type="checkbox" checked>
      <span class="form-check-label">Barbed ends</span>
    </label>
  </div>
  <div id="sp3do-canvas-wrap">
    <canvas id="sp3do-canvas" role="img" aria-label="Interactive three-dimensional view of ACTN2-bound actins and their attached myosins in the successful 1280-actin production run"></canvas>
  </div>
  <div class="viz-row text-small" id="sp3do-status" aria-live="polite"></div>
</div>
<style>
  #successful-production-3d-offline-root {{ width: 100%; color: var(--foreground); }}
  #sp3do-canvas-wrap {{ position: relative; width: 100%; height: 520px; touch-action: none; }}
  #sp3do-canvas {{ display: block; width: 100%; height: 100%; }}
  #sp3do-time {{ min-width: 210px; }}
  #sp3do-status {{ justify-content: space-between; color: var(--muted-foreground); }}
  @media (max-width: 520px) {{
    #sp3do-canvas-wrap {{ height: 390px; }}
    #sp3do-time {{ min-width: 150px; }}
  }}
</style>
<script>
(() => {{
  const root = document.getElementById('successful-production-3d-offline-root');
  const canvas = document.getElementById('sp3do-canvas');
  const wrapElement = document.getElementById('sp3do-canvas-wrap');
  const seed = document.getElementById('sp3do-seed');
  const slider = document.getElementById('sp3do-time');
  const play = document.getElementById('sp3do-play');
  const mature = document.getElementById('sp3do-mature');
  const linksToggle = document.getElementById('sp3do-links');
  const barbedToggle = document.getElementById('sp3do-barbed');
  const status = document.getElementById('sp3do-status');
  const timeLabel = document.getElementById('sp3do-time-label');
  const context = canvas.getContext('2d');
  const DATA = {data};
  let rotationX = -0.34;
  let rotationZ = -0.12;
  let zoom = 1;
  let timer = null;
  let dragging = false;
  let lastX = 0;
  let lastY = 0;

  function color(variable, fallback) {{
    const value = getComputedStyle(root).getPropertyValue(variable).trim();
    return value || fallback;
  }}

  function palette() {{
    return {{
      background: color('--background', 'rgb(255,255,255)'),
      status1: color('--viz-series-2', 'rgb(24,150,112)'),
      status2: color('--viz-series-3', 'rgb(196,61,119)'),
      myosin1: color('--viz-series-1', 'rgb(42,112,164)'),
      myosin2: color('--viz-series-4', 'rgb(222,145,45)'),
      barbed: color('--viz-series-5', 'rgb(232,91,47)'),
      structure: color('--muted-foreground', 'rgb(110,120,128)')
    }};
  }}

  function rotate(point) {{
    const cz = Math.cos(rotationZ), sz = Math.sin(rotationZ);
    const cx = Math.cos(rotationX), sx = Math.sin(rotationX);
    const x1 = cz * point[0] - sz * point[1];
    const y1 = sz * point[0] + cz * point[1];
    const y2 = cx * y1 - sx * point[2];
    const z2 = sx * y1 + cx * point[2];
    return [x1, y2, z2];
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
    const h = box.map(value => value * 0.5);
    const corners = [];
    for (const x of [-h[0], h[0]]) for (const y of [-h[1], h[1]]) for (const z of [-h[2], h[2]]) corners.push([x,y,z]);
    const segments = [];
    for (let i = 0; i < corners.length; i += 1) {{
      for (let j = i + 1; j < corners.length; j += 1) {{
        const differences = [0,1,2].filter(axis => corners[i][axis] !== corners[j][axis]).length;
        if (differences === 1) segments.push([corners[i], corners[j]]);
      }}
    }}
    return segments;
  }}

  function segment(project, first, second, stroke, width, alpha, kind) {{
    const a = project(first), b = project(second);
    return {{a, b, depth: 0.5 * (a[2] + b[2]), stroke, width, alpha, kind}};
  }}

  function draw() {{
    const run = DATA.runs[Number(seed.value)];
    slider.max = String(run.frames.length - 1);
    const frame = run.frames[Math.min(Number(slider.value), run.frames.length - 1)];
    const onlyMature = mature.checked;
    const colors = palette();
    const project = projector(run);
    const commands = [];
    const actins = frame.actins.filter(row => !onlyMature || row[6] === 2);
    const myosins = frame.myosins.filter(row => !onlyMature || row[6] === 2);

    for (const edge of boxSegments(run.box)) commands.push(segment(project, edge[0], edge[1], colors.structure, 1, 0.28, 'box'));
    for (const row of myosins) {{
      const stroke = row[6] === 2 ? colors.myosin2 : colors.myosin1;
      commands.push(segment(project, endpoint(row, run.myosin_length * 0.5, -1), endpoint(row, run.myosin_length * 0.5, 1), stroke, row[6] === 2 ? 6 : 4, 0.78, 'myosin'));
    }}
    for (const row of actins) {{
      const stroke = row[6] === 2 ? colors.status2 : colors.status1;
      commands.push(segment(project, endpoint(row, run.actin_length * 0.5, -1), endpoint(row, run.actin_length * 0.5, 1), stroke, row[6] === 2 ? 5 : 3, 0.94, 'actin'));
    }}
    if (linksToggle.checked) {{
      for (const row of frame.links) {{
        if (onlyMature && row[6] !== 2) continue;
        commands.push(segment(project, row.slice(0,3), row.slice(3,6), row[6] === 2 ? colors.status2 : colors.status1, row[6] === 2 ? 2.2 : 1.3, 0.85, 'link'));
      }}
    }}
    commands.sort((a, b) => b.depth - a.depth);
    context.clearRect(0, 0, canvas.clientWidth, canvas.clientHeight);
    context.lineCap = 'round';
    for (const command of commands) {{
      context.globalAlpha = command.alpha;
      context.strokeStyle = command.stroke;
      context.lineWidth = command.width;
      context.beginPath();
      context.moveTo(command.a[0], command.a[1]);
      context.lineTo(command.b[0], command.b[1]);
      context.stroke();
    }}
    if (barbedToggle.checked) {{
      context.globalAlpha = 1;
      context.fillStyle = colors.barbed;
      for (const row of actins) {{
        const p = project(endpoint(row, run.actin_length * 0.5, -1));
        context.beginPath(); context.arc(p[0], p[1], row[6] === 2 ? 3.8 : 3, 0, 2 * Math.PI); context.fill();
      }}
    }}
    context.globalAlpha = 1;
    timeLabel.textContent = `Time: ${{frame.time_ms.toFixed(1)}} ms`;
    status.textContent = `${{run.label}} · ACTN2-bound actins: status 1 = ${{frame.status1}}, status 2 = ${{frame.status2}} · drag to rotate · wheel to zoom`;
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

  function stop() {{
    if (timer) clearInterval(timer);
    timer = null; play.textContent = 'Play'; play.setAttribute('aria-pressed', 'false');
  }}
  play.addEventListener('click', () => {{
    if (timer) {{ stop(); return; }}
    play.textContent = 'Pause'; play.setAttribute('aria-pressed', 'true');
    timer = setInterval(() => {{
      slider.value = String((Number(slider.value) + 1) % (Number(slider.max) + 1));
      draw();
    }}, 650);
  }});
  seed.addEventListener('change', () => {{ slider.value = slider.max; draw(); }});
  slider.addEventListener('input', draw);
  mature.addEventListener('change', draw);
  linksToggle.addEventListener('change', draw);
  barbedToggle.addEventListener('change', draw);
  new ResizeObserver(resize).observe(wrapElement);
  resize();
}})();
</script>
'''


def page_setup(title: str, subtitle: str | None = None):
    fig = plt.figure(figsize=(8.5, 11), facecolor="white")
    fig.text(0.08, 0.945, title, fontsize=22, weight="bold", color="#183449")
    if subtitle:
        fig.text(0.08, 0.912, subtitle, fontsize=10.5, color="#5c6d78")
    fig.lines.append(plt.Line2D([0.08, 0.92], [0.895, 0.895], color="#d5dde2", lw=1))
    return fig


def paragraph(fig, x, y, text, width=92, size=10.5, color="#263842", leading=0.022):
    lines = []
    for block in text.split("\n"):
        lines.extend(wrap(block, width=width) or [""])
    fig.text(x, y, "\n".join(lines), fontsize=size, color=color, va="top", linespacing=1.35)
    return y - leading * len(lines)


def equation(fig, y, text, size=15):
    fig.text(0.12, y, text, fontsize=size, color="#0c5278", va="center")


def bullet(fig, y, title, body, color="#263842"):
    fig.text(0.095, y, "•", fontsize=15, color="#c43d2b", va="top")
    fig.text(0.12, y, title, fontsize=11, weight="bold", color="#183449", va="top")
    paragraph(fig, 0.12, y - 0.025, body, width=88, size=9.8, color=color, leading=0.020)


def draw_component(fig, x, y, w, h, title, body, edge):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
        transform=fig.transFigure, fc="#f7fafb", ec=edge, lw=1.5
    )
    fig.patches.append(patch)
    fig.text(x + 0.018, y + h - 0.035, title, fontsize=11, weight="bold", color=edge, va="top")
    fig.text(x + 0.018, y + h - 0.073, "\n".join(wrap(body, 31)), fontsize=9.2, color="#263842", va="top", linespacing=1.35)


def make_pdf(path: Path):
    plt.rcdefaults()

    def save_page(pdf, fig, index):
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(
            OUT / f"model-brief-page-{index}.png",
            dpi=135, bbox_inches="tight", facecolor="white"
        )

    with PdfPages(path) as pdf:
        # Page 1: scope and one-page model map.
        fig = page_setup(
            "Current coarse-grained sarcomere model",
            "Mechanistic summary and changes since the roughly 29 July working model · 5 August 2026",
        )
        paragraph(
            fig, 0.08, 0.86,
            "Goal: test whether myosin motors, turnover-capable ACTN2 and individual titin-like tethers can assemble and maintain preferred A-A register while preserving finite, dynamic A-M engagement—without prescribing a Z-disc node or an external load.",
            width=104, size=11.2, leading=0.024,
        )
        draw_component(fig, 0.08, 0.59, 0.23, 0.17, "A-M heads", "Spatial crowns bind only axially reachable material sites. Weak heads slide/turn over; strong heads execute a finite load-dependent stroke.", "#256b9a")
        draw_component(fig, 0.385, 0.59, 0.23, 0.17, "ACTN2", "Multiple explicit molecules bind barbed-end zones, anneal while weak, lock at register, and gain lifetime only from endogenous bilateral loading.", "#1f8f62")
        draw_component(fig, 0.69, 0.59, 0.23, 0.17, "Titin-like tethers", "Individual thick-tip domains dynamically bind individual actin barbed ends. Springs are tension-only: compression is force-free.", "#b94776")
        for x0, x1 in [(0.31, 0.385), (0.615, 0.69)]:
            fig.patches.append(FancyArrowPatch((x0, 0.675), (x1, 0.675), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=12, color="#7c8a93", lw=1.2))
        bullet(fig, 0.51, "Dynamics", "Anisotropic overdamped Brownian translation and rotation of finite rods; all molecular forces and torques are reciprocal; filament sterics prevent interpenetration.")
        bullet(fig, 0.42, "No imposed stall", "The force-velocity curve is sampled from actual ACTN2 connector load and actual sliding. Balanced motors slow because their working strokes load the crosslinked network.")
        bullet(fig, 0.33, "No pre-existing Z-disc", "There is no six-actin centroid, Z-node, prescribed ring radius or centroid-to-myosin spring. Any local bundle/register must emerge from molecular interactions.")
        fig.text(0.08, 0.18, "Evidence included in this packet", fontsize=12, weight="bold", color="#183449")
        paragraph(fig, 0.08, 0.145, "Production-scale 3-D views use only the successful frozen model-19 runs: 1,280 actins, 640 myosins, 100 ms, 35 heads per half, alignment stiffness 120, seeds 2080201 and 2080202. Later collapsed-retention production runs are excluded. These runs support A-A maturation and nematic alignment; they are not presented as proof that A-M placement is solved.", width=106, size=9.8, color="#5c6d78", leading=0.019)
        save_page(pdf, fig, 1)
        plt.close(fig)

        # Page 2: rod mechanics and motors.
        fig = page_setup("Brownian rods and explicit A-M heads")
        paragraph(fig, 0.08, 0.855, "Actin and myosin are finite spherocylindrical rods. Their mobility is anisotropic and derived from rod geometry, solvent viscosity and temperature.", width=103)
        equation(fig, 0.785, r"$d\mathbf{r}=\mathbf{M}_t\mathbf{F}\,dt+\sqrt{2k_BT\mathbf{M}_t}\,d\mathbf{W}$")
        equation(fig, 0.735, r"$d\mathbf{u}=\left[\mathbf{M}_r\boldsymbol{\tau}\,dt+\sqrt{2k_BT\mathbf{M}_r}\,d\mathbf{W}_r\right]\times\mathbf{u}$")
        paragraph(fig, 0.08, 0.68, "Each myosin half owns spatial head anchors. A detached head can bind actin i only when its anchor projects inside the finite actin segment, lies on the productive myosin half, and is within surface capture reach:", width=104)
        equation(fig, 0.605, r"$|s_i^*|\leq L_A/2,\qquad h_{head}=h_{prod},\qquad d_{surf}\leq d_c$")
        equation(fig, 0.555, r"$k_{on}^{hi}=k_{on}^0\,a_i(t)\exp\!\left[-\frac{(d_{surf}-d_0)^2}{2\sigma_d^2}\right]$")
        paragraph(fig, 0.08, 0.505, "The bound actin coordinate is material and remains fixed until detachment. This prevents heads outside an actin from grabbing its endpoint and makes relative sliding generate strain.", width=104)
        equation(fig, 0.425, r"$D\ \longrightarrow\ W\ \longrightarrow\ S,\qquad W\ \longrightarrow\ D,\qquad S\ \longrightarrow\ D$")
        equation(fig, 0.37, r"$F_{\parallel}=\min\!\left[k_m(q-z)_+,F_s\right],\qquad \dot q=v_0\left(1-\frac{F_{resist}}{F_s}\right)_+$")
        equation(fig, 0.315, r"$U_{\perp}=\frac{1}{2} k_{\perp}\,[\rho-\rho_0-\delta_{slack}]_+^2$")
        equation(fig, 0.26, r"$k_{off}=k_{chem}(F)\exp\!\left[\beta(\theta_{\parallel}U_{\parallel}^{np}+\theta_{\perp}U_{\perp})\right]$")
        paragraph(fig, 0.08, 0.205, "The radial spring is a passive, slack tether on each bound head, not a permanent A-M bond. It appears only while that head is weak or strong and disappears on detachment. The final geometric rupture is an emergency contour guard, not the normal turnover mechanism.", width=104, size=9.8, color="#5c6d78", leading=0.020)
        save_page(pdf, fig, 2)
        plt.close(fig)

        # Page 3: ACTN2.
        fig = page_setup("Explicit ACTN2: annealing, register and endogenous load")
        paragraph(fig, 0.08, 0.855, "An ACTN2 molecule stores two material coordinates on antiparallel actins. Its connector acts in full 3-D, so it constrains both axial register and transverse separation.", width=103)
        equation(fig, 0.79, r"$\mathbf{c}_{\ell}=(\mathbf{r}_j+s_j\mathbf{u}_j)-(\mathbf{r}_i+s_i\mathbf{u}_i),\qquad r_{\ell}=|\mathbf{c}_{\ell}|$")
        equation(fig, 0.735, r"$U_{AA,\ell}=\frac{1}{2} k_{r,\ell}(r_{\ell}-r_0)^2+U_{register}+U_{angle}$")
        paragraph(fig, 0.08, 0.675, "Status 1 is a weak encounter: sites can anneal within the 200 nm barbed-end association zones. Status 2 is a registered, material-site lock near the preferred 0.10 micrometre A-A overlap. Association/unbinding and status are distinct variables.", width=104)
        fig.text(0.11, 0.585, "status 0", fontsize=11, weight="bold", color="#6f7c84")
        fig.text(0.39, 0.585, "status 1", fontsize=11, weight="bold", color="#1f8f62")
        fig.text(0.69, 0.585, "status 2", fontsize=11, weight="bold", color="#b94776")
        fig.patches.append(FancyArrowPatch((0.20, 0.59), (0.36, 0.59), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=12, color="#7c8a93"))
        fig.patches.append(FancyArrowPatch((0.49, 0.59), (0.66, 0.59), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=12, color="#7c8a93"))
        fig.text(0.10, 0.545, "unbound", fontsize=9.5, color="#5c6d78")
        fig.text(0.37, 0.545, "bound, sliding", fontsize=9.5, color="#5c6d78")
        fig.text(0.66, 0.545, "registered, loaded", fontsize=9.5, color="#5c6d78")
        paragraph(fig, 0.08, 0.485, "The load used for maturation is not an applied force and is not inferred from overlap. Productive motor force on each actin is allocated across the ACTN2 arms capable of transmitting it:", width=104)
        equation(fig, 0.405, r"$A_{i\ell}=M_i\,\frac{R_{i\ell}}{\sum_m R_{im}},\qquad \widetilde L_{i\ell}=\min(R_{i\ell},A_{i\ell})$")
        equation(fig, 0.35, r"$L_{\ell}=\min(\widetilde L_{i\ell},\widetilde L_{j\ell})$")
        equation(fig, 0.295, r"$\mathrm{status\ 2}=\mathrm{registered}\times\mathrm{productive\ opposed\ motor\ support}$")
        equation(fig, 0.24, r"$\tau_{\ell}=\tau_0+\mathcal{I}_{s=2}\left(\tau_{dir}+\lambda L_{\ell}\right),\qquad P_{off}=\min(dt/\tau_{\ell},1)$")
        paragraph(fig, 0.08, 0.185, "A separate contour test removes a connector if its actual 3-D length exceeds the rupture limit. This split lets us distinguish stochastic chemistry from geometric failure in diagnostics.", width=104, size=9.8, color="#5c6d78", leading=0.020)
        save_page(pdf, fig, 3)
        plt.close(fig)

        # Page 4: titin.
        fig = page_setup("Explicit titin-like molecules—without a Z-node")
        paragraph(fig, 0.08, 0.855, "Each thick-filament tip carries explicit binding domains. A domain binds one individual actin barbed end, turns over, and can rebind. It never refers to a six-actin centroid.", width=104)
        equation(fig, 0.785, r"$x_T=[d_T-\ell_T]_+$")
        equation(fig, 0.735, r"$U_T=\frac{1}{2} k_Tx_T^2,\qquad F_T=k_Tx_T$")
        paragraph(fig, 0.08, 0.69, "Compression therefore produces no force. A tether becomes mechanically active only after its end-to-end distance exceeds the slack length.", width=104)
        equation(fig, 0.61, r"$k_{on,T}=k_{on,T}^0\,\mathcal{I}(d_T\leq d_c)\,\mathcal{I}(\mathrm{live\ ACTN2\ support})$")
        equation(fig, 0.555, r"$k_{off,T}=k_{off,T}^0\exp\!\left[\beta(1-\theta_T)U_T\right]$")
        paragraph(fig, 0.08, 0.505, "The support gate is a formation/rebinding rule, not a centroid constraint: titin is licensed by a locally crosslinked barbed end, then applies equal and opposite force to that actin and the individual myosin-tip domain. Unsupported serial anchors can release quickly.", width=104)
        bullet(fig, 0.405, "What titin fixes", "It removes the common-translation zero mode of an A-A pair by coupling individual barbed ends to the thick-filament lattice while remaining slack under compression.")
        bullet(fig, 0.31, "What titin does not assume", "No Z-disc position, ring radius, six-actin membership, A-A centroid, or preferred azimuth is supplied to the force law.")
        bullet(fig, 0.215, "Remaining coarse graining", "The binding domain and spring stand in for a titin-like axial tether; actin polymerization and the molecular architecture of a mature Z-disc are not modeled.")
        save_page(pdf, fig, 4)
        plt.close(fig)

        # Page 5: changes and evidence.
        fig = page_setup("What changed over the past week")
        rows = [
            ("A-M capture", "Broad pair/end-point reach", "Material crowns; anchor projection must lie inside finite actin and on the productive half"),
            ("Motor mechanics", "Overlap/status proxy", "Weak/strong heads, fixed material sites, finite stroke, radial tether, reciprocal force and torque"),
            ("A-A model", "One pair-level spring", "Multiple explicit ACTN2 molecules; weak annealing and registered material-site locking"),
            ("Status 2", "cb=2 plus substantial overlap / motor identity", "Registry plus productive opposed motor support; load comes from the actual connector"),
            ("Titin", "Occupancy/centroid-triggered spring", "Individual dynamic tip-to-barbed-end molecules; tension-only; no centroid"),
            ("Motion", "Simpler translational noise", "Anisotropic physical rod mobility, rotational diffusion, torques and sterics"),
            ("Balanced F-V", "Applied diagnostic load", "Endogenous ACTN2 load paired with measured sliding velocity"),
            ("Assembly prior", "Implicit Z-node/centroid stabilizers in sandbox", "Removed: local geometry must emerge from molecular interactions"),
        ]
        ax = fig.add_axes([0.075, 0.33, 0.85, 0.54])
        ax.axis("off")
        table = ax.table(
            cellText=rows,
            colLabels=["Mechanism", "~29 July working model", "Current model"],
            colWidths=[0.18, 0.30, 0.52],
            cellLoc="left",
            colLoc="left",
            bbox=[0, 0, 1, 1],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.4)
        for (r, c), cell in table.get_celld().items():
            cell.set_edgecolor("#d5dde2")
            cell.set_linewidth(0.7)
            cell.PAD = 0.025
            if r == 0:
                cell.set_facecolor("#183449")
                cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
            elif c == 0:
                cell.set_facecolor("#eef4f6")
                cell.get_text().set_weight("bold")
            else:
                cell.set_facecolor("white")
        fig.text(0.08, 0.275, "Successful older scale-up used for the 3-D explorer", fontsize=12, weight="bold", color="#183449")
        fig.text(0.08, 0.232, "seed 2080201", fontsize=9.5, weight="bold", color="#183449")
        fig.text(0.25, 0.232, "342 active ACTN2 links · 185 status-2 actins · S_A = 0.838", fontsize=9.5, color="#263842")
        fig.text(0.08, 0.198, "seed 2080202", fontsize=9.5, weight="bold", color="#183449")
        fig.text(0.25, 0.198, "322 active ACTN2 links · 186 status-2 actins · S_A = 0.893", fontsize=9.5, color="#263842")
        paragraph(fig, 0.08, 0.15, "Interpretation: the run demonstrates robust nematic alignment and maturation near the preferred A-A register. Finite A-M placement remains a model question and is deliberately not claimed here. No result from the later collapsed-ACTN2-retention production trajectories appears in this document.", width=105, size=9.6, color="#5c6d78", leading=0.019)
        save_page(pdf, fig, 5)
        plt.close(fig)


def write_readme():
    animation_lines = "\n".join(
        f"- [{name}](animations/{name}) — {description}" for name, _, description in ANIMATIONS
    )
    text = f"""# Collaborator model packet — 2026-08-05

This packet intentionally excludes every later production trajectory with collapsed ACTN2 retention.

## Model brief

- [Current coarse-grained model with equations](sarcomere-model-brief.pdf)

## Interactive 3-D

- [Successful 1280-actin production explorer](successful-production-3d.html)

The explorer embeds ten frames from each of the two successful 100 ms model-19 runs (`PRODALIGN120`, seeds 2080201 and 2080202). It shows only ACTN2-bound actins and myosins attached to those actins. Status 1 and status 2 are separately colored; the mature-only control removes encounter complexes.

## Minimal animations

{animation_lines}

## Production data provenance

- `/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080201.h5`
- `/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080202.h5`

The 3-D explorer is evidence for A-A maturation and nematic alignment. It is not presented as proof that finite A-M placement is solved.
"""
    (OUT / "README.md").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "animations").mkdir(parents=True, exist_ok=True)
    VIS_FRAGMENT.parent.mkdir(parents=True, exist_ok=True)

    missing = [str(path) for _, path in RUNS if not path.exists()]
    missing += [str(path) for _, path, _ in ANIMATIONS if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing packet inputs:\n" + "\n".join(missing))

    datasets = [extract_frames(path, label) for label, path in RUNS]
    fragment = make_offline_fragment(datasets)
    if len(fragment.encode("utf-8")) >= 2_000_000:
        raise RuntimeError("Interactive fragment exceeds the 2 MB limit")
    VIS_FRAGMENT.write_text(fragment)

    make_pdf(OUT / "sarcomere-model-brief.pdf")
    for name, source, _ in ANIMATIONS:
        shutil.copy2(source, OUT / "animations" / name)
    write_readme()

    manifest = {
        "production_sources": [str(path) for _, path in RUNS],
        "excluded_result_class": "later production runs with collapsed ACTN2 retention",
        "interactive_fragment": str(VIS_FRAGMENT),
        "fragment_bytes": len(fragment.encode("utf-8")),
        "animations": [name for name, _, _ in ANIMATIONS],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
