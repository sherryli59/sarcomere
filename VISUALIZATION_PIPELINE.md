# Sarcomere interactive visualization pipeline

This is the model-independent instruction for generating the established
interactive sarcomere snapshot viewer. Codex, Claude Code, or a human can use
the same workflow. The result is a single inline HTML file with no binary side
files: HDF5 geometry is extracted in Python, packed as `float32`, embedded as
ASCII base64, and drawn by a rotatable canvas viewer.

## Canonical implementation

The current source of truth is:

- extractor utilities: `/home/users/xsli/.codex/visualizations/2026/08/26/01a03f41-9578-7e83-9ee5-f379f3493c26/build-simulations.py`
- compact multi-run builder: `/home/users/xsli/.codex/visualizations/2026/08/26/01a03f41-9578-7e83-9ee5-f379f3493c26/build-a2400-latest-snapshots.py`
- interactive template: `/home/users/xsli/.codex/visualizations/2026/08/26/01a03f41-9578-7e83-9ee5-f379f3493c26/a2400-latest-snapshots-template.html`
- working single-run example: `/home/users/xsli/.codex/visualizations/2026/08/26/01a03f41-9578-7e83-9ee5-f379f3493c26/a1600-latest-snapshot.html`

Treat the Python builder and HTML template as source files. Treat the rendered
HTML as a generated artifact; do not hand-edit a large embedded payload.

For a new comparison, copy the builder and template to a new, clearly named
pair, change only the run table/title/output path initially, then extend them
if the requested visualization genuinely needs new layers.

## Required behavior

Every snapshot viewer must provide:

1. the latest **committed** HDF5 frame, never a partially written frame;
2. periodic-boundary-safe filament and ACTN2 segments;
3. run-selection controls when more than one trajectory is shown;
4. independent toggles for every visual layer;
5. drag-to-rotate, wheel-to-zoom, responsive resizing, and depth sorting;
6. the same species colors and thickness scale across all runs being compared;
7. frame, step, time, counts, M order, and Z order in visible text;
8. a self-contained, ASCII-safe embedded data payload.

## HDF5 input contract

Open trajectories read-only with HDF5 file locking disabled:

```python
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
with h5py.File(path, "r", locking=False) as handle:
    ...
```

The core datasets are:

```text
/state/current_step
/actin/center
/actin/direction
/actin/cb_status
/myosin/center
/myosin/direction
/catch_bond/local_link_frame_offsets
/catch_bond/local_link_identity
/catch_bond/local_link_coordinates
```

Read `box_x`, `box_y`, `box_z`, `dt`, `actin_length`, and `myosin_length`
from `/meta` attributes.

Choose the frame count as the minimum of every required dataset's first
dimension and `/meta/committed_frames[0]` when that dataset exists. The latest
frame is `frame_count - 1`. This rule is mandatory for live trajectories.

The local-link offset row is `(step, offset, count)`. Slice both identity and
coordinate arrays with `offset:offset+count`; status is identity column 4 and
the linked actin IDs are columns 2 and 3.

## Geometry extraction

Represent every filament or crosslink as six floats:

```text
x0, y0, z0, x1, y1, z1
```

Build actin and myosin endpoints from center, direction, and half-length.
Split actin into `actin0`, `actin1`, and `actin2` using `actin/cb_status`.
Build `aa1` and `aa2` from the material coordinates stored for each live local
ACTN2 link. Build `barbed` as three-float point records.

For an ACTN2 endpoint pair, minimum-image the second endpoint relative to the
first before drawing:

```python
end = start + (raw_end - start - box * np.rint((raw_end - start) / box))
```

Then pass all filament and ACTN2 segments through
`build-simulations.py::clip_periodic_segments`. Do not draw a long line across
the box when a segment crosses a periodic boundary.

## M and Z order

Keep the definition identical across every plot and snapshot:

```python
def band_order(x, box_x):
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(
        2j * np.pi * wrapped[..., None] * harmonics[None, None, :] / box_x
    ).mean(axis=1)
    return np.max(np.abs(phase), axis=1)
```

M order is `band_order(myosin_center_x, box_x)`. Z order is the maximum of
the same measure evaluated at the two actin-tip populations. Do not substitute
actin centers or only status-2 actins.

## Layer style contract

Use one hue per species. Status is encoded by opacity and width, not by changing
the species hue. Keep these settings unchanged across comparison panels:

```javascript
const layers = [
  {key:'actin0', label:'unlinked actin', color:'var(--viz-series-3)', width:0.9,  alpha:0.15, kind:'segments'},
  {key:'actin1', label:'status-1 actin', color:'var(--viz-series-3)', width:1.2,  alpha:0.40, kind:'segments'},
  {key:'actin2', label:'status-2 actin', color:'var(--viz-series-3)', width:1.7,  alpha:0.90, kind:'segments'},
  {key:'myosin', label:'myosin',          color:'var(--viz-series-4)', width:4.8,  alpha:0.86, kind:'segments'},
  {key:'aa1',    label:'ACTN2 status-1',  color:'var(--viz-series-5)', width:0.8,  alpha:0.25, kind:'segments'},
  {key:'aa2',    label:'ACTN2 status-2',  color:'var(--viz-series-5)', width:1.35, alpha:0.82, kind:'segments'},
  {key:'barbed', label:'barbed ends',     color:'var(--viz-series-6)', width:2.5,  alpha:0.86, kind:'points'}
];
```

In particular, myosin must remain visibly thicker than actin and must have the
same 4.8-pixel width in every run. Both ACTN2 states use the same ACTN2 hue.

Use theme variables (`--foreground`, `--muted-foreground`, `--border`, and
`--viz-series-*`) rather than hard-coded page background or text colors.

**Every theme variable must carry a literal fallback.** The viewer resolves
colors by assigning them to a probe element and reading back the computed
`color`. If a variable is undefined the declaration is invalid at
computed-value time and the probe returns the initial color, so *every layer
renders black* -- the failure reported against the motor-arm viewer on
2026-09-01, and present in the canonical template before then. A host theme
must still win, so do not simply hard-code the hue; give each variable a
fallback the artifact defines under its own prefix:

```css
:root {
  --sarc-series-3: #1baf7a;   /* actin   */
  --sarc-series-4: #eda100;   /* myosin  */
  --sarc-series-5: #e87ba4;   /* ACTN2   */
  --sarc-series-6: #008300;   /* barbed  */
}
@media (prefers-color-scheme: dark) {
  :root { --sarc-series-3: #199e70; --sarc-series-4: #c98500;
          --sarc-series-5: #d55181; --sarc-series-6: #008300; }
}
```

```javascript
color: 'var(--viz-series-3, var(--sarc-series-3))'
```

Those four hues are the dataviz reference palette and pass all six checks of
`validate_palette.js` on both the light and the dark chart surface. The light
surface returns a contrast WARN, which the always-present layer legend
discharges. Re-run the validator if any hue is changed.

## Safe inline data encoding

Pack numeric layers in Python exactly as contiguous little-endian/native
`float32` bytes and convert the result to an ASCII string:

```python
def encode(values: np.ndarray) -> str:
    packed = np.ascontiguousarray(values, dtype=np.float32)
    return base64.b64encode(packed.tobytes()).decode("ascii")
```

Inject the payload with `json.dumps`, not Python `repr`. Require exactly one
template marker before replacement.

Do **not** use browser `atob()` for these large inline payloads. The previous
`a1600-latest-snapshot.html` failure (`InvalidCharacterError: ... characters
outside of the Latin1 range`) came from sending a non-clean string to `atob`.
Use the manual decoder already present in
`a2400-latest-snapshots-template.html` and the fixed A1600 artifact:

```javascript
function decode(base64, Type) {
  const alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
  const clean = base64.replace(/=+$/, '');
  const byteLength = Math.floor(clean.length * 6 / 8);
  const bytes = new Uint8Array(byteLength);
  let buffer = 0, bits = 0, offset = 0;
  for (const character of clean) {
    const value = alphabet.indexOf(character);
    if (value < 0) throw new Error('Invalid embedded snapshot payload');
    buffer = (buffer << 6) | value;
    bits += 6;
    if (bits >= 8) {
      bits -= 8;
      bytes[offset++] = (buffer >> bits) & 255;
    }
  }
  return new Type(bytes.buffer);
}
```

Do not pass the payload through Unicode escaping, URL decoding, HTML entity
decoding, `decodeURIComponent`, or a text encoder/decoder.

## Building a viewer

For the existing generic M/Z plus snapshot page:

```bash
VIS=/home/users/xsli/.codex/visualizations/2026/08/26/01a03f41-9578-7e83-9ee5-f379f3493c26
python3 "$VIS/build-simulations.py" \
  --m320 /absolute/path/to/run-one.h5 \
  --m640 /absolute/path/to/run-two.h5 \
  --template "$VIS/simulations.template.html" \
  --output "$VIS/simulations.html"
```

For the compact, dependency-free comparison viewer, edit only the `RUNS`,
`TEMPLATE`, `OUTPUT`, and target-time declarations in a copy of
`build-a2400-latest-snapshots.py`, then run that copy with `python3`.

When Claude Code is not allowed to write under `.codex/visualizations`, place
the builder, template, and rendered HTML under a new directory such as
`analysis/visualizations/<study-name>/`; the viewer itself is unchanged.

## Validation checklist

Before reporting a visualization as complete:

1. Run `python3 -m py_compile <builder.py>`.
2. Build the HTML and confirm the builder reports nonzero layer counts.
3. Confirm the output contains no unreplaced template marker.
4. Confirm every embedded layer string matches `^[A-Za-z0-9+/]*={0,2}$`.
5. Confirm the template/output does not call `atob(`.
5b. Confirm no layer resolves to black: every `var(--viz-series-*)`,
   `var(--foreground)` and `var(--muted-foreground)` reference carries a
   fallback, and the artifact defines those fallbacks itself. Grep for
   `var(--viz-series-[0-9])` with no comma -- any hit is the black-filament bug.
6. Decode at least one layer in Python and verify its float count is divisible
   by 6 for segments or 3 for points, with all values finite.
7. Verify the displayed step equals the latest committed HDF5 step.
8. Open the HTML and test run selection, every layer toggle, drag rotation,
   wheel zoom, resizing, and readable light/dark theme colors.
9. Compare the visible counts and M/Z values with direct HDF5 calculations.
   The payload is injected as a bare JSON **list** at `const runs = /*__RUNS__*/`,
   and layers live at `run["snapshot"]["layers"]`. A validator that looks
   anywhere else finds nothing and passes vacuously; assert a nonzero run and
   layer count before trusting any check.
10. Return the absolute path to the rendered HTML and state which committed
    frame/step/time it shows.

## Agent instruction

When asked to visualize one or more sarcomere trajectories, follow this file as
an execution specification. Reuse the canonical extraction and viewer code;
do not replace it with a static image unless the user explicitly asks for one.
Inspect the actual HDF5 schema and committed frame count, generate the inline
HTML, validate it with the checklist above, and report the artifact path plus
the snapshot metadata. Preserve the layer style contract across comparisons.

