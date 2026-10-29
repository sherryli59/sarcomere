# Newly mature ACTN2 links entering the excluded region

Source: original production trajectory
`$SCRATCH/sarcomere/gaphunt_20261023/data/r27zt0npg130_seed20260805.h5`.
These clips use the production run, not the diagnostic continuation.

Combined video: `newly_mature_entries.mp4`, 1600 x 900, H.264, 14.75 s.
Both selected molecules were born and matured after the 2.0 s resume.
Each was observed globally clear as status 2 before its first sampled entry
into the finite 130 nm cylinder of the highlighted myosin.

| Case | Link / actins | Myosin | Birth | Maturation | First sampled entry |
|---|---|---|---|---|---|
| Axial entry | 168764 / 1337, 1084 | 32, minus end | step 208820 | step 209057 | step 211000, 19.43 ms later |
| Lateral entry | 171917 / 484, 159 | 255, plus end | step 215131 | step 215133 | step 219500, 43.67 ms later |

The axial clip spans frames 418–436 (2.090–2.180 s).
The lateral clip spans frames 430–453 (2.150–2.265 s).

In the axial case, the complete connector is still 0.322 nm beyond the end
plane in the preceding saved pose; it then intersects the cylinder and
remains mature. In the lateral case, the minimum connector radius within
axial overlap falls from 154.458 nm at 2.190 s to 128.872 nm at 2.195 s,
while the overlap is already deep within the thick filament's axial span.

## Reading the videos

- Blue region/wireframe: finite 130 nm maturation exclusion cylinder.
- Gray body: myosin radius 50 nm.
- Teal/orange: the two tracked actins; dots mark material barbed ends.
- Purple connector: clear of the selected myosin cylinder.
- Red connector: intersects the cylinder. The title reports its actual status.
- X: minimum-radius connector point within the physical axial slab.

The 3D camera follows the selected myosin tip, with a smoothly transported
transverse basis. Filaments are clipped to the local view. The axial–radial
panel shows exact radial distance sqrt(y^2+z^2), rather than a single transverse
projection. The bottom panels show radial eligibility and axial clearance.
Maturation time comes from the exact event log; position at that exact step
was not saved unless it coincided with a sampled frame.

Saved poses are 5 ms apart. Each pose is held for slow playback, with extra
pauses at clip endpoints and first entry. No intermediate coordinates or
interpolated motion are generated. Future sampled trace segments are faint;
the moving cursor marks the displayed pose.

## Reproduction and validation

Script: `../../render_entry_movies.py`.
Rendering was performed in Slurm job 49981811 using the existing sarcomere
environment's matplotlib and ffmpeg. Selection requires at least two
clear mature poses and multiple excluded mature poses for the same myosin.
The script also verifies that the independent geometry agrees with these
selected states.

The per-case JSON files contain source path, size, modification time, source
binary hash, rendering-script hash, identities, event times, selected frame
range, local coordinates, and measured per-frame geometry. `selection.json`
records the candidate ranking; the selected clips favor entries within
100 ms of maturation. The JSON's `deepest_frame` is the selection search
result; `poster_frame` is constrained to the actual displayed clip.

The combined stream was checked with ffprobe: H.264, 1600 x 900, 177 frames,
14.75 seconds. Both generated posters were visually inspected.

Primary scratch outputs:
`$SCRATCH/sarcomere/gaphunt_20261023/torque_out/status2_gate_20261028/entry_movies_v2/`.
This directory is the durable copy of the completed outputs.

For further Sherlock work, read `/etc/agents/AGENTS.md` and relevant topic
files at session start and after context compaction.
