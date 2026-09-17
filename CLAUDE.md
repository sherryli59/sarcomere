# Repository instructions for Claude Code

## Sarcomere visualization

When asked to visualize a trajectory, latest snapshot, simulation comparison,
or M/Z-order trend, read and follow `VISUALIZATION_PIPELINE.md` as the execution
specification.

Use the established interactive inline-HTML pipeline. Preserve its committed-
frame handling, periodic clipping, species color and thickness contract, layer
toggles, M/Z definitions, and manual ASCII base64 decoder. Do not use `atob()`
and do not substitute a static image unless the user explicitly requests one.

Generate and validate the artifact, then report its absolute path and the
committed frame, step, and simulation time shown.
