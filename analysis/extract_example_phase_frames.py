#!/usr/bin/env python3
"""Extract representative actin/myosin center positions for target phase values.

This script chooses, for each target Phi_3D value, the candidate frame with the
largest number of actins in cb_status == 2 among frames within the requested
absolute Phi_3D error tolerance. It then extracts that exact frame from the
corresponding H5 file.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np


DEFAULT_TARGET_PHI = (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
DEFAULT_MAX_ABS_ERROR = 0.02
DEFAULT_EXCLUDED_NX = (6, 7)
DEFAULT_BASE_DIRS = (
    Path("/scratch/users/xsli/sarcomere/data0416"),
    Path("/scratch/users/xsli/sarcomere/data0430"),
)
FILENAME_RE = re.compile(
    r"^A(?P<actins>\d+)_M(?P<myosins>\d+)_lt(?P<lt>.+)_r(?P<replicate>\d+)_L.*\.h5$"
)


@dataclass(frozen=True)
class Candidate:
    filename: str
    path: Path
    actins: int
    myosins: int
    ratio: float
    lt: str
    replicate: str
    frame: int
    phi: float
    theta_deg: float | None
    nx: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each requested Phi_3D value, pick the frame within tolerance "
            "with the most cb_status==2 actins and save myosin centers plus "
            "cb_status==2 actin centers."
        )
    )
    parser.add_argument(
        "--out_csv",
        "--csv",
        dest="out_csv",
        type=Path,
        default=Path("structural_phase_data.csv"),
        help="Path to structural_phase_data.csv. The --csv spelling is kept as an alias.",
    )
    parser.add_argument(
        "--per_frame_csv",
        type=Path,
        default=None,
        help="Path to sampled per-frame structural CSV; defaults to <out_csv stem>_per_frame.csv.",
    )
    parser.add_argument(
        "--base_dirs",
        "--data-root",
        dest="base_dirs",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Base directories to recursively search for H5 files. "
            "The --data-root spelling is kept as an alias."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/example_phase_frame_centers.npz"),
        help="Output .npz path for numeric center arrays.",
    )
    parser.add_argument(
        "--targets",
        type=float,
        nargs="+",
        default=list(DEFAULT_TARGET_PHI),
        help="Target Phi_3D values.",
    )
    parser.add_argument(
        "--max-abs-error",
        type=float,
        default=DEFAULT_MAX_ABS_ERROR,
        help="Maximum allowed absolute Phi_3D error for an eligible frame.",
    )
    parser.add_argument(
        "--nmyosins",
        type=int,
        default=1000,
        help="Myosin count used to compute actin:myosin ratio from A_val.",
    )
    parser.add_argument(
        "--exclude-nx",
        type=int,
        nargs="*",
        default=list(DEFAULT_EXCLUDED_NX),
        help=(
            "Exclude per-frame rows with these winning nx values. "
            "Defaults to 6 7. Use --exclude-nx with no values to disable."
        ),
    )
    parser.add_argument(
        "--frame",
        default="matched",
        help='Frame to extract: "matched" uses the closest per-frame CSV row; also accepts "last", "middle", or an integer index.',
    )
    args = parser.parse_args()
    if args.base_dirs is None:
        args.base_dirs = list(DEFAULT_BASE_DIRS)
    return args


def default_per_frame_csv_path(out_csv: Path) -> Path:
    return out_csv.with_name(f"{out_csv.stem}_per_frame{out_csv.suffix}")


def resolve_csv_path(path: Path) -> Path:
    if path.exists():
        return path
    if not path.is_absolute():
        script_dir = Path(__file__).resolve().parent
        candidates = [
            script_dir / path,
            script_dir.parent / path,
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return path


def discover_h5_files(roots: Iterable[Path]) -> dict[str, Path]:
    files: dict[str, Path] = {}
    for root in roots:
        root = root.expanduser()
        if not root.exists():
            continue
        for path in root.resolve().rglob("A*_M*_lt*_r*_L*.h5"):
            # Keep the first matching basename. Sorting roots/files makes this deterministic.
            files.setdefault(path.name, path)
    return dict(sorted(files.items()))


def resolve_h5_path(row: dict[str, str], csv_path: Path, h5_files: dict[str, Path]) -> Path | None:
    h5_path = row.get("H5_Path") or row.get("h5_path") or row.get("Path")
    if h5_path:
        path = Path(h5_path).expanduser()
        if not path.is_absolute():
            path = csv_path.parent / path
        if path.exists():
            return path.resolve()
    return h5_files.get(row["Filename"])


def parse_filename(filename: str) -> tuple[int | None, int | None, str | None, str | None]:
    match = FILENAME_RE.match(filename)
    if not match:
        return None, None, None, None
    return (
        int(match.group("actins")),
        int(match.group("myosins")),
        match.group("lt"),
        match.group("replicate"),
    )


def load_candidates(
    csv_path: Path,
    h5_files: dict[str, Path],
    nmyosins: int,
    excluded_nx: set[int],
) -> list[Candidate]:
    candidates: list[Candidate] = []
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"Filename", "A_val", "lt_str", "Frame", "Phi_3D"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{csv_path} is missing columns: {sorted(missing)}")

        for row in reader:
            filename = row["Filename"]
            path = resolve_h5_path(row, csv_path, h5_files)
            if path is None:
                continue

            parsed_actins, parsed_myosins, parsed_lt, replicate = parse_filename(path.name)
            if parsed_actins is None:
                parsed_actins, parsed_myosins, parsed_lt, replicate = parse_filename(filename)
            actins = int(row["A_val"]) if row.get("A_val") else parsed_actins
            myosins = parsed_myosins or nmyosins
            if actins is None or myosins is None:
                continue

            theta = None
            if row.get("Theta_deg"):
                theta = float(row["Theta_deg"])
            nx = int(row["nx"]) if row.get("nx") else None
            if nx in excluded_nx:
                continue

            candidates.append(
                Candidate(
                    filename=filename,
                    path=path,
                    actins=int(actins),
                    myosins=int(myosins),
                    ratio=float(actins) / float(nmyosins),
                    lt=row.get("lt_str") or parsed_lt or "",
                    replicate=replicate or "",
                    frame=int(row["Frame"]),
                    phi=float(row["Phi_3D"]),
                    theta_deg=theta,
                    nx=nx,
                )
            )
    return candidates


def count_actin_cb2(path: Path, frame: int) -> int:
    with h5py.File(path, "r") as fh:
        if "/actin/cb_status" not in fh:
            raise KeyError(f"{path} is missing required dataset /actin/cb_status")
        n_frames = fh["/actin/cb_status"].shape[0]
        resolved_frame = resolve_frame(str(frame), n_frames)
        cb_status = np.asarray(fh["/actin/cb_status"][resolved_frame]).reshape(-1)
    return int(np.count_nonzero(cb_status == 2))


def choose_candidates(
    candidates: list[Candidate],
    targets: Iterable[float],
    max_abs_error: float,
) -> list[tuple[float, Candidate]]:
    selected: list[tuple[float, Candidate]] = []
    cb2_count_cache: dict[tuple[Path, int], int] = {}
    for target in targets:
        ranked = []
        for candidate in candidates:
            phi_abs_error = abs(candidate.phi - target)
            if phi_abs_error > max_abs_error:
                continue
            cache_key = (candidate.path, candidate.frame)
            if cache_key not in cb2_count_cache:
                cb2_count_cache[cache_key] = count_actin_cb2(candidate.path, candidate.frame)
            ranked.append(
                (
                    -cb2_count_cache[cache_key],
                    phi_abs_error,
                    candidate.filename,
                    candidate.frame,
                    candidate,
                )
            )
        if not ranked:
            raise RuntimeError(
                f"No H5-backed CSV candidates matched target Phi_3D={target:g} "
                f"within max absolute error {max_abs_error:g}."
            )
        ranked.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
        selected.append((target, ranked[0][4]))
    return selected


def resolve_frame(frame_spec: str, n_frames: int) -> int:
    if n_frames <= 0:
        raise ValueError("Cannot extract from an H5 file with zero frames.")
    if frame_spec == "last":
        return n_frames - 1
    if frame_spec == "middle":
        return n_frames // 2
    frame = int(frame_spec)
    if frame < 0:
        frame = n_frames + frame
    if not 0 <= frame < n_frames:
        raise IndexError(f"Frame {frame_spec!r} is outside available range 0..{n_frames - 1}.")
    return frame


def extract_centers(path: Path, frame_spec: str | int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as fh:
        for dataset in ("/actin/center", "/actin/cb_status", "/myosin/center"):
            if dataset not in fh:
                raise KeyError(f"{path} is missing required dataset {dataset}")

        n_frames = min(
            fh["/actin/center"].shape[0],
            fh["/actin/cb_status"].shape[0],
            fh["/myosin/center"].shape[0],
        )
        frame = resolve_frame(str(frame_spec), n_frames)

        actin_centers = np.asarray(fh["/actin/center"][frame], dtype=np.float64)
        cb_status = np.asarray(fh["/actin/cb_status"][frame]).reshape(-1)
        myosin_centers = np.asarray(fh["/myosin/center"][frame], dtype=np.float64)

    cb2_indices = np.flatnonzero(cb_status == 2).astype(np.int64)
    return frame, actin_centers[cb2_indices], cb2_indices, myosin_centers


def write_outputs(
    output_path: Path,
    selected: list[tuple[float, Candidate]],
    frame_spec: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    metadata = {
        "frame_spec": frame_spec,
        "samples": [],
    }

    for idx, (target_phi, candidate) in enumerate(selected):
        frame_to_extract = candidate.frame if frame_spec == "matched" else frame_spec
        frame, actin_cb2_centers, actin_cb2_indices, myosin_centers = extract_centers(
            candidate.path, frame_to_extract
        )
        prefix = f"sample_{idx:02d}"
        arrays[f"{prefix}_actin_cb2_centers"] = actin_cb2_centers
        arrays[f"{prefix}_actin_cb2_indices"] = actin_cb2_indices
        arrays[f"{prefix}_myosin_centers"] = myosin_centers

        metadata["samples"].append(
            {
                "sample_key": prefix,
                "target_phi": float(target_phi),
                "matched_phi": float(candidate.phi),
                "phi_abs_error": float(abs(candidate.phi - target_phi)),
                "matched_frame_from_csv": int(candidate.frame),
                "filename": candidate.filename,
                "path": str(candidate.path),
                "actins": candidate.actins,
                "myosins": candidate.myosins,
                "actin_myosin_ratio": candidate.ratio,
                "lt": candidate.lt,
                "replicate": candidate.replicate,
                "theta_deg": candidate.theta_deg,
                "nx": candidate.nx,
                "frame_index": int(frame),
                "n_actin_cb2": int(actin_cb2_centers.shape[0]),
                "n_myosin": int(myosin_centers.shape[0]),
            }
        )

    np.savez_compressed(output_path, **arrays)
    metadata_path = output_path.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print(f"saved arrays: {output_path}")
    print(f"saved metadata: {metadata_path}")
    for sample in metadata["samples"]:
        print(
            "{sample_key}: target_phi={target_phi:.3g}, matched_phi={matched_phi:.6g}, "
            "ratio={actin_myosin_ratio:.3g}, frame={frame_index}, "
            "actin_cb2={n_actin_cb2}, myosin={n_myosin}, file={filename}".format(**sample)
        )


def main() -> None:
    args = parse_args()
    out_csv_path = resolve_csv_path(args.out_csv)
    per_frame_csv = args.per_frame_csv or default_per_frame_csv_path(out_csv_path)
    csv_path = resolve_csv_path(per_frame_csv)
    h5_files = discover_h5_files(args.base_dirs)

    candidates = load_candidates(csv_path, h5_files, args.nmyosins, set(args.exclude_nx))
    if not candidates:
        raise RuntimeError(
            f"No per-frame CSV rows in {csv_path} had usable H5 paths. Regenerate the per-frame CSV with H5_Path or pass --base_dirs for fallback discovery."
        )

    selected = choose_candidates(
        candidates,
        args.targets,
        args.max_abs_error,
    )
    write_outputs(args.output, selected, args.frame)


if __name__ == "__main__":
    main()
