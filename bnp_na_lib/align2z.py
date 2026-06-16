#!/usr/bin/env python3
"""Align a nucleic-acid PDB helix to +Z using DSSR's reported helix axis.

This module intentionally depends on ``x3dna-dssr --more`` for axis extraction,
because DSSR reports the helix ``point-one`` and ``point-two`` endpoints more
accurately than a generic PCA estimate for the start/end locations.

The output PDB is translated so that DSSR ``point-one`` is at the origin and
rotated so that the vector ``point-two - point-one`` points along +Z.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from edit_pdb_atom import file2rec, rec2file
from geometry_utils import GeometryError, rotation_matrix_from_to


class Align2ZError(Exception):
    """Raised when DSSR axis extraction or alignment fails."""


def which_or(path_guess: str) -> Optional[str]:
    hit = shutil.which(path_guess)
    if hit:
        return hit
    if Path(path_guess).is_file():
        return path_guess
    return None


def _path_arg(path: Union[str, Path], cwd: Union[str, Path]) -> str:
    p = Path(path)
    cwd_path = Path(cwd).resolve()
    try:
        if p.resolve().parent == cwd_path:
            return p.name
    except Exception:
        pass
    return str(p)


def command_to_text(cmd: Sequence[str]) -> str:
    return " ".join(str(part) for part in cmd)


def run_dssr_more_axis(
    pdb_path: Union[str, Path],
    out_txt: Union[str, Path],
    cwd: Union[str, Path],
) -> Tuple[bool, str, List[str]]:
    """Run DSSR ``--more`` to obtain the global helix axis.

    Equivalent command:
        x3dna-dssr -i=<pdb_path> -o=<out_txt> --more
    """
    exe = which_or("x3dna-dssr") or "/usr/local/bin/x3dna-dssr"
    if not which_or(exe):
        return False, "x3dna-dssr not found in PATH or at /usr/local/bin/x3dna-dssr", []

    cmd = [
        exe,
        f"-i={_path_arg(pdb_path, cwd)}",
        f"-o={_path_arg(out_txt, cwd)}",
        "--more",
    ]
    try:
        cp = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
        output = (cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")
        return cp.returncode == 0, output.strip(), cmd
    except Exception as exc:
        return False, str(exc), cmd


def parse_first_axis_points(out_txt_path: Union[str, Path]) -> Tuple[List[float], List[float]]:
    """Read the first DSSR ``point-one`` and ``point-two`` from ``--more`` output."""
    p1 = None
    p2 = None
    with open(out_txt_path, "r", encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            stripped = line.strip()
            if stripped.startswith("point-one:"):
                parts = stripped.split("point-one:", 1)[1].strip().split()
                if len(parts) >= 3:
                    p1 = [float(parts[0]), float(parts[1]), float(parts[2])]
            elif stripped.startswith("point-two:"):
                parts = stripped.split("point-two:", 1)[1].strip().split()
                if len(parts) >= 3:
                    p2 = [float(parts[0]), float(parts[1]), float(parts[2])]
            if p1 is not None and p2 is not None:
                break

    if p1 is None or p2 is None:
        raise Align2ZError("Failed to find 'point-one'/'point-two' in DSSR --more output.")
    return p1, p2


def _load_pdb_records(pdb_path: Union[str, Path]):
    records = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fin:
        file2rec(fin, records)
    return records


def _is_coordinate_record(record) -> bool:
    return all(hasattr(record, attr) for attr in ("x", "y", "z"))


def _record_coord(record) -> np.ndarray:
    return np.array([float(record.x), float(record.y), float(record.z)], dtype=float)


def _apply_alignment(
    inp_pdb: Union[str, Path],
    out_pdb: Union[str, Path],
    point_one: Iterable[float],
    point_two: Iterable[float],
) -> np.ndarray:
    p1 = np.asarray(list(point_one), dtype=float)
    p2 = np.asarray(list(point_two), dtype=float)
    axis_vec = p2 - p1
    if np.linalg.norm(axis_vec) <= 1e-8:
        raise Align2ZError("DSSR point-one and point-two define a near-zero axis length.")

    try:
        rotation = rotation_matrix_from_to(axis_vec, np.array([0.0, 0.0, 1.0], dtype=float))
    except GeometryError as exc:
        raise Align2ZError(str(exc)) from exc

    records = _load_pdb_records(inp_pdb)
    for rec in records:
        if not _is_coordinate_record(rec):
            continue
        coord = _record_coord(rec)
        coord = rotation @ (coord - p1)
        rec.update_xyz(*coord)

    with open(out_pdb, "w", encoding="utf-8") as fout:
        rec2file(records, fout, reorder_serial=False)

    return rotation


def align_pdb_to_Z(
    inp_pdb: str,
    out_pdb: Optional[str] = None,
    start_point: Optional[Iterable[float]] = None,
    end_point: Optional[Iterable[float]] = None,
    dssr_out: Optional[str] = None,
    cwd: Optional[Union[str, Path]] = None,
    atom_selection: Optional[str] = None,
) -> Dict[str, object]:
    """Align a PDB to +Z and return metadata.

    By default, this function runs ``x3dna-dssr -i=<input> -o=<input>-out --more``,
    parses DSSR ``point-one``/``point-two``, and aligns using those endpoints.

    ``start_point``/``end_point`` are retained only for backward-compatible manual
    use. When both are supplied, DSSR is not called. ``atom_selection`` is accepted
    as a no-op compatibility argument from the earlier PCA implementation.
    """
    inp_path = Path(inp_pdb).expanduser()
    if not inp_path.is_absolute():
        inp_path = inp_path.resolve()
    if not inp_path.exists():
        raise Align2ZError(f"Input PDB not found: {inp_path}")

    if out_pdb is None:
        out_path = inp_path.with_name(inp_path.stem + "_aligned2Z" + inp_path.suffix)
    else:
        out_path = Path(out_pdb)

    workdir = Path(cwd) if cwd is not None else inp_path.resolve().parent
    workdir.mkdir(parents=True, exist_ok=True)

    cmd: List[str] = []
    dssr_log = ""
    method = "provided_points"
    dssr_out_path: Optional[Path] = None

    if start_point is not None and end_point is not None:
        p1 = [float(x) for x in start_point]
        p2 = [float(x) for x in end_point]
    else:
        raw_dssr_out = Path(dssr_out) if dssr_out is not None else Path(str(inp_path) + "-out")
        dssr_out_path = raw_dssr_out if raw_dssr_out.is_absolute() else (workdir / raw_dssr_out)
        ok_more, dssr_log, cmd = run_dssr_more_axis(inp_path, dssr_out_path, cwd=workdir)
        if not ok_more or not dssr_out_path.exists():
            raise Align2ZError(
                "DSSR --more failed during align-to-Z axis extraction.\n"
                f"Command: {command_to_text(cmd)}\n"
                f"Expected output file: {dssr_out_path}\n"
                f"Output:\n{dssr_log}"
            )
        p1, p2 = parse_first_axis_points(dssr_out_path)
        method = "dssr_more"

    rotation = _apply_alignment(inp_path, out_path, p1, p2)

    return {
        "input_pdb": str(inp_path),
        "output_pdb": str(out_path),
        "method": method,
        "dssr_more_out": str(dssr_out_path) if dssr_out_path is not None else None,
        "dssr_command": command_to_text(cmd) if cmd else None,
        "dssr_output": dssr_log,
        "point_one": [float(x) for x in p1],
        "point_two": [float(x) for x in p2],
        "rotation_matrix": rotation.tolist(),
        "atom_selection": atom_selection,
    }


def format_alignment_report(result: Dict[str, object]) -> str:
    lines = [
        "=== Align to +Z (DSSR --more axis) ===",
        f"Input : {result.get('input_pdb')}",
        f"Output: {result.get('output_pdb')}",
        f"Method: {result.get('method')}",
        f"DSSR --more output: {result.get('dssr_more_out')}",
    ]
    if result.get("dssr_command"):
        lines.append(f"DSSR command: {result.get('dssr_command')}")
    lines += [
        f"point-one: {result.get('point_one')}",
        f"point-two: {result.get('point_two')}",
    ]
    dssr_output = str(result.get("dssr_output") or "")
    if dssr_output:
        lines.append(f"DSSR output:\n{dssr_output}")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align a PDB helix to +Z using DSSR --more axis endpoints.")
    parser.add_argument("input_pdb", help="Input PDB file.")
    parser.add_argument("-o", "--output", dest="output_pdb", default=None, help="Output PDB file.")
    parser.add_argument(
        "--dssr-output",
        dest="dssr_out",
        default=None,
        help="Optional path for DSSR --more output. Default: <input_pdb>-out.",
    )
    parser.add_argument(
        "--cwd",
        dest="cwd",
        default=None,
        help="Working directory for DSSR. Default: input PDB folder.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = align_pdb_to_Z(
        args.input_pdb,
        out_pdb=args.output_pdb,
        dssr_out=args.dssr_out,
        cwd=args.cwd,
    )
    print(format_alignment_report(result))


if __name__ == "__main__":
    main()
