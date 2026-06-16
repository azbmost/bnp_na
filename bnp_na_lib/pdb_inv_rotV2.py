#!/usr/bin/env python3
"""
pdb_inv_rotV2.py

Apply coordinate inversion/reflection operations to all ATOM/HETATM records in
a PDB file.

This module is used by bnp_na to generate mirror-image nucleic-acid models
after DSSR align-to-Z and before final placement/orientation.

Input:
    1st command-line argument: input PDB filename
    2nd command-line argument: optional operation instruction

Output:
    A new PDB file with a suffix added before the extension.

Operation modes:
    i mode: point inversion followed by optional 180-degree rotation(s)
        i       : inversion only
        ix      : inversion + 180-degree rotation around x
        iy      : inversion + 180-degree rotation around y
        iz      : inversion + 180-degree rotation around z
        ixy     : inversion + rotations around x and y
        xyz     : same as ixyz, because i mode is default

    o mode: reflection across a coordinate plane
        oxy     : reflection across xy plane
        oyz     : reflection across yz plane
        oxz     : reflection across xz plane

Examples:
    python pdb_inv_rotV2.py model.pdb
        -> model_i.pdb

    python pdb_inv_rotV2.py model.pdb x
        -> model_i_x.pdb

    python pdb_inv_rotV2.py model.pdb ix
        -> model_i_x.pdb

    python pdb_inv_rotV2.py model.pdb oxy
        -> model_o_xy.pdb
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from edit_pdb_atom import file2rec, rec2file


VALID_AXES = set("xyz")
VALID_PLANES = {"xy", "yz", "xz"}
REFLECTION_TO_ROTATION_AXIS = {
    "xy": "z",
    "yz": "x",
    "xz": "y",
}


class InvRotError(Exception):
    """Raised when inversion/reflection of a PDB file fails."""

    def __init__(self, message: str, log_text: str = ""):
        super().__init__(message)
        self.log_text = log_text


def clean_instruction(text: Optional[str]) -> str:
    """Normalize the user-provided operation instruction."""
    if text is None:
        return ""
    return text.strip().lower().replace("_", "").replace("-", "").replace(" ", "")


def normalize_axes(axes: str) -> str:
    """Normalize rotation axes and ignore duplicate axes."""
    clean = []
    seen = set()

    for ch in axes:
        if ch not in VALID_AXES:
            raise ValueError(f"Invalid rotation axis {ch!r}. Allowed axes are x, y, and z.")
        if ch not in seen:
            clean.append(ch)
            seen.add(ch)

    return "".join(clean)


def parse_operation(instruction: Optional[str]) -> Tuple[str, str, str]:
    """Parse an inversion/reflection instruction.

    Returns:
        mode:
            "i" for inversion mode or "o" for reflection mode.
        label:
            File-name-safe operation label.
        rotation_axes:
            Axes for 180-degree rotation after inversion.
    """
    instr = clean_instruction(instruction)

    if instr == "":
        return "i", "i", ""

    if instr.startswith("o"):
        plane = instr[1:]
        if plane not in VALID_PLANES:
            raise ValueError(f"Invalid reflection instruction {instruction!r}. Use oxy, oyz, or oxz.")
        rotation_axes = REFLECTION_TO_ROTATION_AXIS[plane]
        return "o", f"o_{plane}", rotation_axes

    if instr.startswith("i"):
        axes = normalize_axes(instr[1:])
    else:
        axes = normalize_axes(instr)

    return "i", f"i_{axes}" if axes else "i", axes


def rotation_signs(rotation_axes: str) -> Tuple[int, int, int]:
    """Compute sign multipliers for a composition of 180-degree rotations."""
    sx, sy, sz = 1, 1, 1

    for axis in rotation_axes:
        if axis == "x":
            sy *= -1
            sz *= -1
        elif axis == "y":
            sx *= -1
            sz *= -1
        elif axis == "z":
            sx *= -1
            sy *= -1

    return sx, sy, sz


def final_coordinate_signs(rotation_axes: str) -> Tuple[int, int, int]:
    """Return signs after point inversion plus optional 180-degree rotations."""
    rx, ry, rz = rotation_signs(rotation_axes)
    return -rx, -ry, -rz


def make_outname(infile: Union[str, Path], label: str) -> str:
    """Create output filename by appending the operation label before the extension."""
    root, ext = os.path.splitext(str(infile))
    if not ext:
        ext = ".pdb"
    return f"{root}_{label}{ext}"


def update_record_xyz(rec, x: float, y: float, z: float) -> None:
    """Update atom coordinates while preserving PDB line formatting."""
    if hasattr(rec, "update_xyz"):
        rec.update_xyz(x, y, z)
    else:
        rec.x = x
        rec.y = y
        rec.z = z


def process_pdb(infile: Union[str, Path], outfile: Union[str, Path], rotation_axes: str) -> int:
    """Apply inversion plus optional rotation to all ATOM/HETATM records."""
    in_path = Path(infile)
    out_path = Path(outfile)
    rec_list = []

    with open(in_path, "r", encoding="utf-8", errors="ignore") as fin:
        file2rec(fin, rec_list)

    sx, sy, sz = final_coordinate_signs(rotation_axes)
    n_changed = 0

    for rec in rec_list:
        if getattr(rec, "recordName", "") in ("ATOM", "HETATM"):
            update_record_xyz(rec, float(rec.x) * sx, float(rec.y) * sy, float(rec.z) * sz)
            n_changed += 1

    with open(out_path, "w", encoding="utf-8") as fout:
        rec2file(rec_list, fout, reorder_serial=False)

    return n_changed


def describe_operation(mode: str, label: str, rotation_axes: str) -> str:
    """Return a human-readable explanation of an inversion/reflection operation."""
    sx, sy, sz = final_coordinate_signs(rotation_axes)
    lines = [f"Operation label: {label}"]

    if mode == "i":
        lines += [
            "Mode           : inversion",
            f"Rotation axes  : {rotation_axes if rotation_axes else '(none)'}",
        ]
    elif mode == "o":
        plane = label.replace("o_", "")
        lines += [
            "Mode           : reflection",
            f"Reflection     : across {plane} plane",
            f"Implemented as : inversion + 180-degree rotation around {rotation_axes} axis",
        ]

    lines.append(
        "Final transform: (x, y, z) -> ({}x, {}y, {}z)".format(
            "+" if sx > 0 else "-",
            "+" if sy > 0 else "-",
            "+" if sz > 0 else "-",
        )
    )
    return "\n".join(lines)


def apply_inv_rot_to_pdb(
    infile: Union[str, Path],
    outfile: Optional[Union[str, Path]] = None,
    instruction: Optional[str] = "",
) -> Dict[str, object]:
    """Apply an inversion/reflection instruction and return output metadata."""
    logs = []
    in_path = Path(infile)
    if not in_path.exists():
        raise InvRotError(f"Input PDB not found: {in_path}")

    try:
        mode, label, rotation_axes = parse_operation(instruction)
    except ValueError as exc:
        raise InvRotError(str(exc)) from exc

    out_path = Path(outfile) if outfile is not None else Path(make_outname(in_path, label))
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        n_changed = process_pdb(in_path, out_path, rotation_axes)
    except Exception as exc:
        raise InvRotError(f"PDB inv/rot operation failed: {exc}", "\n".join(logs)) from exc

    logs += [
        "=== L-form mirror operation (pdb_inv_rotV2.py) ===",
        f"Input : {in_path}",
        f"Output: {out_path}",
        describe_operation(mode, label, rotation_axes),
        f"Updated ATOM/HETATM coordinate records: {n_changed}",
    ]

    return {
        "pdb_in": in_path,
        "pdb_out": out_path,
        "mode": mode,
        "label": label,
        "rotation_axes": rotation_axes,
        "coordinate_signs": final_coordinate_signs(rotation_axes),
        "n_changed": n_changed,
        "log_text": "\n".join(logs),
    }


def print_operation_summary(mode: str, label: str, rotation_axes: str) -> None:
    """Print a short explanation of the operation."""
    for line in describe_operation(mode, label, rotation_axes).splitlines():
        print("[INFO]", line)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0)

    infile = sys.argv[1]
    instruction = sys.argv[2] if len(sys.argv) >= 3 else ""

    if not os.path.isfile(infile):
        print(f"[ERROR] Cannot find input file: {infile}")
        sys.exit(1)

    try:
        result = apply_inv_rot_to_pdb(infile, instruction=instruction)
    except InvRotError as exc:
        print("[ERROR]", exc)
        if exc.log_text:
            print(exc.log_text)
        sys.exit(1)

    print("[INFO] Input file     :", result["pdb_in"])
    print("[INFO] Output file    :", result["pdb_out"])
    print_operation_summary(str(result["mode"]), str(result["label"]), str(result["rotation_axes"]))
    print(f"[DONE] Updated {result['n_changed']} ATOM/HETATM coordinate records.")
    print("[DONE] Finished writing:", result["pdb_out"])


if __name__ == "__main__":
    main()
