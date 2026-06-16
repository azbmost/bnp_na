"""Placement/orientation after a helix has been aligned to +Z."""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np

from edit_pdb_atom import file2rec, rec2file
from geometry_utils import rotation_matrix_y, rotation_matrix_z


class PlacerError(Exception):
    """Raised when placement/orientation fails."""

    def __init__(self, message: str, log_text: str = ""):
        super().__init__(message)
        self.log_text = log_text


def place_after_Z(
    pdb_in: str,
    pdb_out: str,
    roll_deg: float,
    phi_deg: float,
    theta_deg: float,
    tx: float,
    ty: float,
    tz: float,
    delta_z: float = 0.0,
    overwrite: bool = True,
) -> Dict[str, object]:
    """Orient/place a helix already aligned with start at origin and axis +Z.

    Transform order:
        pre-shift along +Z by delta_z -> roll about local +Z -> Ry(phi) ->
        Rz(theta) -> final translation.
    """
    logs = []
    in_path = Path(pdb_in)
    out_path = Path(pdb_out)

    if not in_path.exists():
        raise PlacerError(f"Input PDB not found: {in_path}")

    values = {
        "roll": roll_deg,
        "phi": phi_deg,
        "theta": theta_deg,
        "tx": tx,
        "ty": ty,
        "tz": tz,
        "delta_z": delta_z,
    }
    numeric = {}
    for name, value in values.items():
        try:
            numeric[name] = float(value)
            if not np.isfinite(numeric[name]):
                raise ValueError
        except Exception as exc:
            raise PlacerError(f"Invalid numeric parameter: {name}") from exc

    if out_path.exists() and overwrite:
        out_path.unlink()

    rz_roll = rotation_matrix_z(numeric["roll"])
    ry_phi = rotation_matrix_y(numeric["phi"])
    rz_theta = rotation_matrix_z(numeric["theta"])
    r_swing = rz_theta @ ry_phi
    rotation = r_swing @ rz_roll
    final_translation = np.array([numeric["tx"], numeric["ty"], numeric["tz"]], dtype=float)
    pre_translation = np.array([0.0, 0.0, numeric["delta_z"]], dtype=float)

    records = []
    with open(in_path, "r", encoding="utf-8", errors="ignore") as fin:
        file2rec(fin, records)

    for atom in records:
        if not all(hasattr(atom, attr) for attr in ("x", "y", "z")):
            continue
        coord = np.array([float(atom.x), float(atom.y), float(atom.z)], dtype=float)
        coord = coord + pre_translation
        coord = rotation @ coord
        coord = coord + final_translation
        atom.update_xyz(*coord)

    with open(out_path, "w", encoding="utf-8") as fout:
        rec2file(records, fout, reorder_serial=False)

    logs += [
        "=== Orient/Place (after +Z alignment) ===",
        f"Input  : {in_path}",
        f"Output : {out_path}",
        f"delta_z: {numeric['delta_z']} A",
        f"roll/phi/theta (deg): {numeric['roll']}, {numeric['phi']}, {numeric['theta']}",
        f"translate (A): ({numeric['tx']}, {numeric['ty']}, {numeric['tz']})",
    ]

    return {
        "pdb_in": in_path,
        "pdb_out": out_path,
        "roll_deg": numeric["roll"],
        "phi_deg": numeric["phi"],
        "theta_deg": numeric["theta"],
        "translate": (numeric["tx"], numeric["ty"], numeric["tz"]),
        "delta_z": numeric["delta_z"],
        "R": rotation.tolist(),
        "log_text": "\n".join(logs),
    }
