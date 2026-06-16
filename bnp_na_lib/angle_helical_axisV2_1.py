#!/usr/bin/env python3
# angle_helical_axisV2_1.py
#
# Fit (or accept a user-defined) helical axis and compute radial vectors from that axis
# to two points. Points can be given as PDB atoms (CHAIN:RESSEQ:ATOM) or explicit xyz.
# Writes a Chimera/ChimeraX .bild file that draws the axis and the two radial vectors.
# The axis drawing margin can be adjusted from the GUI or CLI.
#
# Example (fit axis from PDB):
#   python angle_helical_axisV2_1.py -i helix.pdb --point1 A:5:C1' --point2 B:18:C1' \
#       -o helix_axis_vectors.bild
#
# Example (custom axis, no PDB needed):
#   python angle_helical_axisV2_1.py --axis-point "0 0 0" --axis-vector "0 0 1" \
#       --point1 "1 0 0" --point2 "0 1 0" -o custom_axis_vectors.bild

"""angle_helical_axisV2_1.py

Identify an approximate nucleic-acid helical axis from a PDB (PCA/SVD fit), OR use a
user-supplied axis (point + direction vector). Then compute the perpendicular (radial)
vectors from that axis to two user-specified points.

Reports:
  * the helical axis (point on axis + unit direction)
  * distance of each point to the axis
  * unit radial vectors for each point
  * angle between the two radial vectors (degrees)

Points can be specified either as:
  * explicit coordinates:  "x,y,z"  or "x y z"
  * an atom in a PDB:      "CHAIN:RESSEQ:ATOM" (also accepts / , @ and spaces)

Axis options:
  * Default: fit axis from PDB using --axis-atoms (default C1')
  * Custom axis: provide BOTH --axis-point and --axis-vector (PDB not required unless
    you specify points by atom).

The script also writes a Chimera/ChimeraX .bild file that draws:
  * the axis as an arrow
  * the two radial vectors as arrows from the axis to each point
  * explanatory .comment lines before each drawn object

The displayed axis length can be extended or shortened with an axis drawing margin.

GUI mode:
  * If run with no arguments, or with --gui, a small Tkinter GUI is shown.
  * In the GUI, Point 1 and Point 2 each have a drop-down to choose XYZ vs Atom.

Dependencies:
  * numpy
  * tkinter (optional; only for GUI)

"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# -----------------------------
# PDB parsing (minimal, standalone)
# -----------------------------


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resName: str
    chainID: str
    resSeq: int
    x: float
    y: float
    z: float

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


def load_pdb_atoms(pdb_path: Path) -> List[Atom]:
    """Load ATOM/HETATM records from a PDB (or PDB-like .txt)."""
    atoms: List[Atom] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resName = line[17:20].strip()
                chainID = line[21].strip() or " "
                resSeq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                # Skip malformed lines
                continue
            atoms.append(Atom(serial, name, resName, chainID, resSeq, x, y, z))

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return atoms


def atom_aliases(atom_name: str) -> List[str]:
    """Return acceptable aliases for a PDB atom name.

    Supports prime/star equivalence: C1' <-> C1*
    """
    nm = atom_name.strip()
    if nm.endswith("'"):
        return [nm, nm[:-1] + "*"]
    if nm.endswith("*"):
        return [nm, nm[:-1] + "'"]
    return [nm]


def parse_axis_atoms(text: str) -> List[str]:
    """Parse a user string into a list of atom names (comma/space separated)."""
    if not text:
        return ["C1'"]
    parts = re.split(r"[\s,]+", text.strip())
    return [p for p in parts if p]


def build_atom_index(atoms: Sequence[Atom]) -> Dict[Tuple[str, int, str], Atom]:
    """Index atoms by (chainID, resSeq, atomName). Keep the first occurrence."""
    idx: Dict[Tuple[str, int, str], Atom] = {}
    for a in atoms:
        idx.setdefault((a.chainID, a.resSeq, a.name), a)
    return idx


# -----------------------------
# Geometry: axis fit, projection, angles
# -----------------------------


@dataclass
class Axis:
    point: np.ndarray  # a point on the axis
    direction: np.ndarray  # unit direction vector
    t_min: float
    t_max: float

    def endpoints(self, margin: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
        p1 = self.point + (self.t_min - margin) * self.direction
        p2 = self.point + (self.t_max + margin) * self.direction
        return p1, p2


def _normalize(v: np.ndarray, *, eps: float = 1e-12) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < eps:
        raise ValueError("Cannot normalize near-zero vector")
    return v / n


def fit_axis_pca(
    atoms: Sequence[Atom],
    axis_atom_names: Sequence[str],
    *,
    prefer_positive_correlation: bool = True,
) -> Axis:
    """Fit a best-fit line (axis) by PCA/SVD of selected atom coordinates."""
    alias_set = set()
    for nm in axis_atom_names:
        alias_set.update(atom_aliases(nm))

    pts: List[np.ndarray] = []
    chain_resseq: List[Tuple[str, int]] = []

    for a in atoms:
        if a.name in alias_set:
            pts.append(a.xyz)
            chain_resseq.append((a.chainID, a.resSeq))

    if len(pts) < 3:
        raise ValueError(
            f"Need at least 3 atoms to fit axis; found {len(pts)} atoms matching {sorted(alias_set)}"
        )

    X = np.stack(pts, axis=0)
    mean = X.mean(axis=0)
    Xc = X - mean

    # Principal component via SVD
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    direction = _normalize(vt[0])

    # Choose sign: try to make projection correlate positively with residue numbering
    if prefer_positive_correlation and len(chain_resseq) >= 6:
        counts: Dict[str, int] = {}
        for ch, _rs in chain_resseq:
            counts[ch] = counts.get(ch, 0) + 1
        major_chain = max(counts.items(), key=lambda kv: kv[1])[0]

        resSeqs: List[int] = []
        ts: List[float] = []
        for p, (ch, rs) in zip(pts, chain_resseq):
            if ch != major_chain:
                continue
            resSeqs.append(rs)
            ts.append(float(np.dot(p - mean, direction)))

        if len(resSeqs) >= 3:
            corr = float(np.corrcoef(resSeqs, ts)[0, 1])
            if not math.isnan(corr) and corr < 0:
                direction = -direction

    t_vals = [float(np.dot(p - mean, direction)) for p in pts]
    return Axis(point=mean, direction=direction, t_min=float(min(t_vals)), t_max=float(max(t_vals)))


def axis_from_point_and_vector(
    axis_point: np.ndarray,
    axis_vector: np.ndarray,
    ref_points: Sequence[np.ndarray],
) -> Axis:
    """Create an Axis from a user-supplied point and direction vector.

    t_min/t_max are chosen to span `ref_points` when projected onto the axis.
    """
    p0 = np.array(axis_point, dtype=float)
    u = _normalize(np.array(axis_vector, dtype=float))

    t_vals: List[float] = []
    for p in ref_points:
        t_vals.append(float(np.dot(np.array(p, dtype=float) - p0, u)))

    if not t_vals:
        # Fallback span if no reference points available
        t_min, t_max = -10.0, 10.0
    else:
        t_min, t_max = float(min(t_vals)), float(max(t_vals))
        if abs(t_max - t_min) < 1e-6:
            t_min -= 10.0
            t_max += 10.0

    return Axis(point=p0, direction=u, t_min=t_min, t_max=t_max)


def project_point_to_axis(p: np.ndarray, axis: Axis) -> Tuple[np.ndarray, float]:
    """Return (projection, t) of point p onto axis."""
    t = float(np.dot(p - axis.point, axis.direction))
    proj = axis.point + t * axis.direction
    return proj, t


@dataclass
class RadialResult:
    spec: str
    point: np.ndarray
    proj: np.ndarray
    vector: np.ndarray
    distance: float
    unit: Optional[np.ndarray]


def radial_vector_to_axis(p: np.ndarray, axis: Axis, spec: str = "") -> RadialResult:
    proj, _t = project_point_to_axis(p, axis)
    v = p - proj
    d = float(np.linalg.norm(v))
    u = None if d < 1e-12 else v / d
    return RadialResult(spec=spec, point=p, proj=proj, vector=v, distance=d, unit=u)


def angle_between_unit_vectors(u1: np.ndarray, u2: np.ndarray) -> float:
    dot = float(np.dot(u1, u2))
    dot = max(-1.0, min(1.0, dot))
    return float(math.degrees(math.acos(dot)))


def signed_angle_around_axis(u1: np.ndarray, u2: np.ndarray, axis_dir: np.ndarray) -> float:
    """Signed angle from u1 to u2 around axis_dir, in (-180, 180]."""
    x = float(np.dot(u1, u2))
    y = float(np.dot(axis_dir, np.cross(u1, u2)))
    ang = math.degrees(math.atan2(y, x))
    if ang <= -180:
        ang += 360
    elif ang > 180:
        ang -= 360
    return float(ang)


# -----------------------------
# Parsing helpers (xyz and point specs)
# -----------------------------


def _try_parse_xyz(spec: str) -> Optional[np.ndarray]:
    txt = spec.strip()
    txt2 = re.sub(r"[\s,]+", " ", txt)
    parts = [p for p in txt2.split(" ") if p]
    if len(parts) != 3:
        return None
    try:
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
    except Exception:
        return None
    return np.array([x, y, z], dtype=float)


def parse_xyz_triplet(spec: str, *, what: str) -> np.ndarray:
    xyz = _try_parse_xyz(spec)
    if xyz is None:
        raise ValueError(f"Could not parse {what} as xyz triplet: {spec!r} (expected 'x y z' or 'x,y,z')")
    return xyz


def parse_point_spec(
    spec: str,
    atom_index: Optional[Dict[Tuple[str, int, str], Atom]],
) -> np.ndarray:
    """Parse a point spec into xyz.

    Accepted forms:
      - "x,y,z" or "x y z"
      - "CHAIN:RESSEQ:ATOM" (also accepts / , @ or spaces)

    CHAIN can be omitted as "RESSEQ:ATOM" (treated as any chain, must be unique)
    or explicitly "*:RESSEQ:ATOM".

    If an atom spec is used, `atom_index` must be provided.
    """
    xyz = _try_parse_xyz(spec)
    if xyz is not None:
        return xyz

    if atom_index is None:
        raise ValueError(
            f"Point spec {spec!r} does not look like xyz, but no PDB atom index is available. "
            "Provide --pdb (or use xyz points)."
        )

    toks = [t for t in re.split(r"[:/@,\s]+", spec.strip()) if t]
    if len(toks) == 3:
        chain, resseq_s, atom = toks
    elif len(toks) == 2:
        chain, resseq_s, atom = "*", toks[0], toks[1]
    else:
        raise ValueError(
            f"Could not parse point spec {spec!r}. Expected 'x,y,z' or 'CHAIN:RESSEQ:ATOM'."
        )

    try:
        resseq = int(resseq_s)
    except Exception as e:
        raise ValueError(f"Bad residue number in {spec!r}: {resseq_s}") from e

    atom_alias = atom_aliases(atom)

    if chain == "*":
        matches: List[Atom] = []
        for (ch, rs, an), a in atom_index.items():
            if rs != resseq:
                continue
            if an in atom_alias:
                matches.append(a)

        if len(matches) == 0:
            raise ValueError(
                f"No atom found for spec {spec!r} (searched any chain, resSeq={resseq}, atom in {atom_alias})."
            )
        if len(matches) > 1:
            chains = sorted({m.chainID for m in matches})
            raise ValueError(
                f"Ambiguous spec {spec!r}: found {len(matches)} matches in chains {chains}. "
                "Specify CHAIN explicitly."
            )
        return matches[0].xyz

    chainID = chain.strip() or " "
    if len(chainID) != 1:
        chainID = chainID[0]

    for an in atom_alias:
        a = atom_index.get((chainID, resseq, an))
        if a is not None:
            return a.xyz

    raise ValueError(
        f"No atom found for spec {spec!r} (chain={chainID}, resSeq={resseq}, atom in {atom_alias})."
    )


# -----------------------------
# BILD output (Chimera/ChimeraX)
# -----------------------------


def write_bild(
    out_path: Path,
    axis: Axis,
    r1: RadialResult,
    r2: RadialResult,
    *,
    axis_color: Tuple[float, float, float] = (0.2, 0.2, 1.0),
    v1_color: Tuple[float, float, float] = (1.0, 0.2, 0.2),
    v2_color: Tuple[float, float, float] = (0.2, 0.8, 0.2),
    axis_radius: float = 1.0,
    vector_radius: float = 1.0,
    sphere_radius: float = 1.25,
    arrow_head_radius_scale: float = 2.5,
    arrow_shaft_fraction: float = 0.85,
    axis_margin: float = 5.0,
) -> None:
    """Write a Chimera/ChimeraX .bild file.

    Uses:
      * .comment for comments (more compatible with UCSF Chimera than '#')
      * .arrow for axis and vectors (user preference)

    .arrow syntax (Chimera/ChimeraX):
      .arrow x1 y1 z1 x2 y2 z2 [r1 [r2 [rho]]]
      where r1=cylinder radius, r2=cone base radius, rho=fraction of arrow that is cylinder.
    """

    def _fmt(p: np.ndarray) -> str:
        return f"{p[0]:.3f} {p[1]:.3f} {p[2]:.3f}"

    def _arrow_line(p_from: np.ndarray, p_to: np.ndarray, r1: float) -> str:
        r2 = max(1e-6, float(r1) * float(arrow_head_radius_scale))
        rho = float(arrow_shaft_fraction)
        return f".arrow {_fmt(p_from)} {_fmt(p_to)} {r1:.3f} {r2:.3f} {rho:.3f}"

    a1, a2 = axis.endpoints(margin=axis_margin)

    lines: List[str] = []
    lines.append(".comment Generated by angle_helical_axisV2_1.py")
    lines.append(f".comment Axis point: {_fmt(axis.point)}")
    lines.append(f".comment Axis direction (unit): {_fmt(axis.direction)}")
    lines.append(f".comment Axis drawing margin: {axis_margin:.3f} A")

    # Axis
    lines.append(".comment Draw helical axis arrow")
    lines.append(f".color {axis_color[0]:.3f} {axis_color[1]:.3f} {axis_color[2]:.3f}")
    lines.append(_arrow_line(a1, a2, axis_radius))

    # Vector 1 (radial)
    lines.append(".comment Draw radial vector 1 from the axis to Point 1")
    lines.append(f".color {v1_color[0]:.3f} {v1_color[1]:.3f} {v1_color[2]:.3f}")
    lines.append(_arrow_line(r1.proj, r1.point, vector_radius))
    lines.append(".comment Draw Point 1 marker sphere")
    lines.append(f".sphere {_fmt(r1.point)} {sphere_radius:.3f}")
    lines.append(".comment Draw Point 1 projection marker sphere on the axis")
    lines.append(f".sphere {_fmt(r1.proj)} {(sphere_radius * 0.7):.3f}")

    # Vector 2 (radial)
    lines.append(".comment Draw radial vector 2 from the axis to Point 2")
    lines.append(f".color {v2_color[0]:.3f} {v2_color[1]:.3f} {v2_color[2]:.3f}")
    lines.append(_arrow_line(r2.proj, r2.point, vector_radius))
    lines.append(".comment Draw Point 2 marker sphere")
    lines.append(f".sphere {_fmt(r2.point)} {sphere_radius:.3f}")
    lines.append(".comment Draw Point 2 projection marker sphere on the axis")
    lines.append(f".sphere {_fmt(r2.proj)} {(sphere_radius * 0.7):.3f}")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# -----------------------------
# Reporting
# -----------------------------


def format_vec(v: np.ndarray) -> str:
    return f"({v[0]: .4f}, {v[1]: .4f}, {v[2]: .4f})"


def make_report(
    pdb_path: Optional[Path],
    axis_source: str,
    axis_atoms: Sequence[str],
    axis: Axis,
    r1: RadialResult,
    r2: RadialResult,
) -> str:
    lines: List[str] = []
    lines.append(f"PDB: {pdb_path if pdb_path is not None else '(none)'}")
    lines.append(f"Axis source: {axis_source}")
    if axis_source.lower().startswith("pdb"):
        lines.append(f"Axis atoms: {', '.join(axis_atoms)}")
    lines.append("")

    lines.append("Helical axis (best-fit line / user-defined line):")
    lines.append(f"  point on axis           : {format_vec(axis.point)}")
    lines.append(f"  direction (unit vector) : {format_vec(axis.direction)}")
    lines.append("")

    def _radial_block(r: RadialResult, idx: int) -> None:
        lines.append(f"Point {idx}: {r.spec}")
        lines.append(f"  xyz                  : {format_vec(r.point)}")
        lines.append(f"  closest point on axis: {format_vec(r.proj)}")
        lines.append(f"  distance to axis (A) : {r.distance:.4f}")
        if r.unit is None:
            lines.append("  radial unit vector   : undefined (point is on axis)")
        else:
            lines.append(f"  radial unit vector   : {format_vec(r.unit)}")
        lines.append("")

    _radial_block(r1, 1)
    _radial_block(r2, 2)

    if r1.unit is None or r2.unit is None:
        lines.append("Angle between radial vectors: undefined (one or both points lie on axis).")
    else:
        ang_unsigned = angle_between_unit_vectors(r1.unit, r2.unit)
        ang_signed = signed_angle_around_axis(r1.unit, r2.unit, axis.direction)
        ang_0_360 = ang_signed % 360.0

        lines.append("Angle between radial unit vectors:")
        lines.append(f"  unsigned (0..180)                  : {ang_unsigned:.3f} deg")
        lines.append(f"  signed around axis (point1->point2): {ang_signed:.3f} deg")
        lines.append(f"  wrapped to (0..360)                : {ang_0_360:.3f} deg")
        lines.append("  (Note: signed angle depends on the chosen axis direction.)")

    return "\n".join(lines) + "\n"


# -----------------------------
# Core execution
# -----------------------------


@dataclass
class RunConfig:
    # PDB is optional in V2 (required if fitting axis from PDB, or if points are atom specs)
    pdb: Optional[Path]

    # Axis from PDB
    axis_atoms: List[str]

    # Custom axis (if provided)
    axis_point: Optional[str] = None  # xyz triplet string
    axis_vector: Optional[str] = None  # xyz triplet string

    # Points (xyz triplet or atom spec)
    point1: str = ""
    point2: str = ""

    out_bild: Path = Path("axis_vectors.bild")

    # Visual params
    axis_radius: float = 1.0
    vector_radius: float = 1.0
    sphere_radius: float = 1.25
    axis_margin: float = 5.0


def _needs_pdb_for_point(spec: str) -> bool:
    return _try_parse_xyz(spec) is None


def run(cfg: RunConfig) -> str:
    use_custom_axis = (cfg.axis_point is not None) or (cfg.axis_vector is not None)
    if use_custom_axis and (cfg.axis_point is None or cfg.axis_vector is None):
        raise ValueError("Custom axis requires BOTH axis_point and axis_vector")

    need_pdb = False
    if not use_custom_axis:
        need_pdb = True
    else:
        if _needs_pdb_for_point(cfg.point1) or _needs_pdb_for_point(cfg.point2):
            need_pdb = True

    atoms: Optional[List[Atom]] = None
    atom_index: Optional[Dict[Tuple[str, int, str], Atom]] = None

    pdb_path: Optional[Path] = None
    if need_pdb:
        if cfg.pdb is None:
            raise ValueError("A PDB file is required for this run (axis fit from PDB and/or atom points).")
        pdb_path = cfg.pdb
        atoms = load_pdb_atoms(pdb_path)
        atom_index = build_atom_index(atoms)

    # Parse points (now that we may have atom_index)
    p1 = parse_point_spec(cfg.point1, atom_index)
    p2 = parse_point_spec(cfg.point2, atom_index)

    # Build axis
    if use_custom_axis:
        axis_p = parse_xyz_triplet(cfg.axis_point or "", what="axis point")
        axis_v = parse_xyz_triplet(cfg.axis_vector or "", what="axis vector")
        axis = axis_from_point_and_vector(axis_p, axis_v, ref_points=[p1, p2])
        axis_source = "Custom (user-defined point + vector)"
    else:
        assert atoms is not None
        axis = fit_axis_pca(atoms, cfg.axis_atoms)
        axis_source = "PDB fit (PCA on axis atoms)"

    r1 = radial_vector_to_axis(p1, axis, spec=cfg.point1)
    r2 = radial_vector_to_axis(p2, axis, spec=cfg.point2)

    # Write BILD
    write_bild(
        cfg.out_bild,
        axis,
        r1,
        r2,
        axis_radius=cfg.axis_radius,
        vector_radius=cfg.vector_radius,
        sphere_radius=cfg.sphere_radius,
        axis_margin=cfg.axis_margin,
    )

    report = make_report(pdb_path, axis_source, cfg.axis_atoms, axis, r1, r2)
    report += f"\nBILD written: {cfg.out_bild}\n"
    return report


# -----------------------------
# GUI
# -----------------------------


def launch_gui(parent=None) -> None:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox
        from tkinter.scrolledtext import ScrolledText
    except Exception as e:
        print("Tkinter is not available; cannot start GUI.")
        print(str(e))
        sys.exit(1)

    owns_mainloop = parent is None
    root = tk.Tk() if owns_mainloop else tk.Toplevel(parent)
    root.title("angle_helical_axisV2_1.py")
    if parent is not None:
        root.transient(parent)

    # Variables
    var_pdb = tk.StringVar(value="")
    var_out = tk.StringVar(value="")

    # Axis source
    var_axis_source = tk.StringVar(value="Fit from PDB")  # or "Custom"
    var_axis_atoms = tk.StringVar(value="C1'")

    # Custom axis (xyz)
    var_axis_px = tk.StringVar(value="0")
    var_axis_py = tk.StringVar(value="0")
    var_axis_pz = tk.StringVar(value="0")
    var_axis_vx = tk.StringVar(value="0")
    var_axis_vy = tk.StringVar(value="0")
    var_axis_vz = tk.StringVar(value="1")

    # Point 1 mode + vars
    var_p1_mode = tk.StringVar(value="Atom")
    var_p1_chain = tk.StringVar(value="*")
    var_p1_resseq = tk.StringVar(value="")
    var_p1_atom = tk.StringVar(value="C1'")
    var_p1_x = tk.StringVar(value="")
    var_p1_y = tk.StringVar(value="")
    var_p1_z = tk.StringVar(value="")

    # Point 2 mode + vars
    var_p2_mode = tk.StringVar(value="Atom")
    var_p2_chain = tk.StringVar(value="*")
    var_p2_resseq = tk.StringVar(value="")
    var_p2_atom = tk.StringVar(value="C1'")
    var_p2_x = tk.StringVar(value="")
    var_p2_y = tk.StringVar(value="")
    var_p2_z = tk.StringVar(value="")

    # Radii defaults (per request)
    var_axis_r = tk.StringVar(value="1.0")
    var_vec_r = tk.StringVar(value="1.0")
    var_sph_r = tk.StringVar(value="1.25")
    var_axis_margin = tk.StringVar(value="5.0")

    def browse_pdb() -> None:
        fn = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=[("PDB files", "*.pdb *.ent *.pdb1 *.pdb.txt *.txt"), ("All files", "*")],
        )
        if fn:
            var_pdb.set(fn)
            # Suggest output
            base = os.path.splitext(fn)[0]
            var_out.set(base + "_axis_vectors.bild")

    def browse_out() -> None:
        fn = filedialog.asksaveasfilename(
            title="Save BILD file",
            defaultextension=".bild",
            filetypes=[("BILD files", "*.bild"), ("All files", "*")],
        )
        if fn:
            var_out.set(fn)

    def _point_spec(mode: str, chain: str, resseq: str, atom: str, x: str, y: str, z: str) -> str:
        if mode == "XYZ":
            return f"{x} {y} {z}".strip()
        # Atom
        ch = (chain or "*").strip()
        rs = resseq.strip()
        at = atom.strip()
        if not rs or not at:
            return ""
        return f"{ch}:{rs}:{at}"

    def _validate_xyz_fields(x: str, y: str, z: str, label: str) -> Optional[str]:
        spec = f"{x} {y} {z}".strip()
        if _try_parse_xyz(spec) is None:
            return f"{label} XYZ is invalid. Please enter numeric x, y, z."
        return None

    def do_run() -> None:
        axis_source = var_axis_source.get().strip()

        # Points
        p1_spec = _point_spec(
            var_p1_mode.get(),
            var_p1_chain.get(),
            var_p1_resseq.get(),
            var_p1_atom.get(),
            var_p1_x.get(),
            var_p1_y.get(),
            var_p1_z.get(),
        )
        p2_spec = _point_spec(
            var_p2_mode.get(),
            var_p2_chain.get(),
            var_p2_resseq.get(),
            var_p2_atom.get(),
            var_p2_x.get(),
            var_p2_y.get(),
            var_p2_z.get(),
        )

        if not p1_spec or not p2_spec:
            messagebox.showerror("Error", "Please provide valid Point 1 and Point 2 inputs")
            return

        # Validate xyz fields if needed
        if var_p1_mode.get() == "XYZ":
            err = _validate_xyz_fields(var_p1_x.get(), var_p1_y.get(), var_p1_z.get(), "Point 1")
            if err:
                messagebox.showerror("Error", err)
                return
        if var_p2_mode.get() == "XYZ":
            err = _validate_xyz_fields(var_p2_x.get(), var_p2_y.get(), var_p2_z.get(), "Point 2")
            if err:
                messagebox.showerror("Error", err)
                return

        out = var_out.get().strip()
        if not out:
            messagebox.showerror("Error", "Please choose an output .bild path")
            return

        # Axis config
        axis_point = None
        axis_vector = None
        axis_atoms = parse_axis_atoms(var_axis_atoms.get())

        if axis_source == "Custom axis":
            axis_point = f"{var_axis_px.get()} {var_axis_py.get()} {var_axis_pz.get()}"
            axis_vector = f"{var_axis_vx.get()} {var_axis_vy.get()} {var_axis_vz.get()}"
            # Validate custom axis xyz
            try:
                parse_xyz_triplet(axis_point, what="axis point")
                parse_xyz_triplet(axis_vector, what="axis vector")
            except Exception as e:
                messagebox.showerror("Error", str(e))
                return

        # PDB requirement
        pdb = var_pdb.get().strip()
        pdb_path: Optional[Path] = None
        need_pdb = False
        if axis_source == "Fit from PDB":
            need_pdb = True
        if var_p1_mode.get() == "Atom" or var_p2_mode.get() == "Atom":
            need_pdb = True

        if need_pdb:
            if not pdb or not os.path.isfile(pdb):
                messagebox.showerror(
                    "Error",
                    "A valid PDB file is required (axis fitted from PDB and/or points specified as atoms).",
                )
                return
            pdb_path = Path(pdb)

        try:
            cfg = RunConfig(
                pdb=pdb_path,
                axis_atoms=axis_atoms,
                axis_point=axis_point,
                axis_vector=axis_vector,
                point1=p1_spec,
                point2=p2_spec,
                out_bild=Path(out),
                axis_radius=float(var_axis_r.get()),
                vector_radius=float(var_vec_r.get()),
                sphere_radius=float(var_sph_r.get()),
                axis_margin=float(var_axis_margin.get()),
            )
            report = run(cfg)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        text.delete("1.0", "end")
        text.insert("1.0", report)

    frm = tk.Frame(root, padx=10, pady=10)
    frm.pack(fill="both", expand=True)

    # Row helpers
    def add_entry_row(row: int, label: str, var: tk.StringVar, browse_cmd=None, width: int = 55):
        tk.Label(frm, text=label, anchor="w").grid(row=row, column=0, sticky="w")
        ent = tk.Entry(frm, textvariable=var, width=width)
        ent.grid(row=row, column=1, sticky="we", padx=(5, 5))
        if browse_cmd is not None:
            tk.Button(frm, text="Browse", command=browse_cmd).grid(row=row, column=2, padx=(0, 5))
        return ent

    # PDB file row
    add_entry_row(0, "PDB file (optional if using custom axis + XYZ points)", var_pdb, browse_pdb)

    # Axis source row
    tk.Label(frm, text="Axis source", anchor="w").grid(row=1, column=0, sticky="w")
    opt_axis = tk.OptionMenu(frm, var_axis_source, "Fit from PDB", "Custom axis")
    opt_axis.grid(row=1, column=1, sticky="w", padx=(5, 5))

    # Axis atoms row (PDB fit)
    tk.Label(frm, text="Axis atoms (PDB fit only)", anchor="w").grid(row=2, column=0, sticky="w")
    ent_axis_atoms = tk.Entry(frm, textvariable=var_axis_atoms, width=20)
    ent_axis_atoms.grid(row=2, column=1, sticky="w", padx=(5, 5))

    # Custom axis frame
    frm_custom = tk.Frame(frm)
    frm_custom.grid(row=3, column=0, columnspan=3, sticky="we", pady=(2, 2))

    tk.Label(frm_custom, text="Custom axis point (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_custom, textvariable=var_axis_px, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_custom, textvariable=var_axis_py, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_custom, textvariable=var_axis_pz, width=8).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_custom, text="Custom axis vector (vx vy vz)").grid(row=1, column=0, sticky="w")
    tk.Entry(frm_custom, textvariable=var_axis_vx, width=8).grid(row=1, column=1, padx=(5, 2))
    tk.Entry(frm_custom, textvariable=var_axis_vy, width=8).grid(row=1, column=2, padx=(2, 2))
    tk.Entry(frm_custom, textvariable=var_axis_vz, width=8).grid(row=1, column=3, padx=(2, 2))

    # Point 1 controls
    tk.Label(frm, text="Point 1 type", anchor="w").grid(row=4, column=0, sticky="w")
    tk.OptionMenu(frm, var_p1_mode, "Atom", "XYZ").grid(row=4, column=1, sticky="w", padx=(5, 5))

    frm_p1_atom = tk.Frame(frm)
    frm_p1_xyz = tk.Frame(frm)
    frm_p1_atom.grid(row=5, column=0, columnspan=3, sticky="we")
    frm_p1_xyz.grid(row=5, column=0, columnspan=3, sticky="we")

    tk.Label(frm_p1_atom, text="Point 1 atom (chain, resSeq, atom)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p1_atom, textvariable=var_p1_chain, width=5).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p1_atom, textvariable=var_p1_resseq, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p1_atom, textvariable=var_p1_atom, width=10).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_p1_xyz, text="Point 1 xyz (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p1_xyz, textvariable=var_p1_x, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p1_xyz, textvariable=var_p1_y, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p1_xyz, textvariable=var_p1_z, width=8).grid(row=0, column=3, padx=(2, 2))

    # Point 2 controls
    tk.Label(frm, text="Point 2 type", anchor="w").grid(row=6, column=0, sticky="w")
    tk.OptionMenu(frm, var_p2_mode, "Atom", "XYZ").grid(row=6, column=1, sticky="w", padx=(5, 5))

    frm_p2_atom = tk.Frame(frm)
    frm_p2_xyz = tk.Frame(frm)
    frm_p2_atom.grid(row=7, column=0, columnspan=3, sticky="we")
    frm_p2_xyz.grid(row=7, column=0, columnspan=3, sticky="we")

    tk.Label(frm_p2_atom, text="Point 2 atom (chain, resSeq, atom)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p2_atom, textvariable=var_p2_chain, width=5).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p2_atom, textvariable=var_p2_resseq, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p2_atom, textvariable=var_p2_atom, width=10).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_p2_xyz, text="Point 2 xyz (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p2_xyz, textvariable=var_p2_x, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p2_xyz, textvariable=var_p2_y, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p2_xyz, textvariable=var_p2_z, width=8).grid(row=0, column=3, padx=(2, 2))

    # Output
    add_entry_row(8, "Output .bild", var_out, browse_out)

    # Axis drawing margin and radii
    tk.Label(frm, text="Axis drawing margin", anchor="w").grid(row=9, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_axis_margin, width=10).grid(row=9, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Axis radius", anchor="w").grid(row=10, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_axis_r, width=10).grid(row=10, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Vector radius", anchor="w").grid(row=11, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_vec_r, width=10).grid(row=11, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Sphere radius", anchor="w").grid(row=12, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_sph_r, width=10).grid(row=12, column=1, sticky="w", padx=(5, 5))

    tk.Button(frm, text="Run", command=do_run).grid(row=13, column=0, pady=(8, 8), sticky="w")

    # Output text
    text = ScrolledText(frm, height=18, width=90)
    text.grid(row=14, column=0, columnspan=3, sticky="nsew", pady=(5, 0))

    frm.columnconfigure(1, weight=1)
    frm.rowconfigure(14, weight=1)

    def refresh_visibility(*_args) -> None:
        # Axis controls
        if var_axis_source.get() == "Custom axis":
            frm_custom.grid()
            ent_axis_atoms.configure(state="disabled")
        else:
            frm_custom.grid_remove()
            ent_axis_atoms.configure(state="normal")

        # Point 1
        if var_p1_mode.get() == "XYZ":
            frm_p1_atom.grid_remove()
            frm_p1_xyz.grid()
        else:
            frm_p1_xyz.grid_remove()
            frm_p1_atom.grid()

        # Point 2
        if var_p2_mode.get() == "XYZ":
            frm_p2_atom.grid_remove()
            frm_p2_xyz.grid()
        else:
            frm_p2_xyz.grid_remove()
            frm_p2_atom.grid()

    var_axis_source.trace_add("write", refresh_visibility)
    var_p1_mode.trace_add("write", refresh_visibility)
    var_p2_mode.trace_add("write", refresh_visibility)

    refresh_visibility()

    text.insert(
        "1.0",
        "Point inputs:\n"
        "  - Atom: chain/resSeq/atom (atom names accept prime or star, e.g. C1' or C1*)\n"
        "  - XYZ : numeric x y z\n"
        "Axis:\n"
        "  - Fit from PDB: uses --axis-atoms (default C1')\n"
        "  - Custom axis: provide axis point + axis vector (vector is normalized automatically)\n"
        "  - Axis drawing margin controls how far the displayed axis extends beyond the selected points/fit range\n",
    )

    if owns_mainloop:
        root.mainloop()


# -----------------------------
# CLI
# -----------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="angle_helical_axisV2_1.py",
        description="Fit (or define) a helical axis and compute radial vectors/angle for two points.",
    )

    p.add_argument("--gui", action="store_true", help="Launch GUI mode.")

    p.add_argument("-i", "--pdb", type=str, default=None, help="Input PDB file (required in some modes)")
    p.add_argument(
        "--axis-atoms",
        type=str,
        default="C1'",
        help="Atom name(s) used to fit the axis (comma/space separated). Default: C1'",
    )

    p.add_argument(
        "--axis-point",
        type=str,
        default=None,
        help="Custom axis: point on axis as xyz triplet, e.g. '0 0 0'. Requires --axis-vector.",
    )
    p.add_argument(
        "--axis-vector",
        type=str,
        default=None,
        help="Custom axis: direction vector as xyz triplet, e.g. '0 0 1'. Requires --axis-point.",
    )

    p.add_argument("--point1", type=str, default=None, help="Point 1: xyz or CHAIN:RESSEQ:ATOM")
    p.add_argument("--point2", type=str, default=None, help="Point 2: xyz or CHAIN:RESSEQ:ATOM")

    p.add_argument(
        "-o",
        "--out-bild",
        type=str,
        default=None,
        help="Output .bild file (default: <pdb_stem>_axis_vectors.bild, or axis_vectors.bild if no PDB)",
    )

    # Defaults per request
    p.add_argument("--axis-margin", type=float, default=5.0, help="Extra margin in Angstroms added to each end of the drawn axis")
    p.add_argument("--axis-radius", type=float, default=1.0, help="Axis arrow radius")
    p.add_argument("--vector-radius", type=float, default=1.0, help="Vector arrow radius")
    p.add_argument("--sphere-radius", type=float, default=1.25, help="Sphere radius")

    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    # GUI mode if no args at all OR explicit --gui
    if len(argv) == 0 or "--gui" in argv:
        launch_gui()
        return

    parser = build_arg_parser()
    args = parser.parse_args(list(argv))

    if args.point1 is None or args.point2 is None:
        parser.error("CLI mode requires --point1 and --point2 (or use --gui).")

    use_custom_axis = (args.axis_point is not None) or (args.axis_vector is not None)
    if use_custom_axis and (args.axis_point is None or args.axis_vector is None):
        parser.error("Custom axis requires BOTH --axis-point and --axis-vector")

    pdb_path: Optional[Path] = Path(args.pdb) if args.pdb else None

    # Determine whether a PDB is required
    need_pdb = False
    if not use_custom_axis:
        need_pdb = True
    else:
        if _needs_pdb_for_point(args.point1) or _needs_pdb_for_point(args.point2):
            need_pdb = True

    if need_pdb:
        if pdb_path is None:
            parser.error("This run requires --pdb (axis fit from PDB and/or points are atom specs).")
        if not pdb_path.is_file():
            parser.error(f"PDB file not found: {pdb_path}")

    # Output naming
    if args.out_bild:
        out_bild = Path(args.out_bild)
    else:
        if pdb_path is not None:
            out_bild = pdb_path.with_suffix("").with_name(pdb_path.stem + "_axis_vectors.bild")
        else:
            out_bild = Path("axis_vectors.bild")

    cfg = RunConfig(
        pdb=pdb_path,
        axis_atoms=parse_axis_atoms(args.axis_atoms),
        axis_point=args.axis_point,
        axis_vector=args.axis_vector,
        point1=args.point1,
        point2=args.point2,
        out_bild=out_bild,
        axis_radius=float(args.axis_radius),
        vector_radius=float(args.vector_radius),
        sphere_radius=float(args.sphere_radius),
        axis_margin=float(args.axis_margin),
    )

    report = run(cfg)
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
