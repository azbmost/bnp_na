#!/usr/bin/env python3
# angle_helical_axisV2_2.py
#
# Fit (or accept a user-defined) helical axis and compute radial vectors from that axis
# to two points. Points can be given as PDB atoms (CHAIN:RESSEQ:ATOM) or explicit xyz.
# Optionally determine a region-defined 2-fold symmetry axis perpendicular to the
# helical axis. Writes a Chimera/ChimeraX .bild file that draws the axis, vectors,
# and optional symmetry-axis helper points. The axis drawing margin can be adjusted
# from the GUI or CLI.
#
# Example (fit axis from PDB):
#   python angle_helical_axisV2_2.py -i helix.pdb --point1 A:5:C1' --point2 B:18:C1' \
#       -o helix_axis_vectors.bild
#   python angle_helical_axisV2_2.py -i helix.pdb --axis_range "A1-A35,B60-B26" \
#       --point1 A:5:C1' --point2 B:56:C1' -o selected_axis_vectors.bild
#
# Example (custom axis, no PDB needed):
#   python angle_helical_axisV2_2.py --axis-point "0 0 0" --axis-vector "0 0 1" \
#       --point1 "1 0 0" --point2 "0 1 0" -o custom_axis_vectors.bild
#
# Example (region-defined 2-fold symmetry axis):
#   python angle_helical_axisV2_2.py -i helix.pdb --symmetry-regions "A1-A35, B26-B60" \
#       -o helix_twofold_axis.bild

"""angle_helical_axisV2_2.py

Identify an approximate nucleic-acid helical axis from a PDB (PCA/SVD fit), OR use a
user-supplied axis (point + direction vector). Then compute the perpendicular (radial)
vectors from that axis to two user-specified points. Optionally, determine a
region-defined 2-fold symmetry axis perpendicular to the helical axis.

Reports:
  * the helical axis (point on axis + unit direction)
  * distance of each point to the axis
  * unit radial vectors for each point
  * angle between the two radial vectors (degrees)
  * optional 2-fold symmetry axis, two points on that axis, and two points rotated
    90 degrees around the helical axis

Points can be specified either as:
  * explicit coordinates:  "x,y,z"  or "x y z"
  * an atom in a PDB:      "CHAIN:RESSEQ:ATOM" (also accepts / , @ and spaces)

Axis options:
  * Default: fit axis from PDB using --axis-atoms (default C1')
  * Limit a PDB fit with --axis-range/--axis_range (for example A1-A35,B60-B26).
    The written start-to-end order of the first range sets the positive direction.
  * Custom axis: provide BOTH --axis-point and --axis-vector (PDB not required unless
    you specify points by atom).

The script also writes a Chimera/ChimeraX .bild file that draws:
  * the axis as an arrow
  * the two radial vectors as arrows from the axis to each point
  * explanatory .comment lines before each drawn object

The displayed axis length can be extended or shortened with an axis drawing margin.

2-fold symmetry axis:
  * provide --symmetry-regions as two residue ranges, e.g. "A1-A35, B26-B60"
  * the first range is paired with the second range in reverse residue order
  * residue centers are calculated from --symmetry-atoms (default C1')

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
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# Keep multiple filename patterns as separate Tcl list elements on platforms
# where filtering is safe.  Older Aqua Tk 8.6 builds can abort in native code
# while converting an otherwise harmless extension to a macOS UTType, so the
# dialog helper below deliberately omits filetypes on macOS.
PDB_OPEN_FILETYPES = (
    ("PDB files", ("*.pdb", "*.ent", "*.pdb1", "*.pdb.txt", "*.txt")),
    ("All files", "*"),
)
BILD_SAVE_FILETYPES = (("BILD files", "*.bild"), ("All files", "*"))


def _file_dialog_options(
    title: str,
    filetypes: object,
    *,
    defaultextension: Optional[str] = None,
    platform: Optional[str] = None,
) -> Dict[str, object]:
    """Build options without unsafe native file-type filtering on macOS."""
    options: Dict[str, object] = {"title": title}
    if defaultextension is not None:
        options["defaultextension"] = defaultextension
    if (sys.platform if platform is None else platform) != "darwin":
        options["filetypes"] = filetypes
    return options


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


def build_residue_atom_index(atoms: Sequence[Atom]) -> Dict[Tuple[str, int], List[Atom]]:
    """Index atoms by (chainID, resSeq)."""
    idx: Dict[Tuple[str, int], List[Atom]] = {}
    for a in atoms:
        idx.setdefault((a.chainID, a.resSeq), []).append(a)
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
    residue_ranges: Optional[Sequence[Sequence["ResidueRef"]]] = None,
    prefer_positive_correlation: bool = True,
) -> Axis:
    """Fit a best-fit line (axis) by PCA/SVD of selected atom coordinates."""
    alias_set = set()
    for nm in axis_atom_names:
        alias_set.update(atom_aliases(nm))

    selected_residues = None
    if residue_ranges is not None:
        selected_residues = {
            (ref.chainID, ref.resSeq)
            for residue_range in residue_ranges
            for ref in residue_range
        }

    pts: List[np.ndarray] = []
    chain_resseq: List[Tuple[str, int]] = []

    for a in atoms:
        if a.name in alias_set and (
            selected_residues is None or (a.chainID, a.resSeq) in selected_residues
        ):
            pts.append(a.xyz)
            chain_resseq.append((a.chainID, a.resSeq))

    if len(pts) < 3:
        range_note = " in the selected axis range(s)" if selected_residues is not None else ""
        raise ValueError(
            f"Need at least 3 atoms to fit axis; found {len(pts)} atoms matching "
            f"{sorted(alias_set)}{range_note}"
        )

    X = np.stack(pts, axis=0)
    mean = X.mean(axis=0)
    Xc = X - mean

    # Principal component via SVD
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    direction = _normalize(vt[0])

    # When ranges are explicit, their first written range is authoritative:
    # orient the axis from that range's start toward its end. Average multiple
    # selected axis atoms within a residue before choosing the sign.
    if residue_ranges is not None:
        first_range = residue_ranges[0]
        points_by_residue: Dict[Tuple[str, int], List[np.ndarray]] = {}
        for point, residue_key in zip(pts, chain_resseq):
            points_by_residue.setdefault(residue_key, []).append(point)
        ordered_centers = [
            np.stack(points_by_residue[(ref.chainID, ref.resSeq)], axis=0).mean(axis=0)
            for ref in first_range
            if (ref.chainID, ref.resSeq) in points_by_residue
        ]
        if len(ordered_centers) < 2:
            raise ValueError(
                "The first axis range must contain selected axis atoms in at least two residues "
                "so its start-to-end direction can be determined."
            )
        direction_anchor = ordered_centers[-1] - ordered_centers[0]
        if float(np.linalg.norm(direction_anchor)) < 1e-12:
            raise ValueError(
                "The first axis range has coincident start/end axis-atom centers; "
                "its positive direction cannot be determined."
            )
        if float(np.dot(direction, direction_anchor)) < 0.0:
            direction = -direction

    # Without an explicit range, retain the original residue-number heuristic.
    elif prefer_positive_correlation and len(chain_resseq) >= 6:
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


def rotation_matrix(axis_dir: np.ndarray, angle_rad: float) -> np.ndarray:
    """Return a 3x3 Rodrigues rotation matrix."""
    u = _normalize(np.array(axis_dir, dtype=float))
    ux, uy, uz = float(u[0]), float(u[1]), float(u[2])
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array(
        [
            [c + ux * ux * one_c, ux * uy * one_c - uz * s, ux * uz * one_c + uy * s],
            [uy * ux * one_c + uz * s, c + uy * uy * one_c, uy * uz * one_c - ux * s],
            [uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c],
        ],
        dtype=float,
    )


def rotate_vector_around_axis(v: np.ndarray, axis_dir: np.ndarray, angle_degrees: float) -> np.ndarray:
    """Rotate vector v around axis_dir by angle_degrees."""
    return rotation_matrix(axis_dir, math.radians(angle_degrees)) @ np.array(v, dtype=float)


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
# 2-fold symmetry axis from paired residue regions
# -----------------------------


@dataclass(frozen=True)
class ResidueRef:
    chainID: str
    resSeq: int

    def label(self) -> str:
        chain = self.chainID if self.chainID != " " else "_"
        return f"{chain}{self.resSeq}"


@dataclass
class SymmetryPair:
    first: ResidueRef
    second: ResidueRef
    first_center: np.ndarray
    second_center: np.ndarray
    midpoint: np.ndarray


@dataclass
class SymmetryAxisResult:
    regions: str
    region_source: str
    atom_names: List[str]
    pairs: List[SymmetryPair]
    center: np.ndarray
    direction: np.ndarray
    perpendicular_direction: np.ndarray
    radius: float
    point_plus: np.ndarray
    point_minus: np.ndarray
    rotated_point_plus: np.ndarray
    rotated_point_minus: np.ndarray
    midpoint_line_rmsd: float
    rotation_rmsd: float


_REGION_RE = re.compile(
    r"^\s*([A-Za-z0-9])\s*:?\s*(-?\d+)\s*-\s*(?:([A-Za-z0-9])\s*:?\s*)?(-?\d+)\s*$"
)


def _parse_residue_range(text: str) -> List[ResidueRef]:
    m = _REGION_RE.match(text)
    if not m:
        raise ValueError(
            f"Could not parse residue range {text!r}. Expected forms like A1-A35 or A:1-A:35."
        )

    chain1 = m.group(1)
    start = int(m.group(2))
    chain2 = m.group(3) or chain1
    end = int(m.group(4))
    if chain2 != chain1:
        raise ValueError(f"Residue range {text!r} crosses chains; use one chain per range.")

    step = 1 if end >= start else -1
    return [ResidueRef(chain1, rs) for rs in range(start, end + step, step)]


def parse_symmetry_regions(text: str) -> Tuple[List[ResidueRef], List[ResidueRef]]:
    """Parse two residue ranges used to define the 2-fold axis.

    The second range is returned in reverse order so A1-A35, B26-B60 pairs
    A1 with B60, A2 with B59, ..., A35 with B26.
    """
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if len(parts) != 2:
        raise ValueError(
            "Provide exactly two residue ranges for --symmetry-regions, "
            'for example: "A1-A35, B26-B60".'
        )

    first = _parse_residue_range(parts[0])
    second = list(reversed(_parse_residue_range(parts[1])))
    if len(first) != len(second):
        raise ValueError(
            f"Symmetry ranges must contain the same number of residues; "
            f"found {len(first)} and {len(second)}."
        )
    if len(first) < 2:
        raise ValueError("Symmetry regions need at least two paired residues.")
    return first, second


def parse_axis_ranges(text: str) -> List[List[ResidueRef]]:
    """Parse one or more residue ranges used to select and orient a PDB axis fit."""
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        raise ValueError(
            "Provide at least one residue range for --axis-range/--axis_range, "
            'for example: "A1-A35,B60-B26".'
        )
    ranges = [_parse_residue_range(part) for part in parts]
    if len(ranges[0]) < 2:
        raise ValueError(
            "The first axis range must contain at least two residues because its "
            "start-to-end order sets the positive axis direction."
        )
    return ranges


def _format_residue_range(refs: Sequence[ResidueRef]) -> str:
    if not refs:
        return "(none)"
    if len(refs) == 1:
        return refs[0].label()
    return f"{refs[0].label()}-{refs[-1].label()}"


def _atom_name_alias_set(atom_names: Sequence[str]) -> Tuple[set, bool]:
    use_all = False
    alias_set = set()
    for nm in atom_names:
        clean = nm.strip()
        if not clean:
            continue
        if clean.lower() in {"*", "all"}:
            use_all = True
            continue
        alias_set.update(atom_aliases(clean))
    return alias_set, use_all


def residue_center(
    ref: ResidueRef,
    residue_index: Dict[Tuple[str, int], List[Atom]],
    atom_names: Sequence[str],
) -> np.ndarray:
    atoms = residue_index.get((ref.chainID, ref.resSeq))
    if not atoms:
        raise ValueError(f"No atoms found for residue {ref.label()}.")

    alias_set, use_all = _atom_name_alias_set(atom_names)
    if use_all:
        selected = atoms
    else:
        selected = [a for a in atoms if a.name in alias_set]

    if not selected:
        atom_text = ", ".join(atom_names)
        raise ValueError(f"No symmetry atom(s) {atom_text!r} found for residue {ref.label()}.")

    return np.stack([a.xyz for a in selected], axis=0).mean(axis=0)


def whole_model_symmetry_regions(
    atoms: Sequence[Atom],
    atom_names: Sequence[str],
) -> Tuple[List[ResidueRef], List[ResidueRef], str]:
    """Infer two paired residue ranges from the whole model.

    Blank symmetry-region input uses this path. It expects exactly two chains
    containing the selected symmetry atoms, with the same number of residues.
    """
    alias_set, use_all = _atom_name_alias_set(atom_names)
    chain_order: List[str] = []
    chain_refs: Dict[str, List[ResidueRef]] = {}
    seen_residues = set()

    for atom in atoms:
        if not use_all and atom.name not in alias_set:
            continue
        key = (atom.chainID, atom.resSeq)
        if key in seen_residues:
            continue
        seen_residues.add(key)
        if atom.chainID not in chain_refs:
            chain_order.append(atom.chainID)
            chain_refs[atom.chainID] = []
        chain_refs[atom.chainID].append(ResidueRef(atom.chainID, atom.resSeq))

    chains = [ch for ch in chain_order if chain_refs.get(ch)]
    if len(chains) != 2:
        chain_text = ", ".join(ch if ch != " " else "_" for ch in chains) or "(none)"
        atom_text = ", ".join(atom_names)
        raise ValueError(
            "Blank 2-fold symmetry regions use the whole model and require exactly two chains "
            f"with selected symmetry atoms ({atom_text}); found {len(chains)} chains: {chain_text}. "
            "Specify explicit regions instead."
        )

    first = sorted(chain_refs[chains[0]], key=lambda ref: ref.resSeq)
    second_forward = sorted(chain_refs[chains[1]], key=lambda ref: ref.resSeq)
    if len(first) != len(second_forward):
        raise ValueError(
            "Blank 2-fold symmetry regions use the whole model and require the two selected chains "
            f"to have the same residue count; found {len(first)} and {len(second_forward)}. "
            "Specify explicit regions instead."
        )
    if len(first) < 2:
        raise ValueError("Whole-model 2-fold symmetry needs at least two paired residues per chain.")

    effective = f"{_format_residue_range(first)}, {_format_residue_range(second_forward)}"
    return first, list(reversed(second_forward)), effective


def build_symmetry_pairs(
    atoms: Sequence[Atom],
    regions: Optional[str],
    atom_names: Sequence[str],
) -> Tuple[List[SymmetryPair], str, str]:
    if _has_text(regions):
        first_refs, second_refs = parse_symmetry_regions(regions or "")
        effective_regions = regions or ""
        region_source = "explicit"
    else:
        first_refs, second_refs, effective_regions = whole_model_symmetry_regions(atoms, atom_names)
        region_source = "whole model"

    residue_index = build_residue_atom_index(atoms)
    pairs: List[SymmetryPair] = []
    for ref1, ref2 in zip(first_refs, second_refs):
        p1 = residue_center(ref1, residue_index, atom_names)
        p2 = residue_center(ref2, residue_index, atom_names)
        pairs.append(
            SymmetryPair(
                first=ref1,
                second=ref2,
                first_center=p1,
                second_center=p2,
                midpoint=(p1 + p2) * 0.5,
            )
        )
    return pairs, effective_regions, region_source


def determine_twofold_symmetry_axis(
    atoms: Sequence[Atom],
    axis: Axis,
    regions: Optional[str],
    atom_names: Sequence[str],
    *,
    point_radius: Optional[float] = None,
) -> SymmetryAxisResult:
    """Determine a 2-fold symmetry axis from two paired residue ranges.

    The resulting line is constrained to be perpendicular to `axis`. The first
    residue range is paired with the second range in reverse order. Pair midpoints
    define the perpendicular 2-fold direction; a 180-degree rotation around the
    fitted line is then scored by RMSD between paired residue centers.
    """
    pairs, effective_regions, region_source = build_symmetry_pairs(atoms, regions, atom_names)
    helix_u = axis.direction

    midpoints = np.stack([p.midpoint for p in pairs], axis=0)
    t_vals = np.array([float(np.dot(m - axis.point, helix_u)) for m in midpoints], dtype=float)
    t_center = float(np.mean(t_vals))
    center = axis.point + t_center * helix_u

    offsets = midpoints - center
    offsets_perp = offsets - np.outer(offsets @ helix_u, helix_u)
    if float(np.max(np.linalg.norm(offsets_perp, axis=1))) < 1e-8:
        raise ValueError(
            "Could not determine a 2-fold symmetry-axis direction: paired-region midpoints "
            "fall too close to the helical axis."
        )

    _, _, vt = np.linalg.svd(offsets_perp, full_matrices=False)
    direction = vt[0]
    direction = direction - float(np.dot(direction, helix_u)) * helix_u
    direction = _normalize(direction)

    mean_offset = offsets_perp.mean(axis=0)
    if float(np.linalg.norm(mean_offset)) > 1e-8:
        if float(np.dot(direction, mean_offset)) < 0.0:
            direction = -direction
    else:
        for offset in offsets_perp:
            if float(np.linalg.norm(offset)) > 1e-8:
                if float(np.dot(direction, offset)) < 0.0:
                    direction = -direction
                break

    perp_dir = rotate_vector_around_axis(direction, helix_u, 90.0)
    perp_dir = perp_dir - float(np.dot(perp_dir, helix_u)) * helix_u
    perp_dir = _normalize(perp_dir)

    coords = offsets_perp @ direction
    if point_radius is not None:
        radius = float(point_radius)
        if radius <= 0:
            raise ValueError("--symmetry-point-radius must be positive when provided.")
    else:
        radius = float(np.max(np.abs(coords)))
        if radius < 1e-6:
            selected_centers = [p.first_center for p in pairs] + [p.second_center for p in pairs]
            distances = [radial_vector_to_axis(c, axis).distance for c in selected_centers]
            radius = float(np.mean(distances)) if distances else 10.0
        if radius < 1e-6:
            radius = 10.0

    point_plus = center + radius * direction
    point_minus = center - radius * direction
    rotated_point_plus = center + radius * perp_dir
    rotated_point_minus = center - radius * perp_dir

    line_residuals = []
    for midpoint in midpoints:
        delta = midpoint - center
        on_line = float(np.dot(delta, direction)) * direction
        line_residuals.append(delta - on_line)
    midpoint_line_rmsd = float(np.sqrt(np.mean([float(np.dot(v, v)) for v in line_residuals])))

    rot180 = rotation_matrix(direction, math.pi)
    rotation_diffs = []
    for pair in pairs:
        rotated = center + rot180 @ (pair.first_center - center)
        diff = rotated - pair.second_center
        rotation_diffs.append(float(np.dot(diff, diff)))
    rotation_rmsd = float(np.sqrt(np.mean(rotation_diffs)))

    return SymmetryAxisResult(
        regions=effective_regions,
        region_source=region_source,
        atom_names=list(atom_names),
        pairs=pairs,
        center=center,
        direction=direction,
        perpendicular_direction=perp_dir,
        radius=radius,
        point_plus=point_plus,
        point_minus=point_minus,
        rotated_point_plus=rotated_point_plus,
        rotated_point_minus=rotated_point_minus,
        midpoint_line_rmsd=midpoint_line_rmsd,
        rotation_rmsd=rotation_rmsd,
    )


# -----------------------------
# BILD output (Chimera/ChimeraX)
# -----------------------------


DEFAULT_AXIS_COLOR = (0.2, 0.2, 1.0)
DEFAULT_POINT1_COLOR = (1.0, 0.2, 0.2)
DEFAULT_POINT2_COLOR = (0.2, 0.8, 0.2)
DEFAULT_SYMMETRY_COLOR = (0.9, 0.1, 0.8)
DEFAULT_SYMMETRY_PLANE_COLOR = (0.0, 0.75, 0.9)
DEFAULT_SYMMETRY_POINT_RADIUS = 15.0


def write_bild(
    out_path: Path,
    axis: Axis,
    r1: Optional[RadialResult] = None,
    r2: Optional[RadialResult] = None,
    *,
    symmetry: Optional[SymmetryAxisResult] = None,
    axis_color: Tuple[float, float, float] = DEFAULT_AXIS_COLOR,
    v1_color: Tuple[float, float, float] = DEFAULT_POINT1_COLOR,
    v2_color: Tuple[float, float, float] = DEFAULT_POINT2_COLOR,
    symmetry_color: Tuple[float, float, float] = DEFAULT_SYMMETRY_COLOR,
    symmetry_plane_color: Tuple[float, float, float] = DEFAULT_SYMMETRY_PLANE_COLOR,
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
    lines.append(".comment Generated by angle_helical_axisV2_2.py")
    lines.append(f".comment Axis point: {_fmt(axis.point)}")
    lines.append(f".comment Axis direction (unit): {_fmt(axis.direction)}")
    lines.append(f".comment Axis drawing margin: {axis_margin:.3f} A")

    # Axis
    lines.append(".comment Draw helical axis arrow")
    lines.append(f".color {axis_color[0]:.3f} {axis_color[1]:.3f} {axis_color[2]:.3f}")
    lines.append(_arrow_line(a1, a2, axis_radius))

    if r1 is not None:
        # Vector 1 (radial)
        lines.append(".comment Draw radial vector 1 from the axis to Point 1")
        lines.append(f".color {v1_color[0]:.3f} {v1_color[1]:.3f} {v1_color[2]:.3f}")
        lines.append(_arrow_line(r1.proj, r1.point, vector_radius))
        lines.append(".comment Draw Point 1 marker sphere")
        lines.append(f".sphere {_fmt(r1.point)} {sphere_radius:.3f}")
        lines.append(".comment Draw Point 1 projection marker sphere on the axis")
        lines.append(f".sphere {_fmt(r1.proj)} {(sphere_radius * 0.7):.3f}")

    if r2 is not None:
        # Vector 2 (radial)
        lines.append(".comment Draw radial vector 2 from the axis to Point 2")
        lines.append(f".color {v2_color[0]:.3f} {v2_color[1]:.3f} {v2_color[2]:.3f}")
        lines.append(_arrow_line(r2.proj, r2.point, vector_radius))
        lines.append(".comment Draw Point 2 marker sphere")
        lines.append(f".sphere {_fmt(r2.point)} {sphere_radius:.3f}")
        lines.append(".comment Draw Point 2 projection marker sphere on the axis")
        lines.append(f".sphere {_fmt(r2.proj)} {(sphere_radius * 0.7):.3f}")

    if symmetry is not None:
        sym_r = max(1e-6, vector_radius)
        sym_sphere = max(1e-6, sphere_radius)
        lines.append(".comment Draw 2-fold symmetry axis arrow from - point to + point")
        lines.append(
            f".comment Symmetry center: {_fmt(symmetry.center)}; "
            f"direction: {_fmt(symmetry.direction)}"
        )
        lines.append(
            f".color {symmetry_color[0]:.3f} {symmetry_color[1]:.3f} {symmetry_color[2]:.3f}"
        )
        lines.append(_arrow_line(symmetry.point_minus, symmetry.point_plus, sym_r))
        lines.append(".comment Draw the two requested points on the 2-fold symmetry axis")
        lines.append(f".sphere {_fmt(symmetry.point_plus)} {sym_sphere:.3f}")
        lines.append(f".sphere {_fmt(symmetry.point_minus)} {sym_sphere:.3f}")
        lines.append(".comment Draw symmetry-axis center on the helical axis")
        lines.append(f".sphere {_fmt(symmetry.center)} {(sym_sphere * 0.8):.3f}")

        lines.append(".comment Draw +90-degree rotated arrow from - rotated point to + rotated point")
        lines.append(
            f".color {symmetry_plane_color[0]:.3f} {symmetry_plane_color[1]:.3f} {symmetry_plane_color[2]:.3f}"
        )
        lines.append(_arrow_line(symmetry.rotated_point_minus, symmetry.rotated_point_plus, sym_r))
        lines.append(f".sphere {_fmt(symmetry.rotated_point_plus)} {sym_sphere:.3f}")
        lines.append(f".sphere {_fmt(symmetry.rotated_point_minus)} {sym_sphere:.3f}")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# -----------------------------
# Reporting
# -----------------------------


def format_vec(v: np.ndarray) -> str:
    return f"({v[0]: .4f}, {v[1]: .4f}, {v[2]: .4f})"


def format_color(name: str, rgb: Tuple[float, float, float]) -> str:
    return f"{name} RGB({rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f})"


def make_report(
    pdb_path: Optional[Path],
    axis_source: str,
    axis_atoms: Sequence[str],
    axis_range: Optional[str],
    axis: Axis,
    r1: Optional[RadialResult],
    r2: Optional[RadialResult],
    symmetry: Optional[SymmetryAxisResult] = None,
) -> str:
    lines: List[str] = []
    lines.append(f"PDB: {pdb_path if pdb_path is not None else '(none)'}")
    lines.append(f"Axis source: {axis_source}")
    if axis_source.lower().startswith("pdb"):
        lines.append(f"Axis atoms: {', '.join(axis_atoms)}")
        lines.append(f"Axis residue ranges: {axis_range if _has_text(axis_range) else '(all residues)'}")
        if _has_text(axis_range):
            lines.append("Axis positive direction: start-to-end order of the first axis range")
    lines.append("")

    lines.append("Helical axis (best-fit line / user-defined line):")
    lines.append(f"  point on axis           : {format_vec(axis.point)}")
    lines.append(f"  direction (unit vector) : {format_vec(axis.direction)}")
    lines.append("")

    lines.append("BILD arrow colors/directions:")
    lines.append(
        f"  helical axis arrow      : {format_color('blue', DEFAULT_AXIS_COLOR)}; "
        "travels along the chosen helical-axis direction"
    )
    if r1 is not None:
        lines.append(
            f"  Point 1 radial arrow    : {format_color('red', DEFAULT_POINT1_COLOR)}; "
            "travels from the helical-axis projection to Point 1"
        )
    if r2 is not None:
        lines.append(
            f"  Point 2 radial arrow    : {format_color('green', DEFAULT_POINT2_COLOR)}; "
            "travels from the helical-axis projection to Point 2"
        )
    if symmetry is not None:
        lines.append(
            f"  2-fold symmetry arrow   : {format_color('magenta', DEFAULT_SYMMETRY_COLOR)}; "
            "travels from - point to + point"
        )
        lines.append(
            f"  +90-deg rotated arrow   : {format_color('cyan', DEFAULT_SYMMETRY_PLANE_COLOR)}; "
            "travels from - rotated point to + rotated point"
        )
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

    if r1 is not None and r2 is not None:
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
            lines.append("")
    else:
        lines.append("Point radial-vector angle: not calculated (Point 1/Point 2 not both provided).")
        lines.append("")

    if symmetry is not None:
        plus_radial = radial_vector_to_axis(symmetry.point_plus, axis)
        minus_radial = radial_vector_to_axis(symmetry.point_minus, axis)
        rot_plus_radial = radial_vector_to_axis(symmetry.rotated_point_plus, axis)
        rot_minus_radial = radial_vector_to_axis(symmetry.rotated_point_minus, axis)

        lines.append("Region-defined 2-fold symmetry axis:")
        lines.append(f"  regions                         : {symmetry.regions}")
        lines.append(f"  region source                   : {symmetry.region_source}")
        lines.append(f"  symmetry atoms                  : {', '.join(symmetry.atom_names)}")
        lines.append(f"  paired residues                 : {len(symmetry.pairs)}")
        lines.append("  pairing                         : first range vs reversed second range")
        lines.append(f"  center on helical axis          : {format_vec(symmetry.center)}")
        lines.append(f"  symmetry direction (unit)       : {format_vec(symmetry.direction)}")
        lines.append(f"  90-deg plane direction (unit)   : {format_vec(symmetry.perpendicular_direction)}")
        lines.append(f"  point radius from helical axis  : {symmetry.radius:.4f} A")
        lines.append("")
        lines.append("  points on the 2-fold symmetry axis:")
        lines.append(f"    + point                       : {format_vec(symmetry.point_plus)}")
        lines.append(
            f"      radial unit vector          : "
            f"{format_vec(plus_radial.unit if plus_radial.unit is not None else symmetry.direction)}"
        )
        lines.append(f"    - point                       : {format_vec(symmetry.point_minus)}")
        lines.append(
            f"      radial unit vector          : "
            f"{format_vec(minus_radial.unit if minus_radial.unit is not None else -symmetry.direction)}"
        )
        lines.append("")
        lines.append("  points rotated +90 deg around the helical axis:")
        lines.append(f"    + rotated point               : {format_vec(symmetry.rotated_point_plus)}")
        lines.append(
            f"      radial unit vector          : "
            f"{format_vec(rot_plus_radial.unit if rot_plus_radial.unit is not None else symmetry.perpendicular_direction)}"
        )
        lines.append(f"    - rotated point               : {format_vec(symmetry.rotated_point_minus)}")
        lines.append(
            f"      radial unit vector          : "
            f"{format_vec(rot_minus_radial.unit if rot_minus_radial.unit is not None else -symmetry.perpendicular_direction)}"
        )
        lines.append("")
        lines.append("  fit quality:")
        lines.append(f"    midpoint-to-axis RMSD         : {symmetry.midpoint_line_rmsd:.4f} A")
        lines.append(f"    180-deg paired-center RMSD    : {symmetry.rotation_rmsd:.4f} A")
        lines.append("    RMSD uses residue centers from the selected symmetry atoms; sequence is ignored.")

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
    axis_range: Optional[str] = None

    # Custom axis (if provided)
    axis_point: Optional[str] = None  # xyz triplet string
    axis_vector: Optional[str] = None  # xyz triplet string

    # Points (xyz triplet or atom spec)
    point1: Optional[str] = None
    point2: Optional[str] = None

    # Optional 2-fold symmetry axis from paired residue regions
    symmetry_enabled: bool = False
    symmetry_regions: Optional[str] = None
    symmetry_atoms: Optional[List[str]] = None
    symmetry_point_radius: Optional[float] = None

    out_bild: Path = Path("axis_vectors.bild")

    # Visual params
    axis_radius: float = 1.0
    vector_radius: float = 1.0
    sphere_radius: float = 1.25
    axis_margin: float = 5.0


def _needs_pdb_for_point(spec: str) -> bool:
    return _try_parse_xyz(spec) is None


def _has_text(value: Optional[str]) -> bool:
    return value is not None and bool(value.strip())


def run(cfg: RunConfig) -> str:
    use_custom_axis = (cfg.axis_point is not None) or (cfg.axis_vector is not None)
    if use_custom_axis and (cfg.axis_point is None or cfg.axis_vector is None):
        raise ValueError("Custom axis requires BOTH axis_point and axis_vector")
    if use_custom_axis and _has_text(cfg.axis_range):
        raise ValueError("axis_range applies only when fitting the axis from a PDB")

    axis_ranges = parse_axis_ranges(cfg.axis_range or "") if _has_text(cfg.axis_range) else None

    has_point1 = _has_text(cfg.point1)
    has_point2 = _has_text(cfg.point2)
    has_symmetry = bool(cfg.symmetry_enabled) or cfg.symmetry_regions is not None
    if has_point1 != has_point2:
        raise ValueError("Provide both point1 and point2, or leave both blank.")
    if not has_point1 and not has_symmetry:
        raise ValueError("Provide point1/point2 and/or enable symmetry-axis output.")

    need_pdb = False
    if not use_custom_axis:
        need_pdb = True
    else:
        if has_point1 and (
            _needs_pdb_for_point(cfg.point1 or "") or _needs_pdb_for_point(cfg.point2 or "")
        ):
            need_pdb = True
    if has_symmetry:
        need_pdb = True

    atoms: Optional[List[Atom]] = None
    atom_index: Optional[Dict[Tuple[str, int, str], Atom]] = None

    pdb_path: Optional[Path] = None
    if need_pdb:
        if cfg.pdb is None:
            raise ValueError("A PDB file is required for this run (axis fit from PDB, atom points, and/or symmetry output).")
        pdb_path = cfg.pdb
        atoms = load_pdb_atoms(pdb_path)
        atom_index = build_atom_index(atoms)

    # Parse points (now that we may have atom_index)
    p1: Optional[np.ndarray] = None
    p2: Optional[np.ndarray] = None
    if has_point1:
        p1 = parse_point_spec(cfg.point1 or "", atom_index)
        p2 = parse_point_spec(cfg.point2 or "", atom_index)

    symmetry_pairs_ref_points: List[np.ndarray] = []
    symmetry_atoms = cfg.symmetry_atoms or ["C1'"]
    if has_symmetry:
        assert atoms is not None
        symmetry_pairs, _effective_regions, _region_source = build_symmetry_pairs(
            atoms, cfg.symmetry_regions, symmetry_atoms
        )
        for pair in symmetry_pairs:
            symmetry_pairs_ref_points.extend([pair.first_center, pair.second_center, pair.midpoint])

    # Build axis
    if use_custom_axis:
        axis_p = parse_xyz_triplet(cfg.axis_point or "", what="axis point")
        axis_v = parse_xyz_triplet(cfg.axis_vector or "", what="axis vector")
        ref_points: List[np.ndarray] = []
        if p1 is not None:
            ref_points.append(p1)
        if p2 is not None:
            ref_points.append(p2)
        ref_points.extend(symmetry_pairs_ref_points)
        axis = axis_from_point_and_vector(axis_p, axis_v, ref_points=ref_points)
        axis_source = "Custom (user-defined point + vector)"
    else:
        assert atoms is not None
        axis = fit_axis_pca(atoms, cfg.axis_atoms, residue_ranges=axis_ranges)
        axis_source = (
            "PDB fit (PCA on axis atoms in selected residue ranges)"
            if axis_ranges is not None
            else "PDB fit (PCA on axis atoms)"
        )

    r1: Optional[RadialResult] = None
    r2: Optional[RadialResult] = None
    if p1 is not None and p2 is not None:
        r1 = radial_vector_to_axis(p1, axis, spec=cfg.point1 or "")
        r2 = radial_vector_to_axis(p2, axis, spec=cfg.point2 or "")

    symmetry: Optional[SymmetryAxisResult] = None
    if has_symmetry:
        assert atoms is not None
        symmetry = determine_twofold_symmetry_axis(
            atoms,
            axis,
            cfg.symmetry_regions,
            symmetry_atoms,
            point_radius=cfg.symmetry_point_radius,
        )

    # Write BILD
    write_bild(
        cfg.out_bild,
        axis,
        r1,
        r2,
        symmetry=symmetry,
        axis_radius=cfg.axis_radius,
        vector_radius=cfg.vector_radius,
        sphere_radius=cfg.sphere_radius,
        axis_margin=cfg.axis_margin,
    )

    report = make_report(
        pdb_path,
        axis_source,
        cfg.axis_atoms,
        cfg.axis_range,
        axis,
        r1,
        r2,
        symmetry,
    )
    report += f"\nBILD written: {cfg.out_bild}\n"
    return report


# -----------------------------
# GUI
# -----------------------------


def _set_optional_window_icon(window, tk_module) -> None:
    icon_path = Path(__file__).resolve().parents[1] / "assets" / "bnp_na_icon.png"
    if not icon_path.exists():
        return
    try:
        icon = tk_module.PhotoImage(file=str(icon_path))
        window._bnp_na_window_icon = icon
        window.iconphoto(True, icon)
    except Exception:
        pass


def launch_gui(parent=None, log_callback=None) -> None:
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
    _set_optional_window_icon(root, tk)
    root.title("angle_helical_axisV2_2.py")
    root.geometry("900x800")
    root.minsize(760, 620)

    # Variables
    var_pdb = tk.StringVar(value="")
    var_out = tk.StringVar(value="")

    # Axis source
    var_axis_source = tk.StringVar(value="Fit from PDB")  # or "Custom"
    var_axis_atoms = tk.StringVar(value="C1'")
    var_axis_range = tk.StringVar(value="")

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

    # Optional 2-fold symmetry-axis vars
    var_sym_enabled = tk.StringVar(value="Off")
    var_sym_regions = tk.StringVar(value="")
    var_sym_atoms = tk.StringVar(value="C1'")
    var_sym_radius = tk.StringVar(value=f"{DEFAULT_SYMMETRY_POINT_RADIUS:.1f}")

    # Radii defaults (per request)
    var_axis_r = tk.StringVar(value="1.0")
    var_vec_r = tk.StringVar(value="1.0")
    var_sph_r = tk.StringVar(value="1.25")
    var_axis_margin = tk.StringVar(value="5.0")

    def browse_pdb() -> None:
        fn = filedialog.askopenfilename(
            **_file_dialog_options("Select PDB file", PDB_OPEN_FILETYPES)
        )
        if fn:
            var_pdb.set(fn)
            # Suggest output
            base = os.path.splitext(fn)[0]
            var_out.set(base + "_axis_vectors.bild")

    def browse_out() -> None:
        fn = filedialog.asksaveasfilename(
            **_file_dialog_options(
                "Save BILD file",
                BILD_SAVE_FILETYPES,
                defaultextension=".bild",
            )
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

        sym_enabled = var_sym_enabled.get() == "On"
        sym_regions = var_sym_regions.get().strip() if sym_enabled else ""
        sym_atoms = parse_axis_atoms(var_sym_atoms.get())
        sym_radius: Optional[float] = None
        if sym_enabled and var_sym_radius.get().strip():
            try:
                sym_radius = float(var_sym_radius.get())
            except Exception:
                messagebox.showerror("Error", "2-fold symmetry point radius must be a number or left blank.")
                return

        if bool(p1_spec) != bool(p2_spec):
            messagebox.showerror("Error", "Please provide both Point 1 and Point 2, or leave both blank.")
            return
        if not p1_spec and not sym_enabled:
            messagebox.showerror("Error", "Please provide Point 1/Point 2 and/or turn on 2-fold symmetry output.")
            return

        # Validate xyz fields if needed
        if p1_spec and var_p1_mode.get() == "XYZ":
            err = _validate_xyz_fields(var_p1_x.get(), var_p1_y.get(), var_p1_z.get(), "Point 1")
            if err:
                messagebox.showerror("Error", err)
                return
        if p2_spec and var_p2_mode.get() == "XYZ":
            err = _validate_xyz_fields(var_p2_x.get(), var_p2_y.get(), var_p2_z.get(), "Point 2")
            if err:
                messagebox.showerror("Error", err)
                return
        if sym_enabled and sym_regions:
            try:
                parse_symmetry_regions(sym_regions)
            except Exception as e:
                messagebox.showerror("Error", str(e))
                return

        out = var_out.get().strip()
        if not out:
            messagebox.showerror("Error", "Please choose an output .bild path")
            return

        # Axis config
        axis_point = None
        axis_vector = None
        axis_atoms = parse_axis_atoms(var_axis_atoms.get())
        axis_range = (
            var_axis_range.get().strip() if axis_source == "Fit from PDB" else ""
        )

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
        elif axis_range:
            try:
                parse_axis_ranges(axis_range)
            except Exception as e:
                messagebox.showerror("Error", str(e))
                return

        # PDB requirement
        pdb = var_pdb.get().strip()
        pdb_path: Optional[Path] = None
        need_pdb = False
        if axis_source == "Fit from PDB":
            need_pdb = True
        if p1_spec and (var_p1_mode.get() == "Atom" or var_p2_mode.get() == "Atom"):
            need_pdb = True
        if sym_enabled:
            need_pdb = True

        if need_pdb:
            if not pdb or not os.path.isfile(pdb):
                messagebox.showerror(
                    "Error",
                    "A valid PDB file is required (axis fitted from PDB, points specified as atoms, and/or symmetry regions).",
                )
                return
            pdb_path = Path(pdb)

        try:
            cfg = RunConfig(
                pdb=pdb_path,
                axis_atoms=axis_atoms,
                axis_range=axis_range or None,
                axis_point=axis_point,
                axis_vector=axis_vector,
                point1=p1_spec or None,
                point2=p2_spec or None,
                symmetry_enabled=sym_enabled,
                symmetry_regions=sym_regions or None,
                symmetry_atoms=sym_atoms,
                symmetry_point_radius=sym_radius,
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
        if log_callback is not None:
            cli_args = ["python3", "bnp_na_lib/angle_helical_axisV2_2.py"]
            if pdb_path is not None:
                cli_args.extend(["--pdb", str(pdb_path)])
            if axis_source == "Custom axis":
                cli_args.extend(["--axis-point", axis_point or "", "--axis-vector", axis_vector or ""])
            else:
                cli_args.extend(["--axis-atoms", " ".join(axis_atoms)])
                if axis_range:
                    cli_args.extend(["--axis_range", axis_range])
            if p1_spec and p2_spec:
                cli_args.extend(["--point1", p1_spec, "--point2", p2_spec])
            if sym_enabled:
                cli_args.append("--symmetry-axis")
                if sym_regions:
                    cli_args.extend(["--symmetry-regions", sym_regions])
                cli_args.extend(["--symmetry-atoms", " ".join(sym_atoms)])
                if sym_radius is not None:
                    cli_args.extend(["--symmetry-point-radius", var_sym_radius.get()])
            cli_args.extend(
                [
                    "--out-bild",
                    out,
                    "--axis-margin",
                    var_axis_margin.get(),
                    "--axis-radius",
                    var_axis_r.get(),
                    "--vector-radius",
                    var_vec_r.get(),
                    "--sphere-radius",
                    var_sph_r.get(),
                ]
            )
            main_log = "\n".join(
                [
                    "=== Helical-axis angle tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    report,
                ]
            )
            try:
                log_callback(main_log)
            except Exception:
                pass

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

    # Axis residue-range row (PDB fit)
    tk.Label(frm, text="Axis residue ranges (PDB fit only)", anchor="w").grid(
        row=3, column=0, sticky="w"
    )
    ent_axis_range = tk.Entry(frm, textvariable=var_axis_range, width=32)
    ent_axis_range.grid(row=3, column=1, sticky="w", padx=(5, 5))

    # Custom axis frame
    frm_custom = tk.Frame(frm)
    frm_custom.grid(row=4, column=0, columnspan=3, sticky="we", pady=(2, 2))

    tk.Label(frm_custom, text="Custom axis point (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_custom, textvariable=var_axis_px, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_custom, textvariable=var_axis_py, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_custom, textvariable=var_axis_pz, width=8).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_custom, text="Custom axis vector (vx vy vz)").grid(row=1, column=0, sticky="w")
    tk.Entry(frm_custom, textvariable=var_axis_vx, width=8).grid(row=1, column=1, padx=(5, 2))
    tk.Entry(frm_custom, textvariable=var_axis_vy, width=8).grid(row=1, column=2, padx=(2, 2))
    tk.Entry(frm_custom, textvariable=var_axis_vz, width=8).grid(row=1, column=3, padx=(2, 2))

    # Point 1 controls
    tk.Label(frm, text="Point 1 type", anchor="w").grid(row=5, column=0, sticky="w")
    tk.OptionMenu(frm, var_p1_mode, "Atom", "XYZ").grid(row=5, column=1, sticky="w", padx=(5, 5))

    frm_p1_atom = tk.Frame(frm)
    frm_p1_xyz = tk.Frame(frm)
    frm_p1_atom.grid(row=6, column=0, columnspan=3, sticky="we")
    frm_p1_xyz.grid(row=6, column=0, columnspan=3, sticky="we")

    tk.Label(frm_p1_atom, text="Point 1 atom (chain, resSeq, atom)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p1_atom, textvariable=var_p1_chain, width=5).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p1_atom, textvariable=var_p1_resseq, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p1_atom, textvariable=var_p1_atom, width=10).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_p1_xyz, text="Point 1 xyz (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p1_xyz, textvariable=var_p1_x, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p1_xyz, textvariable=var_p1_y, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p1_xyz, textvariable=var_p1_z, width=8).grid(row=0, column=3, padx=(2, 2))

    # Point 2 controls
    tk.Label(frm, text="Point 2 type", anchor="w").grid(row=7, column=0, sticky="w")
    tk.OptionMenu(frm, var_p2_mode, "Atom", "XYZ").grid(row=7, column=1, sticky="w", padx=(5, 5))

    frm_p2_atom = tk.Frame(frm)
    frm_p2_xyz = tk.Frame(frm)
    frm_p2_atom.grid(row=8, column=0, columnspan=3, sticky="we")
    frm_p2_xyz.grid(row=8, column=0, columnspan=3, sticky="we")

    tk.Label(frm_p2_atom, text="Point 2 atom (chain, resSeq, atom)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p2_atom, textvariable=var_p2_chain, width=5).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p2_atom, textvariable=var_p2_resseq, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p2_atom, textvariable=var_p2_atom, width=10).grid(row=0, column=3, padx=(2, 2))

    tk.Label(frm_p2_xyz, text="Point 2 xyz (x y z)").grid(row=0, column=0, sticky="w")
    tk.Entry(frm_p2_xyz, textvariable=var_p2_x, width=8).grid(row=0, column=1, padx=(5, 2))
    tk.Entry(frm_p2_xyz, textvariable=var_p2_y, width=8).grid(row=0, column=2, padx=(2, 2))
    tk.Entry(frm_p2_xyz, textvariable=var_p2_z, width=8).grid(row=0, column=3, padx=(2, 2))

    # Optional 2-fold symmetry-axis controls
    tk.Label(frm, text="2-fold symmetry output", anchor="w").grid(row=9, column=0, sticky="w")
    frm_sym_toggle = tk.Frame(frm)
    frm_sym_toggle.grid(row=9, column=1, sticky="w", padx=(5, 5))
    tk.Radiobutton(frm_sym_toggle, text="Off", variable=var_sym_enabled, value="Off").grid(row=0, column=0, sticky="w")
    tk.Radiobutton(frm_sym_toggle, text="On", variable=var_sym_enabled, value="On").grid(row=0, column=1, sticky="w", padx=(12, 0))

    lbl_sym_regions = tk.Label(frm, text="2-fold symmetry regions", anchor="w")
    lbl_sym_regions.grid(row=10, column=0, sticky="w")
    ent_sym_regions = tk.Entry(frm, textvariable=var_sym_regions, width=32)
    ent_sym_regions.grid(row=10, column=1, sticky="w", padx=(5, 5))

    lbl_sym_atoms = tk.Label(frm, text="2-fold symmetry atoms", anchor="w")
    lbl_sym_atoms.grid(row=11, column=0, sticky="w")
    ent_sym_atoms = tk.Entry(frm, textvariable=var_sym_atoms, width=20)
    ent_sym_atoms.grid(row=11, column=1, sticky="w", padx=(5, 5))

    lbl_sym_radius = tk.Label(frm, text="2-fold point radius (A)", anchor="w")
    lbl_sym_radius.grid(row=12, column=0, sticky="w")
    ent_sym_radius = tk.Entry(frm, textvariable=var_sym_radius, width=10)
    ent_sym_radius.grid(row=12, column=1, sticky="w", padx=(5, 5))
    sym_dependent_widgets = [
        lbl_sym_regions,
        ent_sym_regions,
        lbl_sym_atoms,
        ent_sym_atoms,
        lbl_sym_radius,
        ent_sym_radius,
    ]

    # Output
    add_entry_row(13, "Output .bild", var_out, browse_out)

    # Axis drawing margin and radii
    tk.Label(frm, text="Axis drawing margin", anchor="w").grid(row=14, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_axis_margin, width=10).grid(row=14, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Axis radius", anchor="w").grid(row=15, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_axis_r, width=10).grid(row=15, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Vector radius", anchor="w").grid(row=16, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_vec_r, width=10).grid(row=16, column=1, sticky="w", padx=(5, 5))

    tk.Label(frm, text="Sphere radius", anchor="w").grid(row=17, column=0, sticky="w")
    tk.Entry(frm, textvariable=var_sph_r, width=10).grid(row=17, column=1, sticky="w", padx=(5, 5))

    tk.Button(frm, text="Run", command=do_run).grid(row=18, column=0, pady=(8, 8), sticky="w")

    # Output text
    text = ScrolledText(frm, height=18, width=90)
    text.grid(row=19, column=0, columnspan=3, sticky="nsew", pady=(5, 0))

    frm.columnconfigure(1, weight=1)
    frm.rowconfigure(19, weight=1)

    def refresh_visibility(*_args) -> None:
        # Axis controls
        if var_axis_source.get() == "Custom axis":
            frm_custom.grid()
            ent_axis_atoms.configure(state="disabled")
            ent_axis_range.configure(state="disabled")
        else:
            frm_custom.grid_remove()
            ent_axis_atoms.configure(state="normal")
            ent_axis_range.configure(state="normal")

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

        sym_state = "normal" if var_sym_enabled.get() == "On" else "disabled"
        for widget in sym_dependent_widgets:
            try:
                widget.configure(state=sym_state)
            except Exception:
                pass

    var_axis_source.trace_add("write", refresh_visibility)
    var_p1_mode.trace_add("write", refresh_visibility)
    var_p2_mode.trace_add("write", refresh_visibility)
    var_sym_enabled.trace_add("write", refresh_visibility)

    refresh_visibility()

    text.insert(
        "1.0",
        "Point inputs:\n"
        "  - Atom: chain/resSeq/atom (atom names accept prime or star, e.g. C1' or C1*)\n"
        "  - XYZ : numeric x y z\n"
        "  - Point 1/Point 2 can be left blank when 2-fold symmetry regions are provided\n"
        "Axis:\n"
        "  - Fit from PDB: uses --axis-atoms (default C1')\n"
        "  - Optional axis ranges select fit residues, e.g. A1-A35,B60-B26\n"
        "  - The first axis range's written start-to-end order sets the positive direction\n"
        "  - Custom axis: provide axis point + axis vector (vector is normalized automatically)\n"
        "  - Axis drawing margin controls how far the displayed axis extends beyond the selected points/fit range\n"
        "2-fold symmetry axis:\n"
        "  - Default is Off; choose On to calculate and draw symmetry-axis information\n"
        "  - Regions use two ranges such as A1-A35, B26-B60\n"
        "  - If regions are left blank while On, the whole two-chain model is used\n"
        "  - The four symmetry-related points default to 15.0 A from the helical axis\n"
        "  - The first range is paired with the reversed second range\n"
        "  - Symmetry atoms define residue centers; use all or * to average all atoms in each residue\n",
    )

    if owns_mainloop:
        root.mainloop()
    else:
        try:
            root.after(50, lambda: (root.deiconify(), root.lift(), root.focus_force()))
        except Exception:
            pass


# -----------------------------
# CLI
# -----------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="angle_helical_axisV2_2.py",
        description=(
            "Fit (or define) a helical axis, compute radial-vector angles for two points, "
            "and/or determine a region-defined 2-fold symmetry axis."
        ),
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
        "--axis-range",
        "--axis_range",
        dest="axis_range",
        type=str,
        default=None,
        help=(
            'Residue range(s) used for the PDB axis fit, e.g. "A1-A35,B60-B26". '
            "The first range's written start-to-end order sets the positive axis direction."
        ),
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
        "--symmetry-axis",
        action="store_true",
        help="Calculate 2-fold symmetry-axis information. If --symmetry-regions is omitted or blank, use the whole two-chain model.",
    )
    p.add_argument(
        "--symmetry-regions",
        type=str,
        default=None,
        help='Two residue ranges defining a 2-fold axis, e.g. "A1-A35, B26-B60". Providing this also enables symmetry-axis output.',
    )
    p.add_argument(
        "--symmetry-atoms",
        type=str,
        default="C1'",
        help="Atom name(s) used for region residue centers. Use all or * for all atoms. Default: C1'",
    )
    p.add_argument(
        "--symmetry-point-radius",
        type=float,
        default=DEFAULT_SYMMETRY_POINT_RADIUS,
        help=f"Radius from the helical axis for reported/drawn 2-fold and 90-degree points. Default: {DEFAULT_SYMMETRY_POINT_RADIUS:.1f} A.",
    )

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

    has_point1 = _has_text(args.point1)
    has_point2 = _has_text(args.point2)
    has_symmetry = bool(args.symmetry_axis) or args.symmetry_regions is not None
    if has_point1 != has_point2:
        parser.error("Provide both --point1 and --point2, or leave both blank.")
    if not has_point1 and not has_symmetry:
        parser.error("CLI mode requires --point1/--point2 and/or --symmetry-axis/--symmetry-regions (or use --gui).")

    use_custom_axis = (args.axis_point is not None) or (args.axis_vector is not None)
    if use_custom_axis and (args.axis_point is None or args.axis_vector is None):
        parser.error("Custom axis requires BOTH --axis-point and --axis-vector")
    if use_custom_axis and _has_text(args.axis_range):
        parser.error("--axis-range/--axis_range applies only to a PDB-fitted axis")
    if _has_text(args.axis_range):
        try:
            parse_axis_ranges(args.axis_range)
        except Exception as exc:
            parser.error(str(exc))
    if has_symmetry and _has_text(args.symmetry_regions):
        try:
            parse_symmetry_regions(args.symmetry_regions)
        except Exception as exc:
            parser.error(str(exc))

    pdb_path: Optional[Path] = Path(args.pdb) if args.pdb else None

    # Determine whether a PDB is required
    need_pdb = False
    if not use_custom_axis:
        need_pdb = True
    else:
        if has_point1 and (
            _needs_pdb_for_point(args.point1 or "") or _needs_pdb_for_point(args.point2 or "")
        ):
            need_pdb = True
    if has_symmetry:
        need_pdb = True

    if need_pdb:
        if pdb_path is None:
            parser.error("This run requires --pdb (axis fit from PDB, atom points, and/or symmetry regions).")
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
        axis_range=args.axis_range,
        axis_point=args.axis_point,
        axis_vector=args.axis_vector,
        point1=args.point1 if has_point1 else None,
        point2=args.point2 if has_point2 else None,
        symmetry_enabled=has_symmetry,
        symmetry_regions=args.symmetry_regions if has_symmetry else None,
        symmetry_atoms=parse_axis_atoms(args.symmetry_atoms),
        symmetry_point_radius=args.symmetry_point_radius,
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
