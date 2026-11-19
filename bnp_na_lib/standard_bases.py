"""Standard nucleotide bases and base grafting for bnp_na.

The four DNA bases below are in the standard base reference frame (Olson et
al., J. Mol. Biol. 313, 229-237, 2001). They were taken from x3dna-dssr
``rebuild`` output for a base pair whose six base-pair parameters are all zero,
so they match the base geometry of the B-DNA, A-DNA and A-RNA models that
bnp_na builds with DSSR.

``graft_base`` replaces the base of one nucleotide by another base type: the
standard copy of the current base is fitted to the observed ring atoms, the
same way 3DNA defines a base reference frame, and the new base is placed in
that frame. The sugar-phosphate backbone, including C1', is kept as it is.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np


# Base atoms in the standard reference frame, in the atom order DSSR writes.
# C1' is listed first for reference; graft_base keeps the residue's own C1'.
STANDARD_BASES: Dict[str, Tuple[Tuple[str, float, float, float], ...]] = {
    "A": (
        ("C1'", -2.484, 5.344, 0.000),
        ("N9", -1.295, 4.497, 0.000),
        ("C8", 0.020, 4.897, 0.000),
        ("N7", 0.874, 3.903, 0.000),
        ("C5", 0.069, 2.771, 0.000),
        ("C6", 0.368, 1.398, 0.000),
        ("N6", 1.610, 0.910, 0.000),
        ("N1", -0.668, 0.531, 0.000),
        ("C2", -1.913, 1.021, 0.000),
        ("N3", -2.322, 2.288, 0.000),
        ("C4", -1.270, 3.123, 0.000),
    ),
    "C": (
        ("C1'", -2.482, 5.400, 0.000),
        ("N1", -1.289, 4.541, 0.000),
        ("C2", -1.475, 3.157, 0.000),
        ("O2", -2.630, 2.707, 0.001),
        ("N3", -0.393, 2.344, 0.000),
        ("C4", 0.834, 2.869, 0.000),
        ("N4", 1.873, 2.029, 0.001),
        ("C5", 1.052, 4.276, 0.000),
        ("C6", -0.027, 5.068, 0.000),
    ),
    "G": (
        ("C1'", -2.482, 5.397, 0.000),
        ("N9", -1.293, 4.550, 0.000),
        ("C8", 0.019, 4.962, 0.000),
        ("N7", 0.867, 3.970, 0.000),
        ("C5", 0.069, 2.833, 0.000),
        ("C6", 0.423, 1.460, 0.000),
        ("O6", 1.553, 0.956, 0.000),
        ("N1", -0.701, 0.640, 0.000),
        ("C2", -2.000, 1.085, 0.000),
        ("N2", -2.949, 0.136, -0.001),
        ("N3", -2.344, 2.362, 0.001),
        ("C4", -1.268, 3.176, 0.000),
    ),
    "T": (
        ("C1'", -2.486, 5.352, 0.000),
        ("N1", -1.288, 4.499, 0.000),
        ("C2", -1.465, 3.134, 0.000),
        ("O2", -2.564, 2.606, 0.000),
        ("N3", -0.300, 2.407, 0.000),
        ("C4", 0.991, 2.898, 0.000),
        ("O4", 1.942, 2.121, 0.000),
        ("C5", 1.102, 4.339, 0.000),
        ("C7", 2.462, 4.963, 0.001),
        ("C6", -0.028, 5.057, 0.000),
    ),
}

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CTU")
PURINE_RING_ATOMS = ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4")
PYRIMIDINE_RING_ATOMS = ("N1", "C2", "N3", "C4", "C5", "C6")

# Sugar-phosphate atoms of a DNA nucleotide, in the order bnp_na writes them.
DNA_BACKBONE_ATOMS = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")

DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
DNA_RESIDUE_NAMES = {"A": "DA", "C": "DC", "G": "DG", "T": "DT"}

# Every base atom name of the four DNA bases (plus the old thymine methyl name).
ALL_BASE_ATOM_NAMES = frozenset(
    name for atoms in STANDARD_BASES.values() for name, *_xyz in atoms if name != "C1'"
) | {"C5M"}


class BaseGraftError(ValueError):
    """Raised when a base cannot be identified or replaced."""


def base_code(residue_name: str) -> str:
    """Return A/C/G/T/U for residue names such as DA, DC, A, ADE or DT."""
    name = (residue_name or "").strip().upper()
    aliases = {"ADE": "A", "CYT": "C", "GUA": "G", "THY": "T", "URA": "U", "URI": "U"}
    if name in aliases:
        return aliases[name]
    for ch in reversed(name):
        if ch in "ACGTU":
            return ch
    raise BaseGraftError(f"Cannot identify the base of residue {residue_name!r}.")


def dna_reverse_complement(seq: str) -> str:
    try:
        return "".join(DNA_COMPLEMENT[base] for base in reversed(seq.upper()))
    except KeyError as exc:
        raise BaseGraftError(f"Cannot complement base {exc.args[0]!r}; only A/C/G/T are allowed.") from exc


def ring_atoms(base: str) -> Tuple[str, ...]:
    code = base.upper()
    if code in PURINES:
        return PURINE_RING_ATOMS
    if code in PYRIMIDINES:
        return PYRIMIDINE_RING_ATOMS
    raise BaseGraftError(f"Unsupported base {base!r}.")


def standard_base_atoms(base: str, include_c1: bool = False) -> List[Tuple[str, np.ndarray]]:
    code = base.upper()
    if code not in STANDARD_BASES:
        raise BaseGraftError(f"No standard geometry for base {base!r}; only A/C/G/T are available.")
    return [
        (name, np.array([x, y, z], dtype=float))
        for name, x, y, z in STANDARD_BASES[code]
        if include_c1 or name != "C1'"
    ]


def kabsch(moving: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return (R, t, rmsd) minimizing |R @ moving + t - target| over paired points."""
    moving = np.asarray(moving, dtype=float)
    target = np.asarray(target, dtype=float)
    if moving.shape != target.shape or moving.ndim != 2 or moving.shape[1] != 3 or len(moving) < 3:
        raise ValueError(f"kabsch needs two matching (N>=3, 3) arrays; got {moving.shape} and {target.shape}")
    pc = moving.mean(axis=0)
    qc = target.mean(axis=0)
    h = (moving - pc).T @ (target - qc)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = qc - rot @ pc
    fitted = moving @ rot.T + trans
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rot, trans, rmsd


def _rotation_about_origin(moving: np.ndarray, target: np.ndarray) -> np.ndarray:
    h = moving.T @ target
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    return vt.T @ np.diag([1.0, 1.0, d]) @ u.T


def base_frame(
    atoms: Mapping[str, Sequence[float]],
    base: str,
    pivot_on_c1: bool = True,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Fit the standard base to a residue's ring atoms.

    Returns (R, t, rmsd) that map standard-frame coordinates onto the residue,
    i.e. observed = R @ standard + t, and the ring-atom RMSD of that map.

    With ``pivot_on_c1`` (and a C1' atom present) the standard C1' is put on the
    residue's C1' and only the rotation about it is fitted to the ring atoms.
    That keeps the standard C1'-N glycosidic bond (1.47 A) and its angles even
    when the residue's own base is distorted, as in low-resolution models; for
    an ideal base it gives the same frame as the plain ring fit.
    """
    names = ring_atoms(base)
    missing = [name for name in names if name not in atoms]
    if missing:
        raise BaseGraftError(f"Base {base} is missing ring atom(s): {', '.join(missing)}")
    standard = dict(standard_base_atoms(base, include_c1=True))
    moving = np.array([standard[name] for name in names], dtype=float)
    target = np.array([np.asarray(atoms[name], dtype=float) for name in names], dtype=float)
    if not (pivot_on_c1 and "C1'" in atoms):
        return kabsch(moving, target)
    c1_std = standard["C1'"]
    c1_obs = np.asarray(atoms["C1'"], dtype=float)
    rot = _rotation_about_origin(moving - c1_std, target - c1_obs)
    trans = c1_obs - rot @ c1_std
    fitted = moving @ rot.T + trans
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rot, trans, rmsd


def placed_base_atoms(base: str, rot: np.ndarray, trans: np.ndarray) -> List[Tuple[str, np.ndarray]]:
    """Standard base atoms (no C1') moved into the frame (rot, trans)."""
    return [(name, rot @ xyz + trans) for name, xyz in standard_base_atoms(base)]


def standard_c1(base: str) -> np.ndarray:
    for name, x, y, z in STANDARD_BASES[base.upper()]:
        if name == "C1'":
            return np.array([x, y, z], dtype=float)
    raise BaseGraftError(f"No standard C1' for base {base!r}.")


def is_base_atom(name: str) -> bool:
    """True for base heavy atoms and base hydrogens (sugar names carry a prime)."""
    if name in ALL_BASE_ATOM_NAMES:
        return True
    return name.startswith("H") and "'" not in name and "*" not in name


def glycosidic_atom(base: str) -> str:
    return "N9" if base.upper() in PURINES else "N1"


def rotate_base_about_glycosidic_bond(
    atoms: Sequence[Tuple[str, np.ndarray]],
    base: str,
    delta_deg: float,
) -> List[Tuple[str, np.ndarray]]:
    """Turn the base atoms about the C1'-N bond; chi changes by ``delta_deg``."""
    names = dict(atoms)
    c1 = np.asarray(names["C1'"], dtype=float)
    axis = np.asarray(names[glycosidic_atom(base)], dtype=float) - c1
    axis = axis / np.linalg.norm(axis)
    angle = np.radians(delta_deg)
    k = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    rot = np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * (k @ k)
    return [
        (name, rot @ (np.asarray(xyz, dtype=float) - c1) + c1 if is_base_atom(name) else np.asarray(xyz, dtype=float))
        for name, xyz in atoms
    ]


def torsion_deg(p0: Sequence[float], p1: Sequence[float], p2: Sequence[float], p3: Sequence[float]) -> float:
    a0, a1, a2, a3 = (np.asarray(p, dtype=float) for p in (p0, p1, p2, p3))
    b0 = a0 - a1
    b1 = a2 - a1
    b2 = a3 - a2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def chi_deg(atoms: Mapping[str, Sequence[float]], base: str) -> float:
    """Glycosidic torsion O4'-C1'-N9-C4 (purines) or O4'-C1'-N1-C2 (pyrimidines)."""
    second = "C4" if base.upper() in PURINES else "C2"
    return torsion_deg(atoms["O4'"], atoms["C1'"], atoms[glycosidic_atom(base)], atoms[second])


def graft_base(
    atoms: Mapping[str, Sequence[float]],
    current_base: str,
    new_base: str,
) -> Tuple[List[Tuple[str, np.ndarray]], float]:
    """Replace the base of one nucleotide.

    ``atoms`` maps atom names to coordinates for a single residue. Returns the
    new atom list (sugar-phosphate atoms first, in the input order, then the new
    base in standard order) and the RMSD of the ring fit that defined the frame.
    Base hydrogens of the old base are dropped; the new base has none.
    """
    rot, trans, rmsd = base_frame(atoms, current_base)
    if "C1'" in atoms:
        # Put the new base's own standard C1' on the residue's C1', so the new
        # base keeps its standard glycosidic bond length and angles.
        trans = np.asarray(atoms["C1'"], dtype=float) - rot @ standard_c1(new_base)
    kept = [
        (name, np.asarray(xyz, dtype=float))
        for name, xyz in atoms.items()
        if not is_base_atom(name)
    ]
    return kept + placed_base_atoms(new_base, rot, trans), rmsd
