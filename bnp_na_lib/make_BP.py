#!/usr/bin/env python3
"""make_BP.py

Build a multi-segment B-DNA / P-DNA construct connected by B-P junction cores.

This is the B-P counterpart of make_BZV2_3.py and follows the same recipe:

    B1.pdb P1.pdb B2.pdb ...      (or P1.pdb B1.pdb P2.pdb ...)

A fresh copy of the bundled B-P junction core (core_BP.py) is inserted between
every adjacent pair of helices. The core is fitted to the first helix end on
sugar atoms (C1', C2', C3', C4', O5', C5'), the next helix is fitted to the
core's other end, overlapping residues are trimmed, and the construct is
ligated into two continuous strands numbered 5' -> 3'.

The core holds two B-DNA base pairs, two transition levels and two P-DNA
levels. Unlike the B-Z core, it carries no sequence of its own:

  * B end: the helix's two terminal base pairs are replaced by the core's two
    B base pairs, which take over their bases (grafted in the standard base
    reference frame), so every B-DNA base pair is kept.
  * P end: the helix's two terminal P-DNA levels are rebuilt as the core's
    transition levels, which take over their bases. The core's two P levels
    only position the P helix (they match the 1PFI template exactly) and the
    helix keeps its own nucleotides from the third level on.

The output strands therefore carry exactly the input sequences. The two
transition levels on each junction are the P-DNA side's first two levels: the
backbone there runs from the B-DNA radius (P at about 9.5 A) into the P-DNA
core (P at 2.5 A). There is no experimental B-P junction structure; see
core_BP.py for how the transition was modeled.

Lengths: the transition nucleotides count as P-DNA, so a P-DNA of N levels
keeps N - 2 P-form levels with one junction and N - 4 between two junctions;
for k P-form levels, enter k + 2 (or k + 4). Hence a P-DNA needs at least 4
levels, or 7 with junctions at both ends (MIN_P_LEVELS). A B-DNA keeps every
base pair and needs at least 3 bp, or 5 with junctions at both ends (MIN_B_BP).

Helix types are detected from the geometry (P-DNA has its phosphates within
about 5 A of the helix axis; B-DNA near 9.5 A) unless --types is given, and
consecutive inputs must alternate between B and P. Use bnp_na-built inputs:
B-DNA from the default pipeline (minimized and regularized) matches the core's
B pairs exactly; other B-DNA fits within the reported RMSD. P-DNA built by
bnp_na always matches.

Axis handling follows make_BZ: codirectional (default) rotates each moving
helix so its axis is parallel to the previous helix axis; collinear also
shifts it onto the same line; none skips the correction. The core's B and P
parts are coaxial, so for bnp_na-built inputs the correction is close to zero.
B-DNA axes come from DSSR (auto falls back to the helix's screw axis); P-DNA
axes always come from the helix's own screw symmetry, because DSSR finds no
helix in P-DNA. Nothing is relaxed or refined: a build is rigid fitting,
trimming, base grafting and ligation, and takes well under a second plus any
DSSR axis calls.

Example:
    python3 bnp_na_lib/make_BP.py --out B_P_B.pdb B1.pdb P1.pdb B2.pdb
"""

from __future__ import annotations

if __package__ in (None, ""):  # run as a script: python3 bnp_na_lib/make_BP.py
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    __package__ = "bnp_na_lib"

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import core_BP
from .edit_pdb_atom import pdb_atom_record, rec2file
from .make_BZV2_3 import (
    AxisLineResult,
    ResidKey,
    Residue,
    Structure,
    _bp_centroid,
    _distance_min_O3P,
    apply_rigid_transform,
    apply_rotation_about_point,
    apply_translation,
    axis_line_distance,
    build_basepairs,
    build_residues,
    connected_components_keys,
    ensure_order_5to3,
    get_axis_line,
    kabsch,
    load_pdb_records,
    make_ter_record,
    normalize_axis_mode,
    order_path_component_keys,
    perpendicular_shift_to_collinear,
    pick_two_main_chains,
    raw_output_path,
    rotation_matrix_from_vectors,
    write_raw_aligned_pdb,
)
from .standard_bases import base_code, graft_base


DEFAULT_OUTPUT_NAME = "make_BP_out.pdb"
SUGAR_ATOMS: List[str] = list(core_BP.SUGAR_ATOMS)
N_TRANSITION = len(core_BP.TRANSITION_LEVELS)
# Residues each helix gives up at a junction end.
B_END_TRIM = len(core_BP.B_END_BASEPAIRS)   # replaced by the core's B pairs
P_END_TRIM = N_TRANSITION                   # rebuilt as the core's transition levels
# Minimum pairs (bp or levels) a helix needs, by number of junction ends.
# B gives up two bp per end and keeps at least one between them. P gives up two
# levels per end and keeps the next two, which the core is fitted on; with
# junctions at both ends it must keep three, or the two transitions reach into
# each other (a 6-level P-DNA between two B-DNAs puts two bases 0.4 A apart).
MIN_B_BP = {1: B_END_TRIM + 1, 2: 2 * B_END_TRIM + 1}
MIN_P_LEVELS = {1: N_TRANSITION + len(core_BP.P_END_LEVELS), 2: 2 * N_TRANSITION + 3}
P_HELIX_MAX_P_RADIUS = 5.0
# Fits within this RMSD of the best one count as equally good when choosing a helix end.
END_TIE_RMSD = 0.1


# ----------------------------------------------------------------------------
# Loading and helix types
# ----------------------------------------------------------------------------

def load_helix(path: Path, seg_id: str) -> Structure:
    records = load_pdb_records(path)
    residues = build_residues(records, seg_id=seg_id)
    chains = pick_two_main_chains(residues)
    basepairs = build_basepairs(residues, chains)
    return Structure(name=path.name, seg_id=seg_id, records=records, residues=residues, chains=chains,
                     basepairs=basepairs)


def new_core(seg_id: str) -> Structure:
    records = core_BP.core_records()
    residues = build_residues(records, seg_id=seg_id)
    chains = [core_BP.STRAND_E, core_BP.STRAND_F]
    # Physical order along the core: B outer pair first, P outer level last.
    pairs = (list(reversed(core_BP.B_END_BASEPAIRS)) + list(core_BP.TRANSITION_LEVELS)
             + list(core_BP.P_END_LEVELS))
    basepairs = [((seg_id, a[0], a[1]), (seg_id, b[0], b[1])) for a, b in pairs]
    return Structure(name=f"core_BP({seg_id})", seg_id=seg_id, records=records, residues=residues, chains=chains,
                     basepairs=basepairs)


def _key(struct: Structure, chain_res: Tuple[str, int]) -> ResidKey:
    return (struct.seg_id, chain_res[0], int(chain_res[1]))


def _pair_keys(struct: Structure, pair) -> Tuple[ResidKey, ResidKey]:
    return _key(struct, pair[0]), _key(struct, pair[1])


def helix_phosphate_radius(struct: Structure) -> float:
    """Median distance of P atoms from the C1' PCA axis line."""
    point, direction = struct.axis_line()
    radii = []
    for res in struct.residues.values():
        rec = res.atom_by_name.get("P")
        if rec is None:
            continue
        v = np.array([rec.x, rec.y, rec.z], dtype=float) - point
        radii.append(float(np.linalg.norm(v - np.dot(v, direction) * direction)))
    if not radii:
        raise ValueError(f"{struct.name}: no P atoms; cannot tell B-DNA from P-DNA")
    return float(np.median(radii))


def detect_helix_type(struct: Structure) -> str:
    return "P" if helix_phosphate_radius(struct) < P_HELIX_MAX_P_RADIUS else "B"


def parse_types(text: str, n_files: int) -> List[str]:
    letters = [ch for ch in (text or "").upper() if not ch.isspace() and ch != ","]
    if len(letters) != n_files or any(ch not in "BP" for ch in letters):
        raise ValueError(f"--types needs one letter B or P per input file ({n_files}), e.g. BPB; got {text!r}")
    return letters


# ----------------------------------------------------------------------------
# Axis lines
# ----------------------------------------------------------------------------

def screw_axis_line(struct: Structure, extra_exclude: Optional[set] = None) -> AxisLineResult:
    """Axis of a regular helix from the screw that maps each level onto the next.

    Always used for P-DNA, where DSSR finds no helix, and for B-DNA when DSSR is
    unavailable or the screw source is chosen. Backbone and sugar atoms shared by
    consecutive levels (both strands) enter one least-squares fit.
    """
    exclude = set(extra_exclude or set()) | struct.trim | struct.pre_trim
    levels = [bp for bp in struct.basepairs if bp[0] not in exclude and bp[1] not in exclude]
    if len(levels) < 3:
        raise ValueError(f"{struct.name}: need at least 3 levels for a screw axis, found {len(levels)}")
    src: List[np.ndarray] = []
    dst: List[np.ndarray] = []
    for (a0, b0), (a1, b1) in zip(levels, levels[1:]):
        for r0, r1 in ((a0, a1), (b0, b1)):
            res0, res1 = struct.residues[r0], struct.residues[r1]
            for name in ("P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"):
                x0, x1 = res0.xyz(name), res1.xyz(name)
                if x0 is not None and x1 is not None:
                    src.append(x0)
                    dst.append(x1)
    rot, trans, rmsd = kabsch(np.array(src), np.array(dst))
    w, v = np.linalg.eig(rot)
    direction = np.real(v[:, int(np.argmin(np.abs(w - 1.0)))])
    direction = direction / np.linalg.norm(direction)
    first = _bp_centroid(struct, levels[0])
    last = _bp_centroid(struct, levels[-1])
    if np.dot(last - first, direction) < 0:
        direction = -direction
    perp = trans - np.dot(trans, direction) * direction
    point, *_rest = np.linalg.lstsq(np.eye(3) - rot, perp, rcond=None)
    point = point + np.dot(first - point, direction) * direction
    twist = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(rot) - 1.0) / 2.0))))
    return AxisLineResult(
        point=point,
        direction=direction,
        source="screw",
        detail=f"rise={np.dot(trans, direction):.3f} A twist={twist:.2f} deg fit_rmsd={rmsd:.4f} A",
    )


AXIS_SOURCES = ("auto", "dssr", "screw", "pca")


def helix_axis_line(struct: Structure, kind: str, source: str, verbose_prefix: str, label: str,
                    extra_exclude: Optional[set] = None) -> AxisLineResult:
    """Axis of one helix. P-DNA always uses its screw axis (DSSR sees no helix).

    For B-DNA: auto tries DSSR and falls back to the screw axis; dssr requires
    DSSR; screw fits the base-pair-to-base-pair screw; pca is the C1' PCA line,
    which is biased for helices shorter than about two turns.
    """
    src = (source or "auto").strip().lower()
    if kind == "P" or src == "screw":
        result = screw_axis_line(struct, extra_exclude=extra_exclude)
    elif src in ("auto", "dssr"):
        try:
            return get_axis_line(struct, source="dssr", extra_exclude=extra_exclude, verbose_prefix=verbose_prefix,
                                 label=label)
        except Exception as exc:
            if src == "dssr":
                raise
            print(f"{verbose_prefix} axis[{label}] DSSR unavailable ({exc}); using the screw axis instead")
            result = screw_axis_line(struct, extra_exclude=extra_exclude)
    elif src == "pca":
        return get_axis_line(struct, source="pca", extra_exclude=extra_exclude, verbose_prefix=verbose_prefix,
                             label=label)
    else:
        raise ValueError(f"Unknown axis source: {source}")
    print(f"{verbose_prefix} axis[{label}] source={result.source}; {result.detail}")
    return result


def enforce_axis_alignment(
    *,
    fixed_helix: Structure,
    fixed_kind: str,
    moving_helix: Structure,
    moving_kind: str,
    pivot: np.ndarray,
    junction_dir: np.ndarray,
    verbose_prefix: str,
    mode: str,
    axis_source: str,
) -> None:
    """Make the moving helix axis co-directional (and optionally collinear) with the fixed one."""
    mode = normalize_axis_mode(mode)
    if mode == "none":
        print(f"{verbose_prefix} axis-mode=none: no post-fit axis correction applied.")
        return

    def axis(struct: Structure, kind: str, label: str) -> Tuple[np.ndarray, np.ndarray]:
        line = helix_axis_line(struct, kind, axis_source, verbose_prefix, label)
        direction = line.direction if np.dot(line.direction, junction_dir) >= 0 else -line.direction
        return line.point, direction

    pf, u = axis(fixed_helix, fixed_kind, "fixed")
    pm, v = axis(moving_helix, moving_kind, "moving-before")
    angle_before = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u, v))))))
    offset_before = axis_line_distance(pf, u, pm, v)
    if angle_before >= 1e-3:
        apply_rotation_about_point(moving_helix.records, rotation_matrix_from_vectors(v, u), pivot)
    pm2, v2 = axis(moving_helix, moving_kind, "moving-after-rotation")
    angle_after = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u, v2))))))
    offset_after = axis_line_distance(pf, u, pm2, v2)
    if mode == "codirectional":
        print(
            f"{verbose_prefix} codirectional: axis angle {angle_before:.3f} -> {angle_after:.3f} deg; "
            f"axis-line offset {offset_before:.3f} -> {offset_after:.3f} A (no axis-line shift)."
        )
        return
    shift = perpendicular_shift_to_collinear(pf, u, pm2)
    shift_mag = float(np.linalg.norm(shift))
    if shift_mag > 1e-8:
        apply_translation(moving_helix.records, shift)
    pm3, v3 = axis(moving_helix, moving_kind, "moving-after-shift")
    print(
        f"{verbose_prefix} collinear: axis angle {angle_before:.3f} -> {angle_after:.3f} deg; axis-line offset "
        f"{offset_before:.3f} -> {offset_after:.3f} -> {axis_line_distance(pf, u, pm3, v3):.3f} A; "
        f"perpendicular shift={shift_mag:.3f} A."
    )


# ----------------------------------------------------------------------------
# Core fitting
# ----------------------------------------------------------------------------

def _coords_for_pairs(struct: Structure, pairs: Sequence[Tuple[ResidKey, ResidKey]], swap: bool) -> np.ndarray:
    pts: List[np.ndarray] = []
    for r1, r2 in pairs:
        if swap:
            r1, r2 = r2, r1
        for rk in (r1, r2):
            res = struct.residues[rk]
            for name in SUGAR_ATOMS:
                xyz = res.xyz(name)
                if xyz is None:
                    raise ValueError(f"{struct.name}: residue {res.chainID}{res.resSeq} lacks atom {name}")
                pts.append(xyz)
    return np.stack(pts, axis=0)


def core_fit_pairs(core: Structure, kind: str) -> List[Tuple[ResidKey, ResidKey]]:
    """Core pairs used for fitting at one end, interior -> outer."""
    if kind == "B":
        return [_pair_keys(core, pair) for pair in core_BP.B_END_BASEPAIRS]
    if kind == "P":
        return [_pair_keys(core, pair) for pair in core_BP.P_END_LEVELS]
    raise ValueError(f"core end must be 'B' or 'P', got {kind!r}")


def helix_fit_pairs(helix: Structure, end: str, kind: str) -> List[Tuple[ResidKey, ResidKey]]:
    """Helix pairs matched to core_fit_pairs, in the same order.

    B: the helix's terminal pair sits on the core's interior B pair and the next
    one on the core's outer B pair. P: the first N_TRANSITION levels become the
    transition, so the core's interior and outer P levels sit on the next two.
    """
    if kind == "B":
        return list(reversed(helix.end_basepairs(end, n_bp=B_END_TRIM)))
    n = N_TRANSITION + len(core_BP.P_END_LEVELS)
    outer_first = list(reversed(helix.end_basepairs(end, n_bp=n)))  # terminal level first
    return outer_first[N_TRANSITION:]


def helix_replaced_pairs(helix: Structure, end: str, kind: str) -> List[Tuple[ResidKey, ResidKey]]:
    """Helix pairs the core replaces, terminal first, matched to core_replacing_pairs."""
    n = B_END_TRIM if kind == "B" else P_END_TRIM
    return list(reversed(helix.end_basepairs(end, n_bp=n)))


def core_replacing_pairs(core: Structure, kind: str) -> List[Tuple[ResidKey, ResidKey]]:
    """Core pairs that take the place of helix_replaced_pairs (same order)."""
    if kind == "B":
        return [_pair_keys(core, pair) for pair in core_BP.B_END_BASEPAIRS]
    return [_pair_keys(core, pair) for pair in core_BP.TRANSITION_LEVELS]


def helix_anchor_pair(helix: Structure, end: str, kind: str) -> Tuple[ResidKey, ResidKey]:
    """First helix pair kept next to the core (bonded to core_link_pair)."""
    n = (B_END_TRIM if kind == "B" else P_END_TRIM) + 1
    return helix.end_basepairs(end, n_bp=n)[0]


def core_link_pair(core: Structure, kind: str) -> Tuple[ResidKey, ResidKey]:
    """Core pair bonded to the helix anchor pair at that end."""
    if kind == "B":
        return _pair_keys(core, core_BP.B_END_BASEPAIRS[-1])
    return _pair_keys(core, core_BP.TRANSITION_LEVELS[-1])


def core_through_vector(core: Structure, kind_from: str, kind_to: str) -> np.ndarray:
    p_from = _bp_centroid(core, core_fit_pairs(core, kind_from)[-1])
    p_to = _bp_centroid(core, core_fit_pairs(core, kind_to)[-1])
    v = p_to - p_from
    return v / np.linalg.norm(v)


def end_outward_vector(struct: Structure, end: str) -> np.ndarray:
    bps = struct.end_basepairs(end, n_bp=2)
    v = _bp_centroid(struct, bps[1]) - _bp_centroid(struct, bps[0])
    return v / np.linalg.norm(v)


def align_core_to_helix(*, helix: Structure, end: str, kind: str, core: Structure,
                        verbose_prefix: str) -> Tuple[float, bool]:
    helix_pairs = helix_fit_pairs(helix, end, kind)
    core_pairs = core_fit_pairs(core, kind)
    best = None
    P = _coords_for_pairs(core, core_pairs, swap=False)
    for swap in (False, True):
        Q = _coords_for_pairs(helix, helix_pairs, swap=swap)
        R, t, rmsd = kabsch(P, Q)
        if best is None or rmsd < best[0]:
            best = (rmsd, swap, R, t)
    rmsd, swap, R, t = best
    apply_rigid_transform(core.records, R, t)
    print(
        f"{verbose_prefix} align core({kind}-end) -> fixed[{helix.name}] {end}-end: rmsd={rmsd:.3f} A, "
        f"strand_swap={'Y' if swap else 'N'}  [{len(core_pairs)} {'bp' if kind == 'B' else 'levels'}]"
    )
    return rmsd, swap


def _pick(candidates: Sequence[tuple], preferred_end: str) -> tuple:
    """Lowest RMSD, but prefer ``preferred_end`` among fits that are as good.

    P-DNA has a dyad at every level, so both of its ends fit a core equally
    well; preferring one end keeps the result deterministic and keeps each
    input's chain A in one output strand for bnp_na-built helices.
    """
    best = min(c[0] for c in candidates)
    return min(candidates, key=lambda c: (c[0] > best + END_TIE_RMSD, c[1] != preferred_end, c[0]))


def align_helix_to_core(*, helix: Structure, kind: str, core: Structure, junction_dir: np.ndarray,
                        verbose_prefix: str) -> Tuple[str, float, bool, np.ndarray, float]:
    """Move a helix onto the core end of its kind; pick the end/strand mapping that extends forward."""
    core_pairs = core_fit_pairs(core, kind)
    Q = _coords_for_pairs(core, core_pairs, swap=False)
    pivot = Q.mean(axis=0)
    good = []
    every = []
    for end in ("start", "end"):
        helix_pairs = helix_fit_pairs(helix, end, kind)
        connect = helix.basepairs[0] if end == "start" else helix.basepairs[-1]
        far = helix.basepairs[-1] if end == "start" else helix.basepairs[0]
        p_connect0 = _bp_centroid(helix, connect)
        p_far0 = _bp_centroid(helix, far)
        for swap in (False, True):
            P = _coords_for_pairs(helix, helix_pairs, swap=swap)
            R, t, rmsd = kabsch(P, Q)
            ext = (R @ p_far0 + t) - (R @ p_connect0 + t)
            ext_len = float(np.linalg.norm(ext))
            dot_unit = float(np.dot(ext, junction_dir)) / ext_len if ext_len > 1e-12 else -1.0
            cand = (rmsd, end, swap, R, t, dot_unit)
            every.append(cand)
            if dot_unit > 0:
                good.append(cand)
    chosen = _pick(good if good else every, preferred_end="start")
    best_good = chosen if good else None
    rmsd, end, swap, R, t, dot_unit = chosen
    apply_rigid_transform(helix.records, R, t)
    msg = (
        f"{verbose_prefix} align helix[{helix.name}] {end}-end -> core({kind}-end): rmsd={rmsd:.3f} A, "
        f"strand_swap={'Y' if swap else 'N'}, dir_cos={dot_unit:.3f}"
    )
    if best_good is None:
        msg += "  WARNING: no candidate extends forward; using the lowest RMSD"
    print(msg)
    return end, rmsd, swap, pivot, dot_unit


# ----------------------------------------------------------------------------
# Base grafting onto core residues
# ----------------------------------------------------------------------------

def _atom_line(serial: int, name: str, res_name: str, chain: str, res_seq: int, xyz: Sequence[float]) -> str:
    field_name = name if len(name) >= 4 else f" {name:<3s}"
    element = name.lstrip("0123456789")[:1]
    return (
        f"ATOM  {serial:5d} {field_name} {res_name:>3s} {chain}{res_seq:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def regraft_residue(struct: Structure, key: ResidKey, new_base: str) -> str:
    """Give one residue a new base in the standard frame; returns the old base letter."""
    res = struct.residues[key]
    old = base_code(res.resName)
    atoms = {name: np.array([rec.x, rec.y, rec.z], dtype=float) for name, rec in res.atom_by_name.items()}
    new_atoms, _rmsd = graft_base(atoms, old, new_base)
    res_name = "D" + new_base
    new_records = [
        pdb_atom_record(_atom_line(0, name, res_name, res.chainID, res.resSeq, xyz)) for name, xyz in new_atoms
    ]
    old_ids = {id(rec) for rec in res.atoms}
    position = next(i for i, rec in enumerate(struct.records) if id(rec) in old_ids)
    struct.records = [rec for rec in struct.records if id(rec) not in old_ids]
    struct.records[position:position] = new_records
    res.atoms = new_records
    res.atom_by_name = {rec.name: rec for rec in new_records}
    res.resName = res_name
    return old


def graft_core_bases(*, helix: Structure, end: str, kind: str, core: Structure, swap: bool,
                     verbose_prefix: str) -> None:
    helix_pairs = helix_replaced_pairs(helix, end, kind)
    core_pairs = core_replacing_pairs(core, kind)
    done = []
    for (h1, h2), (c_e, c_f) in zip(helix_pairs, core_pairs):
        targets = ((h1, c_f), (h2, c_e)) if swap else ((h1, c_e), (h2, c_f))
        for h_key, c_key in targets:
            new_base = base_code(helix.residues[h_key].resName)
            old_base = base_code(core.residues[c_key].resName)
            regraft_residue(core, c_key, new_base)
            done.append(f"{c_key[1]}{c_key[2]}:{old_base}->{new_base}")
    print(f"{verbose_prefix} core bases from helix[{helix.name}] {end}-end ({kind}): {', '.join(done)}")


# ----------------------------------------------------------------------------
# Junction bookkeeping and ligation
# ----------------------------------------------------------------------------

@dataclass
class BPJunction:
    jid: int
    core: Structure
    sides: List[Tuple[Structure, str, str, bool]] = field(default_factory=list)  # (helix, end, kind, swap)


def build_backbone_graph(structs: Sequence[Structure], junctions: Sequence[BPJunction]) -> Dict[ResidKey, set]:
    adj: Dict[ResidKey, set] = {}

    def add(a: ResidKey, b: ResidKey) -> None:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    for st in structs:
        for chain in st.chains:
            keys = sorted(
                (k for k, r in st.residues.items() if r.chainID == chain and k not in st.trim and k not in st.pre_trim),
                key=lambda k: st.residues[k].resSeq,
            )
            for k1, k2 in zip(keys, keys[1:]):
                add(k1, k2)
    for jx in junctions:
        for helix, end, kind, swap in jx.sides:
            h1, h2 = helix_anchor_pair(helix, end, kind)
            c_e, c_f = core_link_pair(jx.core, kind)
            if swap:
                add(h1, c_f)
                add(h2, c_e)
            else:
                add(h1, c_e)
                add(h2, c_f)
    return adj


def junction_bond_report(jx: BPJunction, res_by_key: Dict[ResidKey, Residue]) -> List[str]:
    lines = []
    for helix, end, kind, swap in jx.sides:
        h1, h2 = helix_anchor_pair(helix, end, kind)
        c_e, c_f = core_link_pair(jx.core, kind)
        pairs = ((h1, c_f), (h2, c_e)) if swap else ((h1, c_e), (h2, c_f))
        for a, b in pairs:
            d = _distance_min_O3P(res_by_key[a], res_by_key[b])
            d_text = f"{d:.3f} A" if d is not None else "NA"
            lines.append(f"   {kind}-side [{helix.name}] {a[1]}{a[2]} <-> core {b[1]}{b[2]}: O3'-P {d_text}")
    return lines


def write_ligated_pdb(*, structs: Sequence[Structure], junctions: Sequence[BPJunction], out_path: Path,
                      first_chain_keys: Optional[set] = None) -> None:
    """Write the two ligated strands as chains A and B, numbered 5'->3'.

    Chain A is the strand holding ``first_chain_keys`` (the first input's first
    chain) when given, so chain A of bnp_na-built inputs stays chain A.
    """
    res_by_key: Dict[ResidKey, Residue] = {}
    for st in structs:
        for k, r in st.residues.items():
            if k not in st.trim and k not in st.pre_trim:
                res_by_key[k] = r
    nodes = sorted(res_by_key, key=lambda k: (k[0], k[1], k[2]))
    print(f"Total residues after trimming (all segments + cores): {len(nodes)}")
    for jx in junctions:
        print(f"[J{jx.jid}] junction bonds:")
        for line in junction_bond_report(jx, res_by_key):
            print(line)
    adj = build_backbone_graph(structs, junctions)
    comps = connected_components_keys(adj, nodes)
    comps.sort(key=len, reverse=True)
    print(f"Ligation graph: components={len(comps)} sizes={[len(c) for c in comps]}")
    if len(comps) != 2:
        raise ValueError(f"Ligation expected exactly 2 strands, found {len(comps)} components")
    if first_chain_keys and not (set(comps[0]) & first_chain_keys) and (set(comps[1]) & first_chain_keys):
        comps = [comps[1], comps[0]]
    orders = []
    for label, comp in (("[strand A]", comps[0]), ("[strand B]", comps[1])):
        order = order_path_component_keys(adj, comp)
        orders.append(ensure_order_5to3(order=order, res_by_key=res_by_key, verbose_label=label))
    out_recs: List[object] = []
    for chain_id, order in zip(("A", "B"), orders):
        for new_seq, key in enumerate(order, start=1):
            res = res_by_key[key]
            for rec in res.atoms:
                rec.update_chainID(chain_id)
                rec.update_resSeq(new_seq)
            out_recs.extend(res.atoms)
        last = res_by_key[order[-1]]
        out_recs.append(make_ter_record(chain_id, last.resName, len(order)))
        seq = "".join(base_code(res_by_key[k].resName) for k in order)
        print(f"  Strand {chain_id} ({len(order)} nt, 5'->3'): {seq}")
    with out_path.open("w") as fout:
        rec2file(out_recs, fout, reorder_serial=True)
        fout.write("END\n")
    print(f"Wrote output PDB: {out_path}")


# ----------------------------------------------------------------------------
# Pipeline
# ----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Build multi-segment B-DNA/P-DNA constructs joined by B-P junction cores. Inputs alternate "
            "between B-DNA and P-DNA PDB files, in either order (B P B ... or P B P ...)."
        ),
        epilog=(
            "Lengths: at each junction the P-DNA's 2 terminal levels become the transition (bases kept), so a "
            "P-DNA of N levels keeps N - 2 P-form levels with one junction and N - 4 between two junctions; for "
            "k P-form levels, enter k + 2 (or k + 4). A P-DNA needs at least 4 levels (7 with junctions at both "
            "ends). A B-DNA keeps every base pair and needs at least 3 bp (5 with junctions at both ends)."
        ),
    )
    p.add_argument("dna_files", nargs="+", help="Input PDB files in construct order.")
    p.add_argument("--out", default=DEFAULT_OUTPUT_NAME,
                   help="Output PDB. Default: make_BP_out.pdb in the first input file's folder.")
    p.add_argument("--types", default="",
                   help="Helix type per input, e.g. BPB. Default: detect from the phosphate radius.")
    p.add_argument("--axis-mode", choices=["codirectional", "collinear", "none"], default="codirectional",
                   help="Axis correction across each junction (default codirectional).")
    p.add_argument("--axis-source", choices=list(AXIS_SOURCES), default="auto",
                   help=("B-DNA axis source: auto = DSSR, falling back to the screw axis; dssr; screw = the "
                         "base-pair-to-base-pair screw; pca = C1' PCA (biased for short helices). "
                         "P-DNA always uses its screw axis."))
    return p


def resolve_output_path(out_text: str, input_paths: Sequence[Path]) -> Path:
    out_path = Path(out_text)
    if out_text == DEFAULT_OUTPUT_NAME and not out_path.is_absolute() and input_paths:
        return input_paths[0].expanduser().resolve().parent / DEFAULT_OUTPUT_NAME
    return out_path


def run_pipeline(args: argparse.Namespace) -> Dict[str, object]:
    axis_mode = normalize_axis_mode(getattr(args, "axis_mode", "codirectional"))
    axis_source = (getattr(args, "axis_source", "auto") or "auto").strip().lower()
    if axis_source not in AXIS_SOURCES:
        raise SystemExit(f"ERROR: unknown axis source {axis_source!r}")
    paths = [Path(p) for p in args.dna_files]
    for p in paths:
        if not p.exists():
            raise SystemExit(f"ERROR: input file not found: {p}")
    if len(paths) < 2:
        raise SystemExit("ERROR: need at least 2 files (one B-DNA and one P-DNA)")

    helices = [load_helix(p, seg_id=f"seg{i + 1}") for i, p in enumerate(paths)]
    if getattr(args, "types", ""):
        try:
            types = parse_types(args.types, len(paths))
        except ValueError as exc:
            raise SystemExit(f"ERROR: {exc}") from exc
        source = "given by --types"
    else:
        types = [detect_helix_type(h) for h in helices]
        source = "detected from the phosphate radius"

    print(f"Input files ({source}):")
    for i, (t, p, h) in enumerate(zip(types, paths, helices), start=1):
        print(f"  {i:2d}. {t}  {p}  (P radius {helix_phosphate_radius(h):.2f} A, "
              f"{len(h.basepairs)} {'bp' if t == 'B' else 'levels'})")
    for i, (t1, t2) in enumerate(zip(types, types[1:]), start=1):
        if t1 == t2:
            raise SystemExit(
                f"ERROR: inputs must alternate between B-DNA and P-DNA; files {i} and {i + 1} are both {t1}-type"
            )
    for i, (t, h) in enumerate(zip(types, helices)):
        ends = (1 if i > 0 else 0) + (1 if i < len(helices) - 1 else 0)
        need = (MIN_B_BP if t == "B" else MIN_P_LEVELS)[ends]
        if len(h.basepairs) < need:
            unit = "bp" if t == "B" else "levels"
            raise SystemExit(
                f"ERROR: {h.name} has {len(h.basepairs)} {unit}; at least {need} are needed for "
                f"{ends} junction end(s)"
            )
    print(f"Junctions to build: {len(paths) - 1}")
    print(f"Axis correction mode: {axis_mode}; B-DNA axis source: {axis_source}; P-DNA axis: screw fit")
    print(
        f"Each junction: the B-DNA end gives its {B_END_TRIM} terminal bp to the core (bases kept); "
        f"the P-DNA end's {P_END_TRIM} terminal levels become the transition (bases kept)."
    )

    helices[0].backward_end = "start"
    all_structs: List[Structure] = [helices[0]]
    junctions: List[BPJunction] = []
    for j in range(len(paths) - 1):
        fixed, moving = helices[j], helices[j + 1]
        t_fixed, t_next = types[j], types[j + 1]
        prefix = f"[J{j + 1}]"
        fixed_end = fixed.forward_end()
        if j == 0:
            good = []
            every = []
            for try_end in (fixed_end, "start" if fixed_end == "end" else "end"):
                core_try = new_core(seg_id=f"core{j + 1}")
                rmsd_try, swap_try = align_core_to_helix(
                    helix=fixed, end=try_end, kind=t_fixed, core=core_try, verbose_prefix=prefix + " [try]"
                )
                dot_dir = float(np.dot(core_through_vector(core_try, t_fixed, t_next), end_outward_vector(fixed, try_end)))
                print(f"{prefix} [try] direction check: dot(core_dir, fixed_out)={dot_dir:.3f} (try_end={try_end})")
                cand = (rmsd_try, try_end, core_try, swap_try)
                every.append(cand)
                if dot_dir > 0:
                    good.append(cand)
            _rmsd, fixed_end, core, swap_fixed = _pick(good if good else every, preferred_end="end")
            if not good:
                print(f"{prefix} WARNING: could not find an outward-pointing core placement; using lowest RMSD")
            fixed.backward_end = "start" if fixed_end == "end" else "end"
            print(f"{prefix} selected fixed_end={fixed_end} for the first junction")
        else:
            core = new_core(seg_id=f"core{j + 1}")
            _rmsd, swap_fixed = align_core_to_helix(
                helix=fixed, end=fixed_end, kind=t_fixed, core=core, verbose_prefix=prefix
            )
            dot_dir = float(np.dot(core_through_vector(core, t_fixed, t_next), end_outward_vector(fixed, fixed_end)))
            if dot_dir <= 0:
                print(f"{prefix} WARNING: core does not point outward from the fixed helix end (dot={dot_dir:.3f})")

        junction_dir = core_through_vector(core, t_fixed, t_next)
        chosen_end, _rmsd2, swap_next, pivot, dot_ext = align_helix_to_core(
            helix=moving, kind=t_next, core=core, junction_dir=junction_dir, verbose_prefix=prefix
        )
        moving.backward_end = chosen_end
        if dot_ext <= 0:
            print(f"{prefix} WARNING: next helix does not extend forward along the junction (dir_cos={dot_ext:.3f})")
        if axis_mode != "none":
            enforce_axis_alignment(
                fixed_helix=fixed, fixed_kind=t_fixed, moving_helix=moving, moving_kind=t_next, pivot=pivot,
                junction_dir=junction_dir, verbose_prefix=prefix, mode=axis_mode, axis_source=axis_source,
            )

        graft_core_bases(helix=fixed, end=fixed_end, kind=t_fixed, core=core, swap=swap_fixed, verbose_prefix=prefix)
        graft_core_bases(helix=moving, end=chosen_end, kind=t_next, core=core, swap=swap_next, verbose_prefix=prefix)

        for helix, end, kind in ((fixed, fixed_end, t_fixed), (moving, chosen_end, t_next)):
            pairs = helix_replaced_pairs(helix, end, kind)
            for bp in pairs:
                helix.trim.update(bp)
            unit = "bp" if kind == "B" else "levels"
            print(f"{prefix} trim: [{helix.name}] {end}-end {len(pairs)} {unit} (taken over by the core)")
        for bp in core_fit_pairs(core, "P"):
            core.trim.update(bp)

        junctions.append(BPJunction(
            jid=j + 1, core=core, sides=[(fixed, fixed_end, t_fixed, swap_fixed), (moving, chosen_end, t_next, swap_next)]
        ))
        all_structs.append(core)
        all_structs.append(moving)

    out_path = resolve_output_path(args.out, paths)
    if str(out_path) != str(Path(args.out)):
        print(f"Output path defaulted to the first input folder: {out_path}")
    raw_path = raw_output_path(out_path)
    write_raw_aligned_pdb(structs=all_structs, out_path=raw_path)
    first = helices[0]
    first_chain_keys = {k for k, r in first.residues.items() if r.chainID == first.chains[0]}
    write_ligated_pdb(structs=all_structs, junctions=junctions, out_path=out_path, first_chain_keys=first_chain_keys)
    return {"pdb_out": out_path, "pdb_raw": raw_path, "types": types, "axis_mode": axis_mode,
            "axis_source": axis_source}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(list(sys.argv[1:] if argv is None else argv))
    run_pipeline(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
