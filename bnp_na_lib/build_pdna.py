#!/usr/bin/env python3
"""P-DNA builder for bnp_na.

P-DNA is the Pauling-like DNA of the Pf1 phage (PDB 1PFI): the phosphates of
both strands wind around the helix axis at r = 2.5 A and every base points
outward, so the two strands form no base pairs. The backbone comes from the
1PFI template level in core_PDNA.py, repeated by the template's exact screw
(131.84 deg and 6.10 A per nucleotide), and each base is replaced by the
requested one in the standard base reference frame. The helical parameters
are therefore fixed and cannot be customized.

Because no base pairs form, the two strands take independent sequences. Strand
II defaults to the reverse complement of strand I, which puts complementary
bases on the same level. Both strands must have the same length: level k holds
strand I nucleotide k and strand II nucleotide N + 1 - k.

Pipeline: write the template-derived model with strand I (chain A) running
5'->3' along +Z and level 1 centered on the origin -> normalize PDB naming ->
optional phenix.geometry_minimization -> optional phosphate regularization ->
align to +Z by superposing onto the ideal model. DSSR is not used: it finds no
base pairs, stacks or helix in P-DNA, so DSSR --more cannot supply an axis.
"""
from __future__ import annotations

if __package__ in (None, ""):  # run as a script: python3 bnp_na_lib/build_pdna.py
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    __package__ = "bnp_na_lib"

import argparse
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import core_PDNA
from .build_common import (
    PipelineError,
    command_to_text,
    expand_sequence,
    expected_phenix_minimized_path,
    run_phenix_minimization,
    sanitize_basename,
    stage_params_to_output_dir,
)
from .pdb_name_standard import normalize_pdb_naming as normalize_nucleotide_pdb_naming
from .regularize_phosphates import default_regularized_output_path, regularize_phosphates
from .standard_bases import (
    DNA_RESIDUE_NAMES,
    PURINES,
    chi_deg,
    dna_reverse_complement,
    graft_base,
    is_base_atom,
    kabsch,
    rotate_base_about_glycosidic_bond,
)


BACKBONE = "P-DNA"
TEMPLATE_TAG = "1PFI"
DNA_ALPHABET = "ATCG"
# The template is poly(dC). With chi unchanged, the minor-groove edges of two
# purines (strand I level k, strand II level k+1) come within 2.88 A (C2...N2).
# Turning every purine by -8 deg about its glycosidic bond keeps all base
# contacts at or above the template's own 3.64 A (O2...O2 of two cytosines).
PURINE_CHI_OFFSET_DEG = -8.0

Atom = Tuple[str, np.ndarray]
Residue = Tuple[str, List[Atom]]  # (residue name, atoms)


def resolve_pdna_sequences(strand1: str, strand2: Optional[str] = None) -> Tuple[str, str, bool]:
    """Expand both strand sequences (5'->3') and fill in strand II when blank.

    Returns (strand I, strand II, strand II was defaulted). Strand II defaults
    to the reverse complement of strand I. Raises ValueError for invalid or
    unequal-length sequences.
    """
    seq1 = expand_sequence(strand1, alphabet=DNA_ALPHABET)
    if not seq1:
        raise ValueError("Strand I sequence is empty.")
    text2 = (strand2 or "").strip()
    if not text2:
        return seq1, dna_reverse_complement(seq1), True
    seq2 = expand_sequence(text2, alphabet=DNA_ALPHABET)
    if len(seq2) != len(seq1):
        raise ValueError(
            f"Strand II has {len(seq2)} nt but strand I has {len(seq1)} nt. P-DNA is built level by level, "
            "so both strands need the same length (leave strand II blank for the reverse complement)."
        )
    return seq1, seq2, False


def level_pairs(seq1: str, seq2: str) -> List[Tuple[int, str, int, str]]:
    """(level, strand I base, strand II residue number, strand II base) per level."""
    n = len(seq1)
    return [(k, seq1[k - 1], n + 1 - k, seq2[n - k]) for k in range(1, n + 1)]


def _template_units() -> Dict[str, List[Atom]]:
    units = core_PDNA.level_atoms()
    return {key: [(name, np.array(xyz, dtype=float)) for name, xyz in atoms] for key, atoms in units.items()}


def _screw(steps: int) -> Tuple[np.ndarray, np.ndarray]:
    rot, dz = core_PDNA.screw_matrix(steps)
    return np.array(rot, dtype=float), np.array([0.0, 0.0, dz], dtype=float)


def build_pdna_model(seq1: str, seq2: str) -> Dict[str, List[Residue]]:
    """Return {'A': strand I residues, 'B': strand II residues}, each 5'->3'.

    Level k (1-based) is the template level moved by k-1 screw steps, so level 1
    is centered on the origin and the helix axis is +Z.
    """
    if len(seq1) != len(seq2) or not seq1:
        raise ValueError("Both strands need the same, nonzero length.")
    units = _template_units()
    grafted: Dict[Tuple[str, str], List[Atom]] = {}
    for strand in ("I", "II"):
        atoms = dict(units[strand])
        for base in sorted(set(seq1 if strand == "I" else seq2)):
            new_atoms, _rmsd = graft_base(atoms, core_PDNA.TEMPLATE_BASE, base)
            if base in PURINES:
                new_atoms = rotate_base_about_glycosidic_bond(new_atoms, base, PURINE_CHI_OFFSET_DEG)
            grafted[(strand, base)] = new_atoms

    n = len(seq1)
    strand_a: List[Residue] = []
    for k, base in enumerate(seq1, start=1):
        rot, shift = _screw(k - 1)
        atoms = [(name, rot @ xyz + shift) for name, xyz in grafted[("I", base)]]
        strand_a.append((DNA_RESIDUE_NAMES[base], atoms))
    strand_b: List[Residue] = []
    for j, base in enumerate(seq2, start=1):
        rot, shift = _screw(n - j)
        atoms = [(name, rot @ xyz + shift) for name, xyz in grafted[("II", base)]]
        strand_b.append((DNA_RESIDUE_NAMES[base], atoms))
    return {"A": strand_a, "B": strand_b}


def _atom_name_field(name: str) -> str:
    return name if len(name) >= 4 else f" {name:<3s}"


def format_atom_line(serial: int, name: str, res_name: str, chain: str, res_seq: int, xyz: Sequence[float]) -> str:
    element = name.lstrip("0123456789")[:1]
    return (
        f"ATOM  {serial:5d} {_atom_name_field(name)} {res_name:>3s} {chain}{res_seq:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def model_pdb_lines(model: Dict[str, List[Residue]]) -> List[str]:
    lines: List[str] = []
    serial = 0
    for chain in ("A", "B"):
        residues = model[chain]
        for res_seq, (res_name, atoms) in enumerate(residues, start=1):
            for name, xyz in atoms:
                serial += 1
                lines.append(format_atom_line(serial, name, res_name, chain, res_seq, xyz))
        serial += 1
        last_name = residues[-1][0]
        lines.append(f"TER   {serial:5d}      {last_name:>3s} {chain}{len(residues):4d}\n")
    lines.append("END\n")
    return lines


def write_pdna_pdb(model: Dict[str, List[Residue]], out_path: Union[str, Path]) -> Path:
    path = Path(out_path)
    path.write_text("".join(model_pdb_lines(model)), encoding="utf-8")
    return path


def closest_nonbonded_contact(
    residues: Sequence[Tuple[str, int, Sequence[Atom]]],
    cutoff: float = 4.0,
) -> Tuple[float, Optional[Tuple[str, int, str, str, int, str]]]:
    """Closest heavy-atom contact between residues that are not covalent neighbors.

    ``residues`` holds (chain, index, atoms) with consecutive indices bonded in a
    chain. For bonded neighbors only contacts that involve a base atom count,
    which catches base-base and base-backbone clashes without the linkage.
    """
    entries = []
    for chain, index, atoms in residues:
        for name, xyz in atoms:
            entries.append((float(xyz[2]), chain, index, name, np.asarray(xyz, dtype=float)))
    entries.sort(key=lambda item: item[0])
    best = (float("inf"), None)
    for i, (zi, ci, ii, ni, xi) in enumerate(entries):
        for zj, cj, ij, nj, xj in entries[i + 1:]:
            if zj - zi >= cutoff:
                break
            if ci == cj and ii == ij:
                continue
            if ci == cj and abs(ii - ij) == 1 and not (is_base_atom(ni) or is_base_atom(nj)):
                continue
            d = float(np.linalg.norm(xi - xj))
            if d < best[0]:
                best = (d, (ci, ii, ni, cj, ij, nj))
    return best


def _model_checks(model: Dict[str, List[Residue]]) -> Dict[str, object]:
    """Linkage, glycosidic and contact statistics used in the build log and the tests."""
    linkage: List[float] = []
    for chain in ("A", "B"):
        residues = model[chain]
        for (_n1, atoms1), (_n2, atoms2) in zip(residues, residues[1:]):
            o3 = dict(atoms1)["O3'"]
            p = dict(atoms2)["P"]
            linkage.append(float(np.linalg.norm(o3 - p)))
    glycosidic: List[float] = []
    chis: Dict[str, List[float]] = {"purine": [], "pyrimidine": []}
    flat = []
    for chain in ("A", "B"):
        for index, (res_name, atoms) in enumerate(model[chain], start=1):
            names = dict(atoms)
            base = res_name[-1]
            n_atom = "N9" if base in PURINES else "N1"
            glycosidic.append(float(np.linalg.norm(names["C1'"] - names[n_atom])))
            chis["purine" if base in PURINES else "pyrimidine"].append(chi_deg(names, base))
            flat.append((chain, index, atoms))
    min_contact, contact_atoms = closest_nonbonded_contact(flat)
    return {
        "linkage_min": min(linkage) if linkage else float("nan"),
        "linkage_max": max(linkage) if linkage else float("nan"),
        "glycosidic_min": min(glycosidic),
        "glycosidic_max": max(glycosidic),
        "chi_purine": (min(chis["purine"]), max(chis["purine"])) if chis["purine"] else None,
        "chi_pyrimidine": (min(chis["pyrimidine"]), max(chis["pyrimidine"])) if chis["pyrimidine"] else None,
        "min_contact": min_contact,
        "min_contact_atoms": contact_atoms,
    }


def _format_checks(checks: Dict[str, object]) -> List[str]:
    lines = [
        f"O3'-P linkage: {checks['linkage_min']:.3f}-{checks['linkage_max']:.3f} Å (template geometry)",
        f"C1'-N glycosidic bond after base grafting: {checks['glycosidic_min']:.3f}-{checks['glycosidic_max']:.3f} Å",
    ]
    for label, key in (("pyrimidine", "chi_pyrimidine"), ("purine", "chi_purine")):
        value = checks.get(key)
        if value:
            lines.append(f"chi ({label}s): {value[0]:.1f} to {value[1]:.1f}°")
    if checks.get("chi_purine"):
        lines.append(
            f"Purines are turned {PURINE_CHI_OFFSET_DEG:+.0f}° about the glycosidic bond relative to the poly(dC) "
            "template, which keeps cross-strand purine edges apart."
        )
    atoms = checks.get("min_contact_atoms")
    where = ""
    if atoms:
        c1, i1, n1, c2, i2, n2 = atoms  # type: ignore[misc]
        where = f" ({c1}{i1} {n1} - {c2}{i2} {n2})"
    lines.append(f"Closest non-bonded heavy-atom contact: {checks['min_contact']:.3f} Å{where}")
    return lines


def _read_coordinates(path: Path) -> Dict[Tuple[str, int, str], np.ndarray]:
    coords: Dict[Tuple[str, int, str], np.ndarray] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            key = (line[21], int(line[22:26]), line[12:16].strip())
            coords[key] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
    return coords


def superpose_onto_reference(moving_pdb: Path, reference_pdb: Path, out_pdb: Path) -> Dict[str, float]:
    """Rigidly superpose ``moving_pdb`` onto ``reference_pdb`` over shared atoms.

    The reference is the ideal P-DNA model, whose axis is +Z by construction, so
    this restores the +Z frame after minimization or regularization.
    """
    moving = _read_coordinates(moving_pdb)
    reference = _read_coordinates(reference_pdb)
    keys = [key for key in reference if key in moving]
    if len(keys) < 3:
        raise ValueError(f"Too few shared atoms to superpose {moving_pdb.name} onto {reference_pdb.name}.")
    rot, trans, rmsd = kabsch(
        np.array([moving[key] for key in keys]),
        np.array([reference[key] for key in keys]),
    )
    out_lines = []
    for line in moving_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
            new = rot @ xyz + trans
            line = line[:30] + f"{new[0]:8.3f}{new[1]:8.3f}{new[2]:8.3f}" + line[54:]
        out_lines.append(line)
    out_pdb.write_text("".join(out_lines), encoding="utf-8")
    return {"rmsd": rmsd, "atoms": float(len(keys))}


def describe_parameters() -> str:
    return (
        f"h-Rise={core_PDNA.RISE:.4f} Å/nt, h-Twist={core_PDNA.TWIST_DEG:.4f}°/nt "
        f"({core_PDNA.NT_PER_TURN:.4f} nt/turn, pitch {core_PDNA.PITCH:.4f} Å), "
        f"P radius {core_PDNA.P_RADIUS:.2f} Å, C1' radius {core_PDNA.C1_RADIUS:.2f} Å"
    )


def build_pdna(
    seq_strand1: str,
    base_name: Optional[str],
    output_dir: Union[str, Path],
    seq_strand2: Optional[str] = None,
    deleteH: bool = False,
    run_phenix: bool = False,
    params_file: Optional[Union[str, Path]] = None,
    run_regularize_phosphates: Optional[bool] = None,
) -> Dict[str, object]:
    """Build P-DNA from the 1PFI template and return output paths plus log text."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logs: List[str] = []
    if run_regularize_phosphates is None:
        run_regularize_phosphates = bool(run_phenix)

    try:
        seq1, seq2, seq2_default = resolve_pdna_sequences(seq_strand1, seq_strand2)
    except Exception as exc:
        raise PipelineError(f"Invalid sequence: {exc}", str(exc)) from exc

    n = len(seq1)
    base = sanitize_basename(base_name or f"P-DNA{n}")
    if not base:
        raise PipelineError("Invalid helix name after sanitization.")

    pdb_raw = out_dir / f"{base}-{TEMPLATE_TAG}.pdb"
    try:
        model = build_pdna_model(seq1, seq2)
        write_pdna_pdb(model, pdb_raw)
        checks = _model_checks(model)
    except Exception as exc:
        raise PipelineError(f"Failed to build the P-DNA model: {exc}", "\n".join(logs)) from exc

    pairs = level_pairs(seq1, seq2)
    logs += [
        "=== P-DNA from the 1PFI template (core_PDNA.py) ===",
        f"Template: {core_PDNA.TEMPLATE_NAME}",
        f"Fixed helical parameters: {describe_parameters()}",
        f"Levels: {n} (one nucleotide per strand per level; no base pairs)",
        f"Strand I  (chain A, 5'->3' along +Z): {seq1}",
        f"Strand II (chain B, 5'->3' along -Z): {seq2}"
        + ("  [default: reverse complement of strand I]" if seq2_default else ""),
        "Level register (strand I residue / strand II residue):",
        "  " + ", ".join(f"{k}:A{k}{b1}/B{j}{b2}" for k, b1, j, b2 in pairs[:12])
        + (" ..." if len(pairs) > 12 else ""),
        *_format_checks(checks),
        f"Raw model: {pdb_raw}",
    ]

    pdb_norm_target = out_dir / f"{base}-{TEMPLATE_TAG}_out.pdb"
    try:
        pdb_norm = Path(normalize_nucleotide_pdb_naming(str(pdb_raw), str(pdb_norm_target), deleteH=deleteH))
        logs += [
            "\n=== Normalize PDB names (pdb_name_standard.py) ===",
            f"Input : {pdb_raw}",
            f"Output: {pdb_norm}",
            f"deleteH: {deleteH}",
        ]
    except Exception as exc:
        raise PipelineError(f"PDB naming normalization failed: {exc}", "\n".join(logs)) from exc

    align_source = pdb_norm
    pdb_min = None
    staged_params = None
    if run_phenix:
        if not params_file:
            raise PipelineError(
                "phenix.geometry_minimization was requested, but no params file was specified.", "\n".join(logs)
            )
        try:
            staged_params = stage_params_to_output_dir(params_file, out_dir)
        except Exception as exc:
            raise PipelineError(str(exc), "\n".join(logs)) from exc
        ok_phx, out_phx, cmd_phx = run_phenix_minimization(pdb_norm, staged_params)
        logs += [
            "\n=== phenix.geometry_minimization ===",
            f"Params : {staged_params}",
            f"Command: (cwd={out_dir}) {command_to_text(cmd_phx)}",
            f"Status : {'OK' if ok_phx else 'FAILED'}",
            f"Output :\n{out_phx}",
        ]
        if not ok_phx:
            raise PipelineError("phenix.geometry_minimization failed.", "\n".join(logs))
        pdb_min = expected_phenix_minimized_path(pdb_norm)
        if not pdb_min.exists():
            raise PipelineError(f"Expected minimized PDB was not found: {pdb_min}", "\n".join(logs))
        logs.append(f"Minimized PDB: {pdb_min}")
        align_source = pdb_min
    else:
        logs += ["\n=== phenix.geometry_minimization ===", "Skipped by user option."]

    pdb_regularized = None
    if run_regularize_phosphates:
        try:
            regularize_result = regularize_phosphates(align_source, default_regularized_output_path(align_source))
            pdb_regularized = regularize_result.output_pdb
            align_source = pdb_regularized
            logs += ["\n" + regularize_result.log_text]
        except Exception as exc:
            raise PipelineError(f"Phosphate regularization failed: {exc}", "\n".join(logs)) from exc
    else:
        logs += ["\n=== Regularize phosphates ===", "Skipped by user option."]

    pdb_aligned = align_source.with_name(align_source.stem + "_aligned2Z" + align_source.suffix)
    try:
        if align_source == pdb_norm:
            shutil.copyfile(pdb_norm, pdb_aligned)
            fit = {"rmsd": 0.0, "atoms": float(sum(len(atoms) for res in model.values() for _n, atoms in res))}
        else:
            fit = superpose_onto_reference(align_source, pdb_norm, pdb_aligned)
        logs += [
            "\n=== Align to +Z (superpose onto the ideal P-DNA model) ===",
            "DSSR finds no base pairs or helix in P-DNA, so the axis is taken from the template construction:",
            "helix axis = +Z, level 1 centered on the origin, strand I running 5'->3' along +Z.",
            f"Input : {align_source}",
            f"Reference: {pdb_norm}",
            f"Output: {pdb_aligned}",
            f"Superposition RMSD over {int(fit['atoms'])} atoms: {fit['rmsd']:.3f} Å",
        ]
    except Exception as exc:
        raise PipelineError(f"P-DNA align-to-Z failed: {exc}", "\n".join(logs)) from exc

    return {
        "na_type": BACKBONE,
        "base_name": base,
        "length": n,
        "strand1": seq1,
        "strand2": seq2,
        "strand2_default": seq2_default,
        "pdb_rebuild": pdb_raw,
        "pdb_normalized": pdb_norm,
        "pdb_minimized": pdb_min,
        "pdb_regularized": pdb_regularized,
        "pdb_aligned": pdb_aligned,
        "run_phenix": run_phenix,
        "run_regularize_phosphates": run_regularize_phosphates,
        "params_file": staged_params,
        "checks": checks,
        "log_text": "\n".join(logs),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build P-DNA (Pf1 phage DNA, PDB 1PFI template) aligned to +Z. Helical parameters are fixed: "
            + describe_parameters()
            + "."
        )
    )
    parser.add_argument("strand1", help="Strand I sequence 5'->3' (A/T/C/G; compact form such as A10T5 is accepted).")
    parser.add_argument(
        "--strand2",
        default="",
        help="Strand II sequence 5'->3', same length as strand I. Default: reverse complement of strand I.",
    )
    parser.add_argument("--name", default="", help="Helix name used for output files. Default: P-DNA<N>.")
    parser.add_argument("--out-dir", default=".", help="Output folder. Default: current folder.")
    parser.add_argument("--phenix", action="store_true", help="Run phenix.geometry_minimization.")
    parser.add_argument(
        "--params",
        default=str(Path(__file__).resolve().parent / "min_P_C5.params"),
        help="Params file for phenix.geometry_minimization. Default: bnp_na_lib/min_P_C5.params.",
    )
    parser.add_argument(
        "--regularize-phosphates",
        action="store_true",
        help="Regularize phosphates from C1' helical symmetry (needs at least 4 nt per strand).",
    )
    args = parser.parse_args(argv)
    try:
        result = build_pdna(
            args.strand1,
            args.name or None,
            args.out_dir,
            seq_strand2=args.strand2 or None,
            run_phenix=args.phenix,
            params_file=args.params if args.phenix else None,
            run_regularize_phosphates=args.regularize_phosphates,
        )
    except PipelineError as exc:
        print(getattr(exc, "log_text", "") or "")
        print(f"ERROR: {exc}")
        return 1
    print(result["log_text"])
    print(f"\nP-DNA aligned to +Z: {result['pdb_aligned']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
