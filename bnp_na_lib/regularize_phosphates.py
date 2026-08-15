#!/usr/bin/env python3
"""Regularize nucleic-acid phosphate groups with C1'-derived helical symmetry."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np


TOOL_VERSION = "V13.16"
PHOSPHATE_ATOMS = ("P", "OP1", "OP2")
BACKBONE_ATOMS = ("O5'", "C5'", "O3'")
REGULARIZED_ATOMS = PHOSPHATE_ATOMS + BACKBONE_ATOMS
PHOSPHATE_ALIASES = {"P": "P", "OP1": "OP1", "OP2": "OP2", "O1P": "OP1", "O2P": "OP2"}


class RegularizePhosphatesError(Exception):
    """Raised when phosphate regularization cannot be completed safely."""


@dataclass
class AtomRecord:
    line_index: int
    atom_name: str
    coord: np.ndarray


@dataclass
class Residue:
    model_id: int
    chain_id: str
    resseq: int
    icode: str
    order: int
    atoms: Dict[str, AtomRecord] = field(default_factory=dict)

    @property
    def label(self) -> str:
        chain = self.chain_id if self.chain_id.strip() else "<blank>"
        return f"chain {chain} residue {self.resseq}{self.icode.strip()}"


@dataclass(frozen=True)
class ChainRegularization:
    model_id: int
    chain_id: str
    nucleotide_count: int
    fitted_step_rmsd: float
    regularized_group_count: int
    regularized_atom_count: int
    terminal_group_count: int
    skipped_partial_group_count: int


@dataclass(frozen=True)
class RegularizePhosphatesResult:
    input_pdb: Path
    output_pdb: Path
    chains: Tuple[ChainRegularization, ...]
    log_text: str


def _normalized_atom_name(raw_name: str) -> str:
    name = raw_name.strip().upper().replace("*", "'")
    return PHOSPHATE_ALIASES.get(name, name)


def _parse_pdb(path: Path) -> Tuple[List[str], List[Residue]]:
    if not path.is_file():
        raise RegularizePhosphatesError(f"Input PDB not found: {path}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    residues: List[Residue] = []
    residue_lookup: Dict[Tuple[int, str, int, str], Residue] = {}
    model_id = 1
    order = 0

    for line_index, line in enumerate(lines):
        record = line[:6].strip().upper()
        if record == "MODEL":
            try:
                model_id = int(line[10:14].strip())
            except ValueError:
                model_id += 1
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            raise RegularizePhosphatesError(f"Could not parse ATOM/HETATM line {line_index + 1}: {line.rstrip()}")
        altloc = line[16:17]
        if altloc not in {" ", "A", "1"}:
            continue
        try:
            resseq = int(line[22:26])
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
        except ValueError as exc:
            raise RegularizePhosphatesError(
                f"Could not parse ATOM/HETATM line {line_index + 1}: {line.rstrip()}"
            ) from exc
        chain_id = line[21:22]
        icode = line[26:27]
        key = (model_id, chain_id, resseq, icode)
        residue = residue_lookup.get(key)
        if residue is None:
            residue = Residue(model_id, chain_id, resseq, icode, order)
            residue_lookup[key] = residue
            residues.append(residue)
            order += 1
        atom_name = _normalized_atom_name(line[12:16])
        # Prefer the blank-altloc record when both blank and A/1 are present.
        if atom_name not in residue.atoms or altloc == " ":
            residue.atoms[atom_name] = AtomRecord(line_index, atom_name, coord)

    if not residues:
        raise RegularizePhosphatesError(f"No ATOM/HETATM records found in {path}")
    return lines, residues


def _rigid_fit(source: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3 or len(source) < 3:
        raise RegularizePhosphatesError("At least three paired C1' points are required to fit helical symmetry.")
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u_matrix, _singular_values, vt_matrix = np.linalg.svd(covariance)
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0.0:
        vt_matrix[-1, :] *= -1.0
        rotation = vt_matrix.T @ u_matrix.T
    translation = target_center - rotation @ source_center
    predicted = (rotation @ source.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((predicted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def _homogeneous_transform(rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    transform = np.eye(4, dtype=float)
    transform[:3, :3] = rotation
    transform[:3, 3] = translation
    return transform


def _apply_transform(transform: np.ndarray, coord: np.ndarray) -> np.ndarray:
    return transform[:3, :3] @ coord + transform[:3, 3]


def _format_xyz(line: str, coord: np.ndarray) -> str:
    if not np.all(np.isfinite(coord)):
        raise RegularizePhosphatesError("Regularization produced a non-finite coordinate.")
    fields = tuple(f"{float(value):8.3f}" for value in coord)
    if any(len(field) > 8 for field in fields):
        raise RegularizePhosphatesError("A regularized coordinate cannot fit in the PDB 8.3 field.")
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    body = body.ljust(54)
    return body[:30] + "".join(fields) + body[54:] + newline


def _terminal_index_for_phosphate_only(
    residue: Residue,
    nucleotides: Sequence[Residue],
) -> Optional[int]:
    if residue.order < nucleotides[0].order:
        return -1
    if residue.order > nucleotides[-1].order:
        return len(nucleotides)
    return None


def _regularize_chain(
    chain_residues: Sequence[Residue],
    output_coords: Dict[int, np.ndarray],
) -> ChainRegularization:
    nucleotides = [residue for residue in chain_residues if "C1'" in residue.atoms]
    model_id = chain_residues[0].model_id
    chain_id = chain_residues[0].chain_id
    chain_label = chain_id if chain_id.strip() else "<blank>"
    if len(nucleotides) < 4:
        raise RegularizePhosphatesError(
            f"Model {model_id}, chain {chain_label} has {len(nucleotides)} C1' atoms; at least 4 are required."
        )

    c1_coords = np.array([residue.atoms["C1'"].coord for residue in nucleotides], dtype=float)
    rotation, translation, fit_rmsd = _rigid_fit(c1_coords[:-1], c1_coords[1:])
    step = _homogeneous_transform(rotation, translation)

    nucleotide_indices = {id(residue): index for index, residue in enumerate(nucleotides)}
    indexed_groups: List[Tuple[int, Residue, bool]] = []
    skipped_partial = 0
    for residue in chain_residues:
        present = [name for name in PHOSPHATE_ATOMS if name in residue.atoms]
        if not present:
            continue
        if len(present) != len(PHOSPHATE_ATOMS):
            skipped_partial += 1
            continue
        if id(residue) in nucleotide_indices:
            phosphate_index = nucleotide_indices[id(residue)]
            is_terminal = phosphate_index == 0
        else:
            terminal_index = _terminal_index_for_phosphate_only(residue, nucleotides)
            if terminal_index is None:
                skipped_partial += 1
                continue
            phosphate_index = terminal_index
            is_terminal = True
        indexed_groups.append((phosphate_index, residue, is_terminal))

    phosphate_donors = [(index, residue) for index, residue, terminal in indexed_groups if not terminal]
    if not phosphate_donors:
        raise RegularizePhosphatesError(
            f"Model {model_id}, chain {chain_label} has no complete nonterminal P/OP1/OP2 group to regularize from."
        )

    # Each atom is transformed back to helical index zero and averaged there.
    # P/OP1/OP2/O5'/C5' on the first nucleotide are excluded because the 5'
    # terminal phosphate is minimized differently. O3' on the last nucleotide
    # is excluded because its 3' terminal environment is also different. All
    # end atoms are populated from the corresponding internal consensus later.
    consensus: Dict[str, np.ndarray] = {}
    atom_targets: Dict[str, List[Tuple[int, AtomRecord]]] = {}
    complete_nucleotide_group_ids = {
        id(residue) for _index, residue, _terminal in indexed_groups if id(residue) in nucleotide_indices
    }
    for atom_name in REGULARIZED_ATOMS:
        targets: List[Tuple[int, AtomRecord]] = []
        donors: List[Tuple[int, AtomRecord]] = []
        for nucleotide_index, residue in enumerate(nucleotides):
            if atom_name in PHOSPHATE_ATOMS and id(residue) not in complete_nucleotide_group_ids:
                continue
            atom = residue.atoms.get(atom_name)
            if atom is None:
                continue
            targets.append((nucleotide_index, atom))
            is_terminal_environment = (
                nucleotide_index == 0 if atom_name != "O3'" else nucleotide_index == len(nucleotides) - 1
            )
            if not is_terminal_environment:
                donors.append((nucleotide_index, atom))
        if atom_name in PHOSPHATE_ATOMS:
            for phosphate_index, residue, _is_terminal in indexed_groups:
                if id(residue) in nucleotide_indices:
                    continue
                atom = residue.atoms.get(atom_name)
                if atom is not None:
                    targets.append((phosphate_index, atom))
        if not donors:
            continue
        atom_targets[atom_name] = targets
        reference_coords = []
        for atom_index, atom in donors:
            to_reference = np.linalg.matrix_power(step, -atom_index)
            reference_coords.append(_apply_transform(to_reference, atom.coord))
        consensus[atom_name] = np.mean(np.asarray(reference_coords), axis=0)

    regularized_atom_count = 0
    for atom_name, targets in atom_targets.items():
        for atom_index, atom in targets:
            from_reference = np.linalg.matrix_power(step, atom_index)
            output_coords[atom.line_index] = _apply_transform(from_reference, consensus[atom_name])
            regularized_atom_count += 1

    terminal_count = sum(int(is_terminal) for _index, _residue, is_terminal in indexed_groups)

    return ChainRegularization(
        model_id=model_id,
        chain_id=chain_id,
        nucleotide_count=len(nucleotides),
        fitted_step_rmsd=fit_rmsd,
        regularized_group_count=len(indexed_groups),
        regularized_atom_count=regularized_atom_count,
        terminal_group_count=terminal_count,
        skipped_partial_group_count=skipped_partial,
    )


def default_regularized_output_path(input_pdb: Union[Path, str]) -> Path:
    path = Path(input_pdb)
    suffix = path.suffix or ".pdb"
    return path.with_name(f"{path.stem}_regularized_phosphates{suffix}")


def _provenance_remark() -> str:
    return (
        f"REMARK BNP_NA_REGULARIZE_PHOSPHATES bnp_na {TOOL_VERSION}; "
        "C1'-DERIVED HELICAL SYMMETRY; TERMINAL ATOMS EXCLUDED FROM CONSENSUS.\n"
    )


def regularize_phosphates(
    input_pdb: Union[Path, str],
    output_pdb: Optional[Union[Path, str]] = None,
) -> RegularizePhosphatesResult:
    """Regularize all atoms selected by min_P_C5.params along each PDB chain."""
    input_path = Path(input_pdb).expanduser().resolve()
    output_path = (
        Path(output_pdb).expanduser().resolve()
        if output_pdb is not None
        else default_regularized_output_path(input_path)
    )
    lines, residues = _parse_pdb(input_path)

    chains: Dict[Tuple[int, str], List[Residue]] = {}
    for residue in residues:
        chains.setdefault((residue.model_id, residue.chain_id), []).append(residue)

    output_coords: Dict[int, np.ndarray] = {}
    reports = []
    for chain_residues in chains.values():
        # Ignore non-nucleic-acid chains without C1' atoms and phosphate groups.
        if not any("C1'" in residue.atoms for residue in chain_residues):
            continue
        if not any(any(name in residue.atoms for name in PHOSPHATE_ATOMS) for residue in chain_residues):
            continue
        reports.append(_regularize_chain(chain_residues, output_coords))

    if not reports:
        raise RegularizePhosphatesError("No chain containing both C1' atoms and phosphate groups was found.")

    for line_index, coord in output_coords.items():
        lines[line_index] = _format_xyz(lines[line_index], coord)
    remark = _provenance_remark()
    insert_at = next((index for index, line in enumerate(lines) if line[:6].strip() in {"ATOM", "HETATM", "MODEL"}), 0)
    lines.insert(insert_at, remark)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(lines), encoding="utf-8")

    log_lines = [
        "=== Regularize phosphates ===",
        f"Input : {input_path}",
        f"Output: {output_path}",
        "Method: fit one-residue screw transform from consecutive C1' atoms; average internal",
        "        P/OP1/OP2/O5'/C5'/O3' positions in a common helical frame; propagate each",
        "        consensus position to internal and terminal atoms.",
    ]
    for report in reports:
        chain_label = report.chain_id if report.chain_id.strip() else "<blank>"
        log_lines.append(
            f"Model {report.model_id}, chain {chain_label}: {report.nucleotide_count} nucleotides; "
            f"C1' step RMSD {report.fitted_step_rmsd:.6f} A; "
            f"{report.regularized_group_count} groups regularized "
            f"({report.terminal_group_count} terminal); "
            f"{report.regularized_atom_count} atom coordinates updated; "
            f"{report.skipped_partial_group_count} partial groups skipped."
        )

    return RegularizePhosphatesResult(
        input_pdb=input_path,
        output_pdb=output_path,
        chains=tuple(reports),
        log_text="\n".join(log_lines),
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Regularize P/OP1/OP2/O5'/C5'/O3' using helical symmetry fitted from each chain's C1' atoms."
        )
    )
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("-o", "--out", default=None, help="Output PDB (default: <input>_regularized_phosphates.pdb)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        result = regularize_phosphates(args.input_pdb, args.out)
    except RegularizePhosphatesError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(result.log_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
