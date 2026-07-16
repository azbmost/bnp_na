#!/usr/bin/env python3
"""Add missing terminal phosphate groups and their preceding O3' atoms."""
from __future__ import annotations

import argparse
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np


TOOL_VERSION = "V13.10"


class AddPhosphateError(Exception):
    """Raised when terminal phosphate analysis or placement fails."""


PHOSPHATE_ALIASES = {
    "P": "P",
    "OP1": "OP1",
    "O1P": "OP1",
    "OP2": "OP2",
    "O2P": "OP2",
}
PHOSPHATE_ATOMS = ("P", "OP1", "OP2")
PHOSPHATE_FRAME_ATOMS = PHOSPHATE_ATOMS + ("O5'",)
SUGAR_ALIGN_ATOMS = ("C1'", "C2'", "C3'", "C4'", "C5'", "O3'", "O4'", "O5'")
SUGAR_MARKER_ATOMS = set(SUGAR_ALIGN_ATOMS) | {"O2'"}


@dataclass
class AtomRecord:
    line: str
    index: int
    record_name: str
    serial: int
    name: str
    canonical_name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: np.ndarray


@dataclass
class Residue:
    chain: str
    resseq: int
    icode: str
    resname: str
    atoms: List[AtomRecord] = field(default_factory=list)

    @property
    def key(self) -> Tuple[str, int, str, str]:
        return self.chain, self.resseq, self.icode, self.resname

    @property
    def first_index(self) -> int:
        return min(atom.index for atom in self.atoms)

    @property
    def last_index(self) -> int:
        return max(atom.index for atom in self.atoms)

    def atom_by_name(self) -> Dict[str, AtomRecord]:
        out: Dict[str, AtomRecord] = {}
        for atom in self.atoms:
            out.setdefault(atom.canonical_name, atom)
        return out

    def phosphate_atoms(self) -> Dict[str, AtomRecord]:
        atoms = self.atom_by_name()
        return {name: atoms[name] for name in PHOSPHATE_ATOMS if name in atoms}

    def has_complete_phosphate(self) -> bool:
        atoms = self.phosphate_atoms()
        return all(name in atoms for name in PHOSPHATE_ATOMS)

    def is_nucleotide(self) -> bool:
        atoms = self.atom_by_name()
        return not self.is_o3prime_only() and any(name in atoms for name in SUGAR_MARKER_ATOMS)

    def is_o3prime_only(self) -> bool:
        atoms = self.atom_by_name()
        other_sugar_atoms = SUGAR_MARKER_ATOMS - {"O3'"}
        return (
            "O3'" in atoms
            and not any(name in atoms for name in other_sugar_atoms)
            and not self.has_complete_phosphate()
        )

    def is_phosphate_only(self) -> bool:
        atoms = self.atom_by_name()
        return self.has_complete_phosphate() and not any(name in atoms for name in SUGAR_MARKER_ATOMS)


@dataclass
class ChainPhosphateStatus:
    chain_id: str
    nucleotide_count: int
    first_residue: Optional[str]
    last_residue: Optional[str]
    has_5prime_phosphate: bool
    has_5prime_o3: bool
    has_3prime_phosphate: bool
    can_add_5prime: bool
    can_add_5prime_o3: bool
    can_add_3prime: bool
    notes: List[str] = field(default_factory=list)


@dataclass
class PhosphateAdditionResult:
    input_pdb: Path
    output_pdb: Path
    statuses_before: List[ChainPhosphateStatus]
    statuses_after: List[ChainPhosphateStatus]
    added: List[str]
    skipped: List[str]
    log_text: str


def _is_atom_line(line: str) -> bool:
    return line.startswith("ATOM") or line.startswith("HETATM")


def _canonical_atom_name(name: str) -> str:
    return PHOSPHATE_ALIASES.get(name.strip().upper(), name.strip().upper())


def _parse_atom_line(line: str, index: int) -> AtomRecord:
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21] if len(line) > 21 else " "
        resseq = int(line[22:26])
        icode = line[26] if len(line) > 26 else " "
        xyz = np.array(
            [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ],
            dtype=float,
        )
    except Exception as exc:
        raise AddPhosphateError(f"Could not parse ATOM/HETATM line {index + 1}: {line.rstrip()}") from exc
    return AtomRecord(
        line=line,
        index=index,
        record_name=line[:6].strip() or "ATOM",
        serial=serial,
        name=name,
        canonical_name=_canonical_atom_name(name),
        resname=resname,
        chain=chain,
        resseq=resseq,
        icode=icode,
        xyz=xyz,
    )


def _read_pdb(path: Path) -> Tuple[List[str], List[AtomRecord], List[Residue], Dict[str, List[Residue]]]:
    if not path.exists():
        raise AddPhosphateError(f"Input PDB not found: {path}")
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    atoms: List[AtomRecord] = []
    residues: List[Residue] = []
    residue_by_key: Dict[Tuple[str, int, str, str], Residue] = {}
    for index, line in enumerate(lines):
        if not _is_atom_line(line):
            continue
        atom = _parse_atom_line(line, index)
        atoms.append(atom)
        key = (atom.chain, atom.resseq, atom.icode, atom.resname)
        residue = residue_by_key.get(key)
        if residue is None:
            residue = Residue(atom.chain, atom.resseq, atom.icode, atom.resname)
            residue_by_key[key] = residue
            residues.append(residue)
        residue.atoms.append(atom)
    if not atoms:
        raise AddPhosphateError(f"No ATOM/HETATM records found in {path}")
    chains: Dict[str, List[Residue]] = {}
    for residue in residues:
        chains.setdefault(residue.chain, []).append(residue)
    return lines, atoms, residues, chains


def _max_pdb_record_serial(lines: Sequence[str]) -> int:
    max_serial = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "TER")):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                continue
    return max_serial


def _insert_provenance_remark(lines: List[str], phosphate_count: int, o3_count: int) -> List[str]:
    if phosphate_count <= 0 and o3_count <= 0:
        return lines
    additions: List[str] = []
    if phosphate_count:
        additions.append(f"{phosphate_count} terminal phosphate group(s)")
    if o3_count:
        additions.append(f"{o3_count} preceding O3' atom(s)")
    remark = f"REMARK BNP_NA_ADD_PHOSPHATES bnp_na {TOOL_VERSION} added {' and '.join(additions)}.\n"
    if remark in lines:
        return lines
    insert_at = 0
    while insert_at < len(lines) and lines[insert_at].startswith("REMARK"):
        insert_at += 1
    return lines[:insert_at] + [remark] + lines[insert_at:]


def _line_serial(line: str) -> Optional[int]:
    try:
        return int(line[6:11])
    except ValueError:
        return None


def _replace_line_serial(line: str, serial: int) -> str:
    if serial > 99999:
        raise AddPhosphateError("PDB atom serial number would exceed 99999 after renumbering.")
    padded = line.rstrip("\n").ljust(11)
    return padded[:6] + f"{serial:5d}" + padded[11:] + ("\n" if line.endswith("\n") else "")


def _renumber_pdb_records(lines: List[str]) -> List[str]:
    """Renumber ATOM/HETATM/TER records in file order and update CONECT refs."""
    serial_map: Dict[int, int] = {}
    renumbered: List[str] = []
    next_serial = 1
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "TER")):
            old_serial = _line_serial(line)
            new_line = _replace_line_serial(line, next_serial)
            if old_serial is not None:
                if not line.startswith("TER") or old_serial not in serial_map:
                    serial_map[old_serial] = next_serial
            renumbered.append(new_line)
            next_serial += 1
        else:
            renumbered.append(line)

    updated: List[str] = []
    for line in renumbered:
        if not line.startswith("CONECT"):
            updated.append(line)
            continue
        numbers: List[int] = []
        for part in line[6:].split():
            try:
                numbers.append(int(part))
            except ValueError:
                pass
        if not numbers:
            updated.append(line)
            continue
        mapped = [serial_map.get(number, number) for number in numbers]
        updated.append("CONECT" + "".join(f"{number:5d}" for number in mapped) + "\n")
    return updated


def _chain_display(chain: str) -> str:
    return chain if chain.strip() else "(blank)"


def _residue_label(residue: Optional[Residue]) -> Optional[str]:
    if residue is None:
        return None
    icode = residue.icode.strip()
    suffix = icode if icode else ""
    return f"{_chain_display(residue.chain)}:{residue.resname}{residue.resseq}{suffix}"


def _nucleotide_residues(chain_residues: Sequence[Residue]) -> List[Residue]:
    return [residue for residue in chain_residues if residue.is_nucleotide()]


def _terminal_3prime_phosphate_residue(chain_residues: Sequence[Residue], last_nt: Residue) -> Optional[Residue]:
    after_last = [residue for residue in chain_residues if residue.first_index > last_nt.last_index]
    for residue in after_last:
        if residue.is_phosphate_only():
            return residue
    return None


def _terminal_5prime_o3_residue(chain_residues: Sequence[Residue], first_nt: Residue) -> Optional[Residue]:
    expected_resseq = first_nt.resseq - 1
    for residue in chain_residues:
        if residue.first_index >= first_nt.first_index:
            continue
        if residue.resseq == expected_resseq and "O3'" in residue.atom_by_name():
            return residue
    return None


def _common_sugar_atom_names(source: Residue, target: Residue) -> List[str]:
    source_atoms = source.atom_by_name()
    target_atoms = target.atom_by_name()
    return [name for name in SUGAR_ALIGN_ATOMS if name in source_atoms and name in target_atoms]


def _can_fit(source: Residue, target: Residue) -> bool:
    return len(_common_sugar_atom_names(source, target)) >= 3


def _common_phosphate_frame_atom_names(source: Residue, target: Residue) -> List[str]:
    source_atoms = source.atom_by_name()
    target_atoms = target.atom_by_name()
    return [name for name in PHOSPHATE_FRAME_ATOMS if name in source_atoms and name in target_atoms]


def _can_fit_phosphate_frame(source: Residue, target: Residue) -> bool:
    return len(_common_phosphate_frame_atom_names(source, target)) >= 3


def _status_for_chain(chain_id: str, chain_residues: Sequence[Residue]) -> ChainPhosphateStatus:
    nts = _nucleotide_residues(chain_residues)
    notes: List[str] = []
    if not nts:
        return ChainPhosphateStatus(
            chain_id=chain_id,
            nucleotide_count=0,
            first_residue=None,
            last_residue=None,
            has_5prime_phosphate=False,
            has_5prime_o3=False,
            has_3prime_phosphate=False,
            can_add_5prime=False,
            can_add_5prime_o3=False,
            can_add_3prime=False,
            notes=["No nucleotide residues with sugar atoms were detected."],
        )

    first = nts[0]
    last = nts[-1]
    has_5 = first.has_complete_phosphate()
    has_5_o3 = _terminal_5prime_o3_residue(chain_residues, first) is not None
    has_3 = _terminal_3prime_phosphate_residue(chain_residues, last) is not None
    can_add_5 = False
    can_add_5_o3 = False
    can_add_3 = False

    if len(nts) < 2:
        notes.append("At least two nucleotide residues are required for neighbor-based phosphate placement.")
    else:
        second = nts[1]
        previous = nts[-2]
        can_add_5 = (not has_5) and second.has_complete_phosphate() and _can_fit(second, first)
        first_o3 = first.atom_by_name().get("O3'")
        if not has_5_o3 and first_o3 is not None:
            if has_5:
                can_add_5_o3 = second.has_complete_phosphate() and _can_fit_phosphate_frame(second, first)
            else:
                can_add_5_o3 = can_add_5
        can_add_3 = (not has_3) and last.has_complete_phosphate() and _can_fit(previous, last)
        if not has_5 and not second.has_complete_phosphate():
            notes.append(f"Cannot add 5' phosphate: neighbor {_residue_label(second)} lacks P/OP1/OP2.")
        if not has_3 and not last.has_complete_phosphate():
            notes.append(f"Cannot add 3' phosphate: terminal {_residue_label(last)} lacks P/OP1/OP2.")
        if not has_5 and second.has_complete_phosphate() and not _can_fit(second, first):
            notes.append("Cannot add 5' phosphate: first two residues share fewer than three sugar alignment atoms.")
        if not has_5_o3 and first_o3 is None:
            notes.append("Cannot add preceding O3': the first nucleotide lacks a donor O3' atom.")
        elif not has_5_o3 and has_5 and not can_add_5_o3:
            notes.append("Cannot add preceding O3': the first two residues lack a compatible phosphate frame.")
        elif not has_5_o3 and not has_5 and can_add_5_o3:
            notes.append("The preceding O3' can be added together with the missing 5' phosphate.")
        if not has_3 and last.has_complete_phosphate() and not _can_fit(previous, last):
            notes.append("Cannot add 3' phosphate: last two residues share fewer than three sugar alignment atoms.")

    return ChainPhosphateStatus(
        chain_id=chain_id,
        nucleotide_count=len(nts),
        first_residue=_residue_label(first),
        last_residue=_residue_label(last),
        has_5prime_phosphate=has_5,
        has_5prime_o3=has_5_o3,
        has_3prime_phosphate=has_3,
        can_add_5prime=can_add_5,
        can_add_5prime_o3=can_add_5_o3,
        can_add_3prime=can_add_3,
        notes=notes,
    )


def analyze_phosphate_termini(pdb_path: Union[Path, str]) -> List[ChainPhosphateStatus]:
    """Return 5' and 3' terminal phosphate status for each detected chain."""
    _lines, _atoms, _residues, chains = _read_pdb(Path(pdb_path).expanduser())
    return [_status_for_chain(chain_id, residues) for chain_id, residues in chains.items()]


def format_phosphate_report(statuses: Sequence[ChainPhosphateStatus]) -> str:
    lines = ["=== Terminal phosphate status ==="]
    if not statuses:
        lines.append("No chains were detected.")
        return "\n".join(lines)
    for status in statuses:
        lines.append(
            f"Chain {_chain_display(status.chain_id)}: {status.nucleotide_count} nucleotide residues"
        )
        if status.first_residue and status.last_residue:
            lines.append(f"  First nucleotide: {status.first_residue}")
            lines.append(f"  Last nucleotide : {status.last_residue}")
        lines.append(
            "  5' phosphate  : "
            + ("present" if status.has_5prime_phosphate else "missing")
            + ("; can add" if status.can_add_5prime else "")
        )
        lines.append(
            "  Preceding O3' : "
            + ("present" if status.has_5prime_o3 else "missing")
            + ("; can add" if status.can_add_5prime_o3 else "")
        )
        lines.append(
            "  3' phosphate  : "
            + ("present" if status.has_3prime_phosphate else "missing")
            + ("; can add" if status.can_add_3prime else "")
        )
        for note in status.notes:
            lines.append(f"  Note: {note}")
    return "\n".join(lines)


def parse_chain_selection(text: str) -> Optional[List[str]]:
    value = (text or "").strip()
    if not value:
        return None
    parts = [part for part in value.replace(",", " ").split() if part]
    chains: List[str] = []
    for part in parts:
        if part.lower() in {"blank", "(blank)", "_"}:
            chain = " "
        elif len(part) == 1:
            chain = part
        else:
            raise AddPhosphateError(f"Chain ID {part!r} is not a single-character PDB chain ID.")
        if chain not in chains:
            chains.append(chain)
    return chains


def default_add_phosphate_output_path(input_pdb: Union[Path, str]) -> Path:
    path = Path(input_pdb).expanduser()
    suffix = path.suffix or ".pdb"
    return path.with_name(f"{path.stem}_add_phosphates{suffix}")


def _fit_transform(source_points: np.ndarray, target_points: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    if source_points.shape != target_points.shape or source_points.shape[0] < 3:
        raise AddPhosphateError("Rigid fit requires at least three paired source/target points.")
    source_centroid = source_points.mean(axis=0)
    target_centroid = target_points.mean(axis=0)
    source_centered = source_points - source_centroid
    target_centered = target_points - target_centroid
    cov = source_centered.T @ target_centered
    u, _s, vt = np.linalg.svd(cov)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = target_centroid - rotation @ source_centroid
    fitted = (rotation @ source_points.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target_points) ** 2, axis=1))))
    return rotation, translation, rmsd


def _fit_residue_sugars(source: Residue, target: Residue) -> Tuple[np.ndarray, np.ndarray, float, List[str]]:
    names = _common_sugar_atom_names(source, target)
    if len(names) < 3:
        raise AddPhosphateError(
            f"Residues {_residue_label(source)} and {_residue_label(target)} share fewer than three sugar atoms."
        )
    source_atoms = source.atom_by_name()
    target_atoms = target.atom_by_name()
    source_points = np.array([source_atoms[name].xyz for name in names], dtype=float)
    target_points = np.array([target_atoms[name].xyz for name in names], dtype=float)
    rotation, translation, rmsd = _fit_transform(source_points, target_points)
    return rotation, translation, rmsd, names


def _fit_phosphate_frames(source: Residue, target: Residue) -> Tuple[np.ndarray, np.ndarray, float, List[str]]:
    names = _common_phosphate_frame_atom_names(source, target)
    if len(names) < 3:
        raise AddPhosphateError(
            f"Residues {_residue_label(source)} and {_residue_label(target)} share fewer than three "
            "phosphate-frame atoms."
        )
    source_atoms = source.atom_by_name()
    target_atoms = target.atom_by_name()
    source_points = np.array([source_atoms[name].xyz for name in names], dtype=float)
    target_points = np.array([target_atoms[name].xyz for name in names], dtype=float)
    rotation, translation, rmsd = _fit_transform(source_points, target_points)
    return rotation, translation, rmsd, names


def _transform_xyz(xyz: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return rotation @ xyz + translation


def _atom_field(atom_name: str) -> str:
    atom = atom_name.strip()
    if len(atom) == 1:
        return f" {atom:<3s}"
    if len(atom) < 4:
        return f" {atom:<3s}"[:4]
    return atom[:4]


def _element_for_atom(atom_name: str) -> str:
    for char in atom_name.strip():
        if char.isalpha():
            return char.upper()
    return ""


def _updated_atom_line(
    template: AtomRecord,
    *,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    icode: str,
    xyz: np.ndarray,
) -> str:
    if serial > 99999:
        raise AddPhosphateError("PDB atom serial number would exceed 99999.")
    if not (-999 <= resseq <= 9999):
        raise AddPhosphateError(f"Residue number {resseq} cannot be written in PDB columns 23-26.")
    base = template.line.rstrip("\n").ljust(80)
    record = (template.record_name or "ATOM")[:6].ljust(6)
    res_field = resname.strip()[:3].rjust(3)
    chain_field = (chain or " ")[:1]
    icode_field = (icode or " ")[:1]
    element = _element_for_atom(atom_name)
    new_line = (
        record
        + f"{serial:5d}"
        + base[11:12]
        + _atom_field(atom_name)
        + base[16:17]
        + res_field
        + base[20:21]
        + chain_field
        + f"{resseq:4d}"
        + icode_field
        + base[27:30]
        + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        + base[54:]
    )
    if len(new_line) < 80:
        new_line = new_line.ljust(80)
    new_line = new_line[:76] + element.rjust(2) + new_line[78:]
    return new_line.rstrip() + "\n"


def _next_resseq(last: Residue, chain_residues: Sequence[Residue]) -> int:
    candidate = last.resseq + 1
    occupied = {(residue.resseq, residue.icode.strip()) for residue in chain_residues}
    if (candidate, "") in occupied:
        raise AddPhosphateError(
            f"Cannot create 3' phosphate residue after {_residue_label(last)}: residue number {candidate} already exists."
        )
    return candidate


def _previous_resseq(first: Residue, chain_residues: Sequence[Residue]) -> int:
    candidate = first.resseq - 1
    occupied = {(residue.resseq, residue.icode.strip()) for residue in chain_residues}
    if (candidate, "") in occupied:
        raise AddPhosphateError(
            f"Cannot create preceding O3' residue before {_residue_label(first)}: "
            f"residue number {candidate} already exists."
        )
    return candidate


def _new_phosphate_lines(
    donor_atoms: Dict[str, AtomRecord],
    *,
    rotation: np.ndarray,
    translation: np.ndarray,
    serial_start: int,
    target_residue: Residue,
    target_resseq: Optional[int] = None,
) -> List[str]:
    lines: List[str] = []
    for offset, atom_name in enumerate(PHOSPHATE_ATOMS):
        donor = donor_atoms[atom_name]
        xyz = _transform_xyz(donor.xyz, rotation, translation)
        lines.append(
            _updated_atom_line(
                donor,
                serial=serial_start + offset,
                atom_name=atom_name,
                resname=target_residue.resname,
                chain=target_residue.chain,
                resseq=target_residue.resseq if target_resseq is None else target_resseq,
                icode=target_residue.icode if target_resseq is None else " ",
                xyz=xyz,
            )
        )
    return lines


def _new_5prime_o3_line(
    donor: AtomRecord,
    *,
    rotation: np.ndarray,
    translation: np.ndarray,
    serial: int,
    target_residue: Residue,
    target_resseq: int,
) -> str:
    return _updated_atom_line(
        donor,
        serial=serial,
        atom_name="O3'",
        resname=target_residue.resname,
        chain=target_residue.chain,
        resseq=target_resseq,
        icode=" ",
        xyz=_transform_xyz(donor.xyz, rotation, translation),
    )


def add_terminal_phosphates(
    input_pdb: Union[Path, str],
    output_pdb: Union[Path, str],
    *,
    chain_ids: Optional[Sequence[str]] = None,
    add_5prime: bool = True,
    add_3prime: bool = True,
    add_5prime_o3: bool = False,
) -> PhosphateAdditionResult:
    """Add selected terminal phosphates and optionally the O3' preceding a 5' phosphate."""
    if not add_5prime and not add_3prime and not add_5prime_o3:
        raise AddPhosphateError("Select at least one item to add: 5' phosphate, 3' phosphate, or preceding O3'.")
    input_path = Path(input_pdb).expanduser()
    output_path = Path(output_pdb).expanduser()
    lines, atoms, _residues, chains = _read_pdb(input_path)
    statuses_before = [_status_for_chain(chain_id, residues) for chain_id, residues in chains.items()]

    selected_chains = list(chain_ids) if chain_ids is not None else list(chains.keys())
    unknown = [chain for chain in selected_chains if chain not in chains]
    if unknown:
        names = ", ".join(_chain_display(chain) for chain in unknown)
        raise AddPhosphateError(f"Selected chain(s) not found in the PDB: {names}")

    before_insertions: Dict[int, List[str]] = {}
    after_insertions: Dict[int, List[str]] = {}
    added: List[str] = []
    skipped: List[str] = []
    phosphate_added_count = 0
    o3_added_count = 0
    next_serial = max(max(atom.serial for atom in atoms), _max_pdb_record_serial(lines)) + 1

    for chain_id in selected_chains:
        chain_residues = chains[chain_id]
        nts = _nucleotide_residues(chain_residues)
        if len(nts) < 2:
            skipped.append(f"Chain {_chain_display(chain_id)}: fewer than two nucleotide residues.")
            continue
        first = nts[0]
        second = nts[1]
        previous = nts[-2]
        last = nts[-1]
        five_prime_transform: Optional[Tuple[np.ndarray, np.ndarray, float, List[str]]] = None
        five_prime_present = first.has_complete_phosphate()

        if add_5prime:
            if five_prime_present:
                skipped.append(f"Chain {_chain_display(chain_id)} 5': already present.")
            else:
                donor_phosphate = second.phosphate_atoms()
                if not all(name in donor_phosphate for name in PHOSPHATE_ATOMS):
                    skipped.append(f"Chain {_chain_display(chain_id)} 5': neighbor lacks P/OP1/OP2.")
                else:
                    rotation, translation, rmsd, names = _fit_residue_sugars(second, first)
                    new_lines = _new_phosphate_lines(
                        donor_phosphate,
                        rotation=rotation,
                        translation=translation,
                        serial_start=next_serial,
                        target_residue=first,
                    )
                    next_serial += len(new_lines)
                    before_insertions.setdefault(first.first_index, []).extend(new_lines)
                    five_prime_transform = (rotation, translation, rmsd, names)
                    five_prime_present = True
                    phosphate_added_count += 1
                    added.append(
                        f"Chain {_chain_display(chain_id)} 5': added P/OP1/OP2 to {_residue_label(first)} "
                        f"from {_residue_label(second)} after sugar fit RMSD {rmsd:.4f} A "
                        f"({', '.join(names)})."
                    )

        if add_5prime_o3:
            existing_o3 = _terminal_5prime_o3_residue(chain_residues, first)
            if existing_o3 is not None:
                skipped.append(
                    f"Chain {_chain_display(chain_id)} preceding O3': already present in {_residue_label(existing_o3)}."
                )
            elif not five_prime_present:
                skipped.append(
                    f"Chain {_chain_display(chain_id)} preceding O3': requires an existing or newly added 5' phosphate."
                )
            else:
                donor_o3 = first.atom_by_name().get("O3'")
                if donor_o3 is None:
                    skipped.append(f"Chain {_chain_display(chain_id)} preceding O3': first nucleotide lacks O3'.")
                else:
                    if five_prime_transform is None:
                        rotation, translation, rmsd, names = _fit_phosphate_frames(second, first)
                        fit_description = "phosphate-frame fit"
                    else:
                        rotation, translation, rmsd, names = five_prime_transform
                        fit_description = "sugar fit used for the new 5' phosphate"
                    new_resseq = _previous_resseq(first, chain_residues)
                    new_line = _new_5prime_o3_line(
                        donor_o3,
                        rotation=rotation,
                        translation=translation,
                        serial=next_serial,
                        target_residue=first,
                        target_resseq=new_resseq,
                    )
                    next_serial += 1
                    before_insertions.setdefault(first.first_index, []).insert(0, new_line)
                    o3_added_count += 1
                    added.append(
                        f"Chain {_chain_display(chain_id)} 5': added O3' as "
                        f"{_chain_display(chain_id)}:{first.resname}{new_resseq}, preceding {_residue_label(first)}, "
                        f"after {fit_description} RMSD {rmsd:.4f} A ({', '.join(names)})."
                    )

        if add_3prime:
            if _terminal_3prime_phosphate_residue(chain_residues, last) is not None:
                skipped.append(f"Chain {_chain_display(chain_id)} 3': already present.")
            else:
                donor_phosphate = last.phosphate_atoms()
                if not all(name in donor_phosphate for name in PHOSPHATE_ATOMS):
                    skipped.append(f"Chain {_chain_display(chain_id)} 3': terminal residue lacks P/OP1/OP2.")
                else:
                    rotation, translation, rmsd, names = _fit_residue_sugars(previous, last)
                    new_resseq = _next_resseq(last, chain_residues)
                    new_lines = _new_phosphate_lines(
                        donor_phosphate,
                        rotation=rotation,
                        translation=translation,
                        serial_start=next_serial,
                        target_residue=last,
                        target_resseq=new_resseq,
                    )
                    next_serial += len(new_lines)
                    after_insertions.setdefault(last.last_index, []).extend(new_lines)
                    phosphate_added_count += 1
                    added.append(
                        f"Chain {_chain_display(chain_id)} 3': added phosphate-only residue "
                        f"{_chain_display(chain_id)}:{last.resname}{new_resseq} from {_residue_label(last)} "
                        f"after {_residue_label(previous)}->{_residue_label(last)} sugar fit RMSD {rmsd:.4f} A "
                        f"({', '.join(names)})."
                    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_lines: List[str] = []
    for index, line in enumerate(lines):
        out_lines.extend(before_insertions.get(index, []))
        out_lines.append(line)
        out_lines.extend(after_insertions.get(index, []))
    out_lines = _insert_provenance_remark(out_lines, phosphate_added_count, o3_added_count)
    out_lines = _renumber_pdb_records(out_lines)
    output_path.write_text("".join(out_lines), encoding="utf-8")

    statuses_after = analyze_phosphate_termini(output_path)
    log_parts = [
        "=== Add phosphates ===",
        f"Input PDB : {input_path}",
        f"Output PDB: {output_path}",
        "",
        format_phosphate_report(statuses_before),
        "",
        "Additions:",
    ]
    log_parts.extend([f"  {line}" for line in added] or ["  None"])
    if skipped:
        log_parts.append("")
        log_parts.append("Skipped:")
        log_parts.extend(f"  {line}" for line in skipped)
    log_parts.extend(["", "After writing:", format_phosphate_report(statuses_after)])

    return PhosphateAdditionResult(
        input_pdb=input_path,
        output_pdb=output_path,
        statuses_before=statuses_before,
        statuses_after=statuses_after,
        added=added,
        skipped=skipped,
        log_text="\n".join(log_parts),
    )


def _ends_from_text(text: str) -> Tuple[bool, bool]:
    value = (text or "both").strip().lower()
    if value in {"both", "all", "5,3", "3,5"}:
        return True, True
    if value in {"5", "5prime", "5'", "five"}:
        return True, False
    if value in {"3", "3prime", "3'", "three"}:
        return False, True
    if value in {"none", "neither"}:
        return False, False
    raise AddPhosphateError("--ends must be one of: both, 5, 3, none")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Report and add terminal 5'/3' phosphate groups and optionally the O3' atom preceding a 5' phosphate."
        )
    )
    parser.add_argument("input_pdb", help="Input PDB file.")
    parser.add_argument("-o", "--out", default=None, help="Output PDB file. Default: <input>_add_phosphates.pdb")
    parser.add_argument("--chains", default="", help="Chain IDs to edit, separated by spaces or commas. Blank means all chains.")
    parser.add_argument("--ends", default="both", help="Ends to add: both, 5, 3, or none. Default: both.")
    parser.add_argument(
        "--add-5prime-o3",
        action="store_true",
        help="Add O3' as residue n-1 before an existing or newly added 5' phosphate on residue n.",
    )
    parser.add_argument("--report-only", action="store_true", help="Only report terminal phosphate status; do not write a PDB.")
    args = parser.parse_args(argv)

    input_path = Path(args.input_pdb).expanduser()
    if args.report_only:
        print(format_phosphate_report(analyze_phosphate_termini(input_path)))
        return 0

    add_5, add_3 = _ends_from_text(args.ends)
    output_path = Path(args.out).expanduser() if args.out else default_add_phosphate_output_path(input_path)
    chains = parse_chain_selection(args.chains)
    result = add_terminal_phosphates(
        input_path,
        output_path,
        chain_ids=chains,
        add_5prime=add_5,
        add_3prime=add_3,
        add_5prime_o3=args.add_5prime_o3,
    )
    print(result.log_text)
    print()
    cli = ["python3", "bnp_na_lib/add_phosphates.py", str(input_path), "-o", str(output_path), "--ends", args.ends]
    if chains:
        cli.extend(["--chains", args.chains])
    if args.add_5prime_o3:
        cli.append("--add-5prime-o3")
    print("Equivalent CLI command:")
    print("  " + " ".join(shlex.quote(part) for part in cli))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
