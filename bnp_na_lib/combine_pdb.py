#!/usr/bin/env python3
"""Combine PDB files while assigning consecutive, unique chain IDs."""
from __future__ import annotations

import argparse
import re
import string
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple, Union


TOOL_VERSION = "V13.17"
CHAIN_IDS = string.ascii_uppercase
ATOM_RECORDS = ("ATOM  ", "HETATM")
SERIAL_RECORDS = ATOM_RECORDS + ("TER   ",)
ATOM_COMPANION_RECORDS = ("ANISOU", "SIGATM", "SIGUIJ")
ALL_CHAINS_WORDS = frozenset({"all", "*", "-"})
BLANK_CHAIN_WORDS = frozenset({"blank", "(blank)", "_"})

# One per input file: a string such as "A B", an explicit list of chain IDs, or
# None/blank/"all" for every chain in that file.
ChainSelection = Union[str, Sequence[str], None]


class CombinePDBError(Exception):
    """Raised when PDB inputs cannot be combined safely."""


@dataclass(frozen=True)
class InputChainMapping:
    input_pdb: Path
    chain_map: Dict[str, str]
    atom_count: int
    excluded_chains: Tuple[str, ...] = ()


@dataclass(frozen=True)
class CombinePDBResult:
    input_pdbs: List[Path]
    output_pdb: Path
    mappings: List[InputChainMapping]
    atom_count: int
    chain_count: int
    remark_count: int
    link_count: int
    log_text: str


def _record_name(line: str) -> str:
    return line[:6].ljust(6)


def _line_chain(line: str) -> str:
    return line[21] if len(line) > 21 else " "


def format_chain_id(chain_id: str) -> str:
    """Return a human-readable label for a PDB chain ID."""
    return chain_id if chain_id.strip() else "(blank)"


def format_chain_token(chain_id: str) -> str:
    """Return a chain ID in the form accepted by ``parse_chain_selection``."""
    return "_" if not chain_id.strip() else chain_id


def _replace_columns(line: str, start: int, end: int, value: str) -> str:
    body = line.rstrip("\r\n")
    body = body.ljust(end)
    return body[:start] + value + body[end:] + "\n"


def _parse_serial(line: str, path: Path, line_number: int) -> int:
    try:
        return int(line[6:11])
    except (ValueError, IndexError) as exc:
        raise CombinePDBError(
            f"Could not read the atom serial on line {line_number} of {path}: {line.rstrip()}"
        ) from exc


def _set_serial(line: str, serial: int) -> str:
    if serial > 99999:
        raise CombinePDBError("The combined PDB would exceed the 99999-record PDB serial-number limit.")
    return _replace_columns(line, 6, 11, f"{serial:5d}")


def _set_chain(line: str, chain_id: str) -> str:
    return _replace_columns(line, 21, 22, chain_id)


class _ChainRemapper:
    """Rename source chain IDs and remember references to excluded chains.

    ``saw_excluded`` lets the REMARK handler drop a record whose chain no longer
    exists in the combined file instead of leaving a stale ID behind that now
    names a different chain.
    """

    def __init__(self, chain_map: Dict[str, str], excluded_chains: Set[str]) -> None:
        self.chain_map = chain_map
        self.excluded_chains = excluded_chains
        self.saw_excluded = False

    def reset(self) -> None:
        self.saw_excluded = False

    def __call__(self, chain_id: str) -> str:
        source_chain = " " if chain_id == "_" else chain_id
        if source_chain in self.excluded_chains:
            self.saw_excluded = True
            return chain_id
        return self.chain_map.get(source_chain, chain_id)


def _replace_remark_field(line: str, field: str, replacement: Callable[[str], str]) -> str:
    pattern = re.compile(rf"(?<!\S)({re.escape(field)}=)(\S+)")

    def replace(match: re.Match[str]) -> str:
        return match.group(1) + replacement(match.group(2))

    return pattern.sub(replace, line)


def _replace_residue_labels(text: str, remap: _ChainRemapper) -> str:
    """Update parse-friendly ``A:12:DA`` and DSSR-style ``A.DA12`` labels."""
    colon_pattern = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z0-9_])(?=:[+-]?\d+(?::[A-Za-z0-9]+)?)")
    dot_pattern = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z0-9_])(?=\.[A-Za-z0-9]{1,3}[+-]?\d+)")

    def replace(match: re.Match[str]) -> str:
        return remap(match.group(1))

    return dot_pattern.sub(replace, colon_pattern.sub(replace, text))


def _update_remark_record(line: str, remap: _ChainRemapper) -> Optional[str]:
    """Update chain references in known structured and common REMARK formats.

    Returns ``None`` when the record describes a chain that was not selected.
    """
    if line.startswith("REMARK 950 RE_SCRIPT COMMAND"):
        return line if line.endswith(("\n", "\r")) else line + "\n"

    remap.reset()
    if line.startswith("REMARK 950 RE_SCRIPT CHAIN_RANGE"):
        updated = _replace_remark_field(line, "chain", remap)
        for field in ("start", "end"):
            updated = _replace_remark_field(updated, field, lambda value: _replace_residue_labels(value, remap))
    elif line.startswith("REMARK 950 RE_SCRIPT CHAIN_RESIDUES"):
        updated = _replace_remark_field(line, "chain", remap)
        updated = _replace_remark_field(
            updated, "residues", lambda value: _replace_residue_labels(value, remap)
        )
    elif line.startswith("REMARK 950 RE_SCRIPT JUNCTION"):
        updated = line
        for field in ("residues", "core", "excluded_nick_ends"):
            updated = _replace_remark_field(
                updated, field, lambda value: _replace_residue_labels(value, remap)
            )
    elif line.startswith("REMARK 950 RE_SCRIPT SPECIAL"):
        updated = _replace_remark_field(line, "chain", remap)
        updated = _replace_remark_field(
            updated, "residue", lambda value: _replace_residue_labels(value, remap)
        )
    elif line.startswith("REMARK 950 RE_SCRIPT"):
        # Unrecognized RE_SCRIPT records are provenance. Preserve them exactly
        # rather than risk changing original/source residue identities.
        updated = line
    else:
        updated = _replace_remark_field(line, "chain", remap)
        updated = re.sub(
            r"(?<!\S)(CHAIN\s+)([A-Za-z0-9_])(?=\s|$)",
            lambda match: match.group(1) + remap(match.group(2)),
            updated,
        )
        updated = _replace_residue_labels(updated, remap)

    if remap.saw_excluded:
        return None
    return updated.rstrip("\r\n") + "\n"


def _update_link_record(
    line: str, chain_map: Dict[str, str], excluded_chains: Set[str]
) -> Optional[str]:
    """Update both fixed-column chain IDs in a PDB LINK record.

    Returns ``None`` when either endpoint sits on a chain that was not selected.
    """
    updated = line
    for index in (21, 51):
        old_chain = updated[index] if len(updated.rstrip("\r\n")) > index else " "
        if old_chain in excluded_chains:
            return None
        if old_chain in chain_map:
            updated = _replace_columns(updated, index, index + 1, chain_map[old_chain])
    return updated.rstrip("\r\n") + "\n"


def _update_het_record(
    line: str, chain_map: Dict[str, str], excluded_chains: Set[str]
) -> Optional[str]:
    """Update the fixed-column chain ID in a PDB HET residue record."""
    old_chain = line[12] if len(line.rstrip("\r\n")) > 12 else " "
    if old_chain in excluded_chains:
        return None
    if old_chain in chain_map:
        return _replace_columns(line, 12, 13, chain_map[old_chain])
    return line.rstrip("\r\n") + "\n"


def _ordered_chains(lines: Sequence[str], path: Path) -> Tuple[List[str], Dict[str, int]]:
    chains: List[str] = []
    atom_counts: Dict[str, int] = {}
    for line in lines:
        if _record_name(line) not in ATOM_RECORDS:
            continue
        chain_id = _line_chain(line)
        if chain_id not in atom_counts:
            chains.append(chain_id)
            atom_counts[chain_id] = 0
        atom_counts[chain_id] += 1
    if not chains:
        raise CombinePDBError(f"No ATOM/HETATM records were found in {path}.")
    return chains, atom_counts


def _ter_chain(line: str, source_chains: Set[str], previous_atom_chain: str) -> str:
    """Return the chain a TER record belongs to, tolerating a blank chain column."""
    chain_id = _line_chain(line)
    if chain_id in source_chains:
        return chain_id
    return previous_atom_chain


def _conect_serials(line: str, path: Path, line_number: int) -> List[int]:
    values: List[int] = []
    for field_start in range(6, len(line.rstrip("\r\n")), 5):
        field = line[field_start : field_start + 5]
        if not field.strip():
            continue
        try:
            values.append(int(field))
        except ValueError as exc:
            raise CombinePDBError(
                f"Could not read CONECT record on line {line_number} of {path}: {line.rstrip()}"
            ) from exc
    return values


def _normalize_chain_token(part: str) -> str:
    if part == " ":
        return " "
    token = part.strip()
    if token.lower() in BLANK_CHAIN_WORDS:
        return " "
    if len(token) == 1:
        return token
    raise CombinePDBError(f"Chain ID {part!r} is not a single-character PDB chain ID.")


def parse_chain_selection(text: str) -> Optional[List[str]]:
    """Parse ``"A B"``/``"A,B"`` into chain IDs; blank or ``all`` means every chain.

    Use ``_`` or ``blank`` to name a chain whose PDB chain column is empty.
    """
    value = (text or "").strip()
    if not value or value.lower() in ALL_CHAINS_WORDS:
        return None
    chains: List[str] = []
    for part in value.replace(",", " ").split():
        chain = _normalize_chain_token(part)
        if chain not in chains:
            chains.append(chain)
    return chains or None


def _normalize_selection(value: ChainSelection) -> Optional[List[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        return parse_chain_selection(value)
    chains: List[str] = []
    for part in value:
        chain = _normalize_chain_token(str(part))
        if chain not in chains:
            chains.append(chain)
    return chains or None


def _select_chains(path: Path, chains: Sequence[str], selection: Optional[Sequence[str]]) -> List[str]:
    """Keep the requested chains of one input, in their first-appearance order."""
    if selection is None:
        return list(chains)
    available = set(chains)
    missing = [chain for chain in selection if chain not in available]
    if missing:
        raise CombinePDBError(
            f"{path} has no chain "
            + ", ".join(format_chain_id(chain) for chain in missing)
            + "; its chains are "
            + ", ".join(format_chain_id(chain) for chain in chains)
            + "."
        )
    requested = set(selection)
    return [chain for chain in chains if chain in requested]


def list_pdb_chains(pdb_path: Union[Path, str]) -> List[str]:
    """Return a PDB's ATOM/HETATM chain IDs in first-appearance order."""
    path = Path(pdb_path).expanduser()
    if not path.is_file():
        raise CombinePDBError(f"Input PDB not found: {path}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    chains, _atom_counts = _ordered_chains(lines, path)
    return chains


def default_combine_pdb_output_path(output_dir: Union[Path, str]) -> Path:
    """Return the GUI's default output path for the Combine_PDB tool."""
    return Path(output_dir).expanduser() / "combine_PDB_out.pdb"


def combine_pdb_files(
    input_pdbs: Sequence[Union[Path, str]],
    output_pdb: Union[Path, str],
    chain_selections: Optional[Sequence[ChainSelection]] = None,
) -> CombinePDBResult:
    """Combine two or more PDBs and assign chains A, B, C... in input order.

    Source chains are discovered from ATOM/HETATM records in first-appearance
    order. ``chain_selections`` optionally restricts each input to some of its
    chains; it holds one entry per input file, and ``None`` (the default for
    every input) keeps all of that file's chains. Selected chains keep their
    first-appearance order, so the listed order does not reorder them.

    ATOM/HETATM/TER serials are renumbered globally, ANISOU/SIGATM/SIGUIJ
    companion serials and CONECT references are updated, LINK/REMARK/HET metadata
    follows the new chains, and MODEL wrappers are removed so all input
    coordinates form one combined structure. Records belonging to unselected
    chains are dropped along with their coordinates.
    """
    if len(input_pdbs) < 2:
        raise CombinePDBError("Choose at least two input PDB files to combine.")

    if chain_selections is None:
        selections: List[Optional[List[str]]] = [None] * len(input_pdbs)
    elif len(chain_selections) != len(input_pdbs):
        raise CombinePDBError(
            f"Chain selections were given for {len(chain_selections)} inputs, but "
            f"{len(input_pdbs)} input PDB files were chosen."
        )
    else:
        selections = [_normalize_selection(value) for value in chain_selections]

    paths = [Path(value).expanduser().resolve() for value in input_pdbs]
    output_path = Path(output_pdb).expanduser().resolve()
    if output_path in paths:
        raise CombinePDBError("The output PDB must not overwrite an input PDB.")

    all_lines: List[List[str]] = []
    chain_lists: List[List[str]] = []
    excluded_lists: List[List[str]] = []
    atom_counts: List[int] = []
    model_counts: List[int] = []
    for path, selection in zip(paths, selections):
        if not path.is_file():
            raise CombinePDBError(f"Input PDB not found: {path}")
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        chains, chain_atom_counts = _ordered_chains(lines, path)
        kept_chains = _select_chains(path, chains, selection)
        all_lines.append(lines)
        chain_lists.append(kept_chains)
        excluded_lists.append([chain for chain in chains if chain not in set(kept_chains)])
        atom_counts.append(sum(chain_atom_counts[chain] for chain in kept_chains))
        model_counts.append(sum(1 for line in lines if _record_name(line) == "MODEL "))

    multi_model_paths = [str(path) for path, count in zip(paths, model_counts) if count > 1]
    if multi_model_paths:
        raise CombinePDBError(
            "Multi-model PDB files are not supported; extract one model first: " + ", ".join(multi_model_paths)
        )

    total_chains = sum(len(chains) for chains in chain_lists)
    if total_chains > len(CHAIN_IDS):
        raise CombinePDBError(
            f"The inputs contain {total_chains} chains, but the PDB format supports only "
            f"{len(CHAIN_IDS)} uppercase alphabetic chain IDs (A-Z)."
        )

    mappings: List[InputChainMapping] = []
    chain_offset = 0
    for path, chains, excluded, atom_count in zip(paths, chain_lists, excluded_lists, atom_counts):
        chain_map = {old: CHAIN_IDS[chain_offset + index] for index, old in enumerate(chains)}
        chain_offset += len(chains)
        mappings.append(InputChainMapping(path, chain_map, atom_count, tuple(excluded)))

    source_remark_records: List[str] = []
    source_het_records: List[str] = []
    source_hetnam_records: List[str] = []
    source_link_records: List[str] = []
    output_records: List[str] = []
    conect_records: List[str] = []
    next_serial = 1
    for path, lines, mapping in zip(paths, all_lines, mappings):
        chain_map = mapping.chain_map
        excluded_chains = set(mapping.excluded_chains)
        source_chains = set(chain_map) | excluded_chains
        remap = _ChainRemapper(chain_map, excluded_chains)
        serial_map: Dict[int, int] = {}
        source_serials: Set[int] = set()
        retained_het_ids: Set[str] = set()
        pending_hetnam: List[Tuple[str, str]] = []

        # Assign the new serials first so companion and CONECT records can refer
        # to atoms that occur anywhere in this source file.
        planned_serial = next_serial
        previous_atom_chain = " "
        for line_number, line in enumerate(lines, start=1):
            record = _record_name(line)
            if record in ATOM_RECORDS:
                old_serial = _parse_serial(line, path, line_number)
                if old_serial in source_serials:
                    raise CombinePDBError(f"Duplicate atom serial {old_serial} in {path}.")
                source_serials.add(old_serial)
                previous_atom_chain = _line_chain(line)
                if previous_atom_chain in chain_map:
                    serial_map[old_serial] = planned_serial
                    planned_serial += 1
            elif record == "TER   ":
                if _ter_chain(line, source_chains, previous_atom_chain) in chain_map:
                    planned_serial += 1

        previous_atom_chain = " "
        for line_number, line in enumerate(lines, start=1):
            record = _record_name(line)
            if record == "REMARK":
                updated_remark = _update_remark_record(line, remap)
                if updated_remark is not None:
                    source_remark_records.append(updated_remark)
                continue
            if record == "HET   ":
                updated_het = _update_het_record(line, chain_map, excluded_chains)
                if updated_het is not None:
                    retained_het_ids.add(line[7:10].strip())
                    source_het_records.append(updated_het)
                continue
            if record == "HETNAM":
                pending_hetnam.append((line[11:14].strip(), line.rstrip("\r\n") + "\n"))
                continue
            if record == "LINK  ":
                updated_link = _update_link_record(line, chain_map, excluded_chains)
                if updated_link is not None:
                    source_link_records.append(updated_link)
                continue
            if record in SERIAL_RECORDS:
                if record in ATOM_RECORDS:
                    old_serial = _parse_serial(line, path, line_number)
                    old_chain = _line_chain(line)
                    previous_atom_chain = old_chain
                    if old_serial not in serial_map:
                        continue
                    new_serial = serial_map[old_serial]
                    if record == "HETATM":
                        retained_het_ids.add(line[17:20].strip())
                else:
                    old_chain = _ter_chain(line, source_chains, previous_atom_chain)
                    if old_chain not in chain_map:
                        continue
                    new_serial = next_serial
                updated = _set_serial(line, new_serial)
                if old_chain in chain_map:
                    updated = _set_chain(updated, chain_map[old_chain])
                output_records.append(updated)
                next_serial += 1
                continue

            if record in ATOM_COMPANION_RECORDS:
                old_serial = _parse_serial(line, path, line_number)
                if old_serial not in source_serials:
                    raise CombinePDBError(
                        f"{record.strip()} record on line {line_number} of {path} has no matching atom serial."
                    )
                if old_serial not in serial_map:
                    continue
                updated = _set_serial(line, serial_map[old_serial])
                old_chain = _line_chain(line)
                if old_chain in chain_map:
                    updated = _set_chain(updated, chain_map[old_chain])
                output_records.append(updated)
                continue

            if record == "CONECT":
                old_serials = _conect_serials(line, path, line_number)
                if not old_serials:
                    continue
                missing = [serial for serial in old_serials if serial not in source_serials]
                if missing:
                    raise CombinePDBError(
                        f"CONECT record on line {line_number} of {path} refers to missing atom serial(s): "
                        + ", ".join(str(value) for value in missing)
                    )
                if old_serials[0] not in serial_map:
                    continue
                partners = [serial for serial in old_serials[1:] if serial in serial_map]
                if len(old_serials) > 1 and not partners:
                    continue
                kept_serials = [old_serials[0], *partners]
                conect_records.append(
                    "CONECT" + "".join(f"{serial_map[value]:5d}" for value in kept_serials) + "\n"
                )

        source_hetnam_records.extend(
            text
            for het_id, text in pending_hetnam
            if not excluded_chains or het_id in retained_het_ids
        )

    remarks = [f"REMARK BNP_NA_COMBINE_PDB bnp_na {TOOL_VERSION}\n"]
    for index, mapping in enumerate(mappings, start=1):
        chain_summary = ", ".join(
            f"{format_chain_id(old)}->{new}" for old, new in mapping.chain_map.items()
        )
        remark = (
            f"REMARK BNP_NA_COMBINE_PDB INPUT {index} {mapping.input_pdb.name} CHAINS {chain_summary}"
        )
        if mapping.excluded_chains:
            remark += " SKIPPED " + ", ".join(
                format_chain_id(chain) for chain in mapping.excluded_chains
            )
        remarks.append(remark + "\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "".join(
            remarks
            + source_remark_records
            + source_het_records
            + source_hetnam_records
            + source_link_records
            + output_records
            + conect_records
            + ["END\n"]
        ),
        encoding="utf-8",
    )

    log_lines = [
        "=== Combine_PDB ===",
        f"Output PDB: {output_path}",
        f"Combined {len(paths)} files, {total_chains} chains, and {sum(atom_counts)} atoms.",
        f"Preserved/updated {len(source_remark_records)} REMARK and {len(source_link_records)} LINK records.",
        "Chain mappings:",
    ]
    for index, mapping in enumerate(mappings, start=1):
        chain_summary = ", ".join(
            f"{format_chain_id(old)} -> {new}" for old, new in mapping.chain_map.items()
        )
        if mapping.excluded_chains:
            chain_summary += "; skipped " + ", ".join(
                format_chain_id(chain) for chain in mapping.excluded_chains
            )
        log_lines.append(f"  Input {index}: {mapping.input_pdb} ({mapping.atom_count} atoms): {chain_summary}")

    return CombinePDBResult(
        input_pdbs=paths,
        output_pdb=output_path,
        mappings=mappings,
        atom_count=sum(atom_counts),
        chain_count=total_chains,
        remark_count=len(source_remark_records),
        link_count=len(source_link_records),
        log_text="\n".join(log_lines),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Combine PDB files and reassign source chains consecutively as A, B, C, ..."
    )
    parser.add_argument("input_pdbs", nargs="+", help="Input PDB files in the desired chain order")
    parser.add_argument("-o", "--output", required=True, help="Output combined PDB path")
    parser.add_argument(
        "-c",
        "--chains",
        action="append",
        metavar="IDS",
        help=(
            "Chains to keep from one input PDB, for example 'A B' or 'A,B'. Repeat the option "
            "once per input file, in the same order as the inputs, and use 'all' for a file "
            "that contributes every chain. Omit the option to combine all chains of every "
            "file. Use '_' to name a chain with a blank chain ID."
        ),
    )
    parser.add_argument("-v", "--version", action="version", version=f"Combine_PDB {TOOL_VERSION}")
    args = parser.parse_args(argv)
    if args.chains is not None and len(args.chains) != len(args.input_pdbs):
        parser.error(
            f"--chains was given {len(args.chains)} time(s) for {len(args.input_pdbs)} input PDB "
            "files; repeat it once per input file, using 'all' where every chain is wanted."
        )
    try:
        result = combine_pdb_files(args.input_pdbs, args.output, args.chains)
    except CombinePDBError as exc:
        parser.error(str(exc))
    print(result.log_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
