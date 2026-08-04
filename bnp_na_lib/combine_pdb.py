#!/usr/bin/env python3
"""Combine PDB files while assigning consecutive, unique chain IDs."""
from __future__ import annotations

import argparse
import re
import string
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union


TOOL_VERSION = "V13.13"
CHAIN_IDS = string.ascii_uppercase
ATOM_RECORDS = ("ATOM  ", "HETATM")
SERIAL_RECORDS = ATOM_RECORDS + ("TER   ",)
ATOM_COMPANION_RECORDS = ("ANISOU", "SIGATM", "SIGUIJ")


class CombinePDBError(Exception):
    """Raised when PDB inputs cannot be combined safely."""


@dataclass(frozen=True)
class InputChainMapping:
    input_pdb: Path
    chain_map: Dict[str, str]
    atom_count: int


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


def _chain_display(chain_id: str) -> str:
    return chain_id if chain_id.strip() else "(blank)"


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


def _remark_chain_id(chain_id: str, chain_map: Dict[str, str]) -> str:
    source_chain = " " if chain_id == "_" else chain_id
    return chain_map.get(source_chain, chain_id)


def _replace_remark_field(line: str, field: str, replacement: Callable[[str], str]) -> str:
    pattern = re.compile(rf"(?<!\S)({re.escape(field)}=)(\S+)")

    def replace(match: re.Match[str]) -> str:
        return match.group(1) + replacement(match.group(2))

    return pattern.sub(replace, line)


def _replace_residue_labels(text: str, chain_map: Dict[str, str]) -> str:
    """Update parse-friendly ``A:12:DA`` and DSSR-style ``A.DA12`` labels."""
    colon_pattern = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z0-9_])(?=:[+-]?\d+(?::[A-Za-z0-9]+)?)")
    dot_pattern = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z0-9_])(?=\.[A-Za-z0-9]{1,3}[+-]?\d+)")

    def replace(match: re.Match[str]) -> str:
        return _remark_chain_id(match.group(1), chain_map)

    return dot_pattern.sub(replace, colon_pattern.sub(replace, text))


def _update_remark_record(line: str, chain_map: Dict[str, str]) -> str:
    """Update chain references in known structured and common REMARK formats."""
    if line.startswith("REMARK 950 RE_SCRIPT COMMAND"):
        return line if line.endswith(("\n", "\r")) else line + "\n"

    if line.startswith("REMARK 950 RE_SCRIPT CHAIN_RANGE"):
        updated = _replace_remark_field(line, "chain", lambda value: _remark_chain_id(value, chain_map))
        for field in ("start", "end"):
            updated = _replace_remark_field(updated, field, lambda value: _replace_residue_labels(value, chain_map))
    elif line.startswith("REMARK 950 RE_SCRIPT CHAIN_RESIDUES"):
        updated = _replace_remark_field(line, "chain", lambda value: _remark_chain_id(value, chain_map))
        updated = _replace_remark_field(
            updated, "residues", lambda value: _replace_residue_labels(value, chain_map)
        )
    elif line.startswith("REMARK 950 RE_SCRIPT JUNCTION"):
        updated = line
        for field in ("residues", "core", "excluded_nick_ends"):
            updated = _replace_remark_field(
                updated, field, lambda value: _replace_residue_labels(value, chain_map)
            )
    elif line.startswith("REMARK 950 RE_SCRIPT SPECIAL"):
        updated = _replace_remark_field(line, "chain", lambda value: _remark_chain_id(value, chain_map))
        updated = _replace_remark_field(
            updated, "residue", lambda value: _replace_residue_labels(value, chain_map)
        )
    elif line.startswith("REMARK 950 RE_SCRIPT"):
        # Unrecognized RE_SCRIPT records are provenance. Preserve them exactly
        # rather than risk changing original/source residue identities.
        updated = line
    else:
        updated = _replace_remark_field(line, "chain", lambda value: _remark_chain_id(value, chain_map))
        updated = re.sub(
            r"(?<!\S)(CHAIN\s+)([A-Za-z0-9_])(?=\s|$)",
            lambda match: match.group(1) + _remark_chain_id(match.group(2), chain_map),
            updated,
        )
        updated = _replace_residue_labels(updated, chain_map)

    return updated.rstrip("\r\n") + "\n"


def _update_link_record(line: str, chain_map: Dict[str, str]) -> str:
    """Update both fixed-column chain IDs in a PDB LINK record."""
    updated = line
    for index in (21, 51):
        old_chain = updated[index] if len(updated.rstrip("\r\n")) > index else " "
        if old_chain in chain_map:
            updated = _replace_columns(updated, index, index + 1, chain_map[old_chain])
    return updated.rstrip("\r\n") + "\n"


def _update_het_record(line: str, chain_map: Dict[str, str]) -> str:
    """Update the fixed-column chain ID in a PDB HET residue record."""
    old_chain = line[12] if len(line.rstrip("\r\n")) > 12 else " "
    if old_chain in chain_map:
        return _replace_columns(line, 12, 13, chain_map[old_chain])
    return line.rstrip("\r\n") + "\n"


def _ordered_chains(lines: Sequence[str], path: Path) -> Tuple[List[str], int]:
    chains: List[str] = []
    atom_count = 0
    for line in lines:
        if _record_name(line) not in ATOM_RECORDS:
            continue
        atom_count += 1
        chain_id = line[21] if len(line) > 21 else " "
        if chain_id not in chains:
            chains.append(chain_id)
    if not atom_count:
        raise CombinePDBError(f"No ATOM/HETATM records were found in {path}.")
    return chains, atom_count


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


def default_combine_pdb_output_path(output_dir: Union[Path, str]) -> Path:
    """Return the GUI's default output path for the combine_PDB tool."""
    return Path(output_dir).expanduser() / "combine_PDB_out.pdb"


def combine_pdb_files(
    input_pdbs: Sequence[Union[Path, str]],
    output_pdb: Union[Path, str],
) -> CombinePDBResult:
    """Combine two or more PDBs and assign chains A, B, C... in input order.

    Source chains are discovered from ATOM/HETATM records in first-appearance
    order. ATOM/HETATM/TER serials are renumbered globally, ANISOU/SIGATM/SIGUIJ
    companion serials and CONECT references are updated, LINK/REMARK/HET metadata
    follows the new chains, and MODEL wrappers are removed so all input
    coordinates form one combined structure.
    """
    if len(input_pdbs) < 2:
        raise CombinePDBError("Choose at least two input PDB files to combine.")

    paths = [Path(value).expanduser().resolve() for value in input_pdbs]
    output_path = Path(output_pdb).expanduser().resolve()
    if output_path in paths:
        raise CombinePDBError("The output PDB must not overwrite an input PDB.")

    all_lines: List[List[str]] = []
    chain_lists: List[List[str]] = []
    atom_counts: List[int] = []
    model_counts: List[int] = []
    for path in paths:
        if not path.is_file():
            raise CombinePDBError(f"Input PDB not found: {path}")
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        chains, atom_count = _ordered_chains(lines, path)
        all_lines.append(lines)
        chain_lists.append(chains)
        atom_counts.append(atom_count)
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
    for path, chains, atom_count in zip(paths, chain_lists, atom_counts):
        chain_map = {old: CHAIN_IDS[chain_offset + index] for index, old in enumerate(chains)}
        chain_offset += len(chains)
        mappings.append(InputChainMapping(path, chain_map, atom_count))

    source_remark_records: List[str] = []
    source_het_records: List[str] = []
    source_hetnam_records: List[str] = []
    source_link_records: List[str] = []
    output_records: List[str] = []
    conect_records: List[str] = []
    next_serial = 1
    for path, lines, mapping in zip(paths, all_lines, mappings):
        serial_map: Dict[int, int] = {}

        # Assign the new serials first so companion and CONECT records can refer
        # to atoms that occur anywhere in this source file.
        planned_serial = next_serial
        for line_number, line in enumerate(lines, start=1):
            record = _record_name(line)
            if record in ATOM_RECORDS:
                old_serial = _parse_serial(line, path, line_number)
                if old_serial in serial_map:
                    raise CombinePDBError(f"Duplicate atom serial {old_serial} in {path}.")
                serial_map[old_serial] = planned_serial
                planned_serial += 1
            elif record == "TER   ":
                planned_serial += 1

        for line_number, line in enumerate(lines, start=1):
            record = _record_name(line)
            if record == "REMARK":
                source_remark_records.append(_update_remark_record(line, mapping.chain_map))
                continue
            if record == "HET   ":
                source_het_records.append(_update_het_record(line, mapping.chain_map))
                continue
            if record == "HETNAM":
                source_hetnam_records.append(line.rstrip("\r\n") + "\n")
                continue
            if record == "LINK  ":
                source_link_records.append(_update_link_record(line, mapping.chain_map))
                continue
            if record in SERIAL_RECORDS:
                if record in ATOM_RECORDS:
                    old_serial = _parse_serial(line, path, line_number)
                    new_serial = serial_map[old_serial]
                else:
                    new_serial = next_serial
                updated = _set_serial(line, new_serial)
                old_chain = line[21] if len(line) > 21 else " "
                if old_chain in mapping.chain_map:
                    updated = _set_chain(updated, mapping.chain_map[old_chain])
                output_records.append(updated)
                next_serial += 1
                continue

            if record in ATOM_COMPANION_RECORDS:
                old_serial = _parse_serial(line, path, line_number)
                if old_serial not in serial_map:
                    raise CombinePDBError(
                        f"{record.strip()} record on line {line_number} of {path} has no matching atom serial."
                    )
                updated = _set_serial(line, serial_map[old_serial])
                old_chain = line[21] if len(line) > 21 else " "
                if old_chain in mapping.chain_map:
                    updated = _set_chain(updated, mapping.chain_map[old_chain])
                output_records.append(updated)
                continue

            if record == "CONECT":
                old_serials = _conect_serials(line, path, line_number)
                if not old_serials:
                    continue
                missing = [serial for serial in old_serials if serial not in serial_map]
                if missing:
                    raise CombinePDBError(
                        f"CONECT record on line {line_number} of {path} refers to missing atom serial(s): "
                        + ", ".join(str(value) for value in missing)
                    )
                conect_records.append("CONECT" + "".join(f"{serial_map[value]:5d}" for value in old_serials) + "\n")

    remarks = [f"REMARK BNP_NA_COMBINE_PDB bnp_na {TOOL_VERSION}\n"]
    for index, mapping in enumerate(mappings, start=1):
        chain_summary = ", ".join(
            f"{_chain_display(old)}->{new}" for old, new in mapping.chain_map.items()
        )
        remarks.append(
            f"REMARK BNP_NA_COMBINE_PDB INPUT {index} {mapping.input_pdb.name} CHAINS {chain_summary}\n"
        )

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
        "=== combine_PDB ===",
        f"Output PDB: {output_path}",
        f"Combined {len(paths)} files, {total_chains} chains, and {sum(atom_counts)} atoms.",
        f"Preserved/updated {len(source_remark_records)} REMARK and {len(source_link_records)} LINK records.",
        "Chain mappings:",
    ]
    for index, mapping in enumerate(mappings, start=1):
        chain_summary = ", ".join(
            f"{_chain_display(old)} -> {new}" for old, new in mapping.chain_map.items()
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
    parser.add_argument("-v", "--version", action="version", version=f"combine_PDB {TOOL_VERSION}")
    args = parser.parse_args(argv)
    try:
        result = combine_pdb_files(args.input_pdbs, args.output)
    except CombinePDBError as exc:
        parser.error(str(exc))
    print(result.log_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
