#!/usr/bin/env python3
"""
pdb_name_standard.py

Normalize nucleotide residue and atom naming in a PDB file.

This helper was changed from the previous script named
`pdb_make_dna_v3_2.py`. The public normalization behavior is preserved, but
the module name now describes the task more directly.

v3.2 history:
  - Previous script name: pdb_make_dna_v3_2.py.
  - Current module name: pdb_name_standard.py.
  - Adds public API function normalize_pdb_naming(...) so other scripts can do:
        from pdb_name_standard import normalize_pdb_naming
  - CLI main() now uses the same public API function.
  - Preserves v3.1 behavior:
      * Canonicalize standard 3-letter nucleotide residue names.
      * Convert DNA A/T/C/G to DA/DT/DC/DG.
      * Preserve RNA residues containing O2'/O2*.
      * Normalize DT methyl atom names.
      * Normalize old RNA sugar atom name O2* -> O2'.
      * Normalize phosphate atom names O1P/O2P -> OP1/OP2.
      * Optional hydrogen deletion.

Usage as command line:
    python pdb_name_standard.py input.pdb [--deleteH]

Usage as module:
    from pdb_name_standard import normalize_pdb_naming

    out_pdb = normalize_pdb_naming("input.pdb")
    out_pdb = normalize_pdb_naming("input.pdb", "output.pdb", deleteH=True)
"""

import argparse
import os
import sys
from collections import Counter
from typing import Dict, Optional, Tuple, Set, List, Any


# ---------------------------------------------------------------------------
# Residue / atom naming maps
# ---------------------------------------------------------------------------

# DNA residue-name normalization
RES_MAP = {
    "A": "DA",
    "T": "DT",
    "C": "DC",
    "G": "DG",
}

# Standard 3-letter nucleotide residue names -> 1-letter canonical names
THREE_TO_ONE_RES_MAP = {
    "ADE": "A",
    "GUA": "G",
    "CYT": "C",
    "THY": "T",
    "URI": "U",
    "URA": "U",
}

# DT methyl atom-name normalization: C5M/H5M* -> C7/H7*
DT_ATOM_MAP = {
    "C5M": "C7",
    "H5M1": "H71",
    "H5M2": "H72",
    "H5M3": "H73",
}

# Old sugar atom-name normalization, applied to all nucleotides
SUGAR_ATOM_MAP = {
    "O2*": "O2'",
}

# Backbone phosphate normalization, applied to all nucleotides
BACKBONE_OP_MAP = {
    "O1P": "OP1",
    "O2P": "OP2",
}

# Atom names that indicate a 2'-oxygen, i.e. RNA sugar
RNA_O2_ATOMS = {
    "O2'",
    "O2*",
}

ResidueID = Tuple[str, str, str]  # (chain_id, resseq, icode)


# ---------------------------------------------------------------------------
# Small PDB helpers
# ---------------------------------------------------------------------------

def _is_atom_record(line: str) -> bool:
    """Return True for ATOM/HETATM records."""
    return line.startswith("ATOM") or line.startswith("HETATM")


def _safe_slice(line: str, start: int, end: int) -> str:
    """
    Slice a line safely even if it is shorter than expected.
    Python slicing is already safe, but this helper makes intent explicit.
    """
    return line[start:end]


def _get_residue_id(line: str) -> ResidueID:
    """
    Extract residue identifier from a PDB ATOM/HETATM line.

    PDB columns:
      chainID: column 22      -> Python index 21
      resSeq : columns 23-26  -> Python slice [22:26]
      iCode  : column 27      -> Python index 26
    """
    chain_id = line[21] if len(line) > 21 else " "
    resseq = _safe_slice(line, 22, 26)
    icode = line[26] if len(line) > 26 else " "
    return chain_id, resseq, icode


def _is_hydrogen_atom(line: str) -> bool:
    """
    Return True if an ATOM/HETATM record is a hydrogen atom.

    Prefer the element column, columns 77-78 -> Python slice [76:78].
    If that is blank or absent, fall back to the first alphabetic character
    in the atom name field, columns 13-16 -> Python slice [12:16].
    This covers atom names such as 1H2'.
    """
    if not _is_atom_record(line):
        return False

    element = _safe_slice(line, 76, 78).strip().upper()
    if element:
        return element == "H"

    atom_field = _safe_slice(line, 12, 16)
    atom = atom_field.strip().upper()
    for char in atom:
        if char.isalpha():
            return char == "H"

    return False


# ---------------------------------------------------------------------------
# Detection / line-level replacement functions
# ---------------------------------------------------------------------------

def find_rna_residues(lines: List[str]) -> Set[ResidueID]:
    """
    Scan all lines once to find residues that contain a 2'-oxygen atom.

    Any residue with an O2'/O2* atom is treated as RNA and excluded from
    DNA-style residue renaming.
    """
    rna_residues: Set[ResidueID] = set()

    for line in lines:
        if not _is_atom_record(line):
            continue

        atom_field = _safe_slice(line, 12, 16)
        atom = atom_field.strip()

        if atom in RNA_O2_ATOMS:
            res_id = _get_residue_id(line)
            rna_residues.add(res_id)

    return rna_residues


def replace_resname(
    line: str,
    rna_residues: Set[ResidueID],
) -> Tuple[str, Optional[Tuple[str, str]]]:
    """
    Normalize residue name in a PDB ATOM/HETATM line.

    Returns:
        (new_line, changed_pair)

    Where changed_pair is:
        (old_resname, new_resname)

    Behavior:
      1. Standard 3-letter nucleotide residue names are canonicalized to
         1-letter residue names first.
      2. Residues whose residue-id appears in rna_residues are treated as RNA
         and therefore keep the canonical 1-letter residue name.
      3. Non-RNA A/T/C/G residues are converted to DA/DT/DC/DG.
    """
    if not _is_atom_record(line):
        return line, None

    res_id = _get_residue_id(line)

    # resName is columns 18-20, Python slice [17:20]
    resname_field = _safe_slice(line, 17, 20)
    resname = resname_field.strip()

    canonical_resname = THREE_TO_ONE_RES_MAP.get(resname, resname)

    if res_id in rna_residues:
        new_res = canonical_resname
    else:
        new_res = RES_MAP.get(canonical_resname, canonical_resname)

    if new_res != resname:
        new_res_field = new_res.rjust(3)
        line = line[:17] + new_res_field + line[20:]
        return line, (resname, new_res)

    return line, None


def replace_dt_methyl(line: str) -> Tuple[str, Optional[str]]:
    """
    If residue is DT, normalize atom names for the methyl group:

        C5M  -> C7
        H5M1 -> H71
        H5M2 -> H72
        H5M3 -> H73

    Returns:
        (new_line, old_atomname_if_changed)
    """
    if not _is_atom_record(line):
        return line, None

    resname = _safe_slice(line, 17, 20).strip()
    if resname != "DT":
        return line, None

    atom_field = _safe_slice(line, 12, 16)
    atom = atom_field.strip()

    new_atom = DT_ATOM_MAP.get(atom)
    if new_atom is None:
        return line, None

    atom_formatted = new_atom.rjust(4)
    line = line[:12] + atom_formatted + line[16:]
    return line, atom


def replace_sugar_old_names(line: str) -> Tuple[str, Optional[str]]:
    """
    Normalize old sugar atom names, currently:

        O2* -> O2'

    Applied to all ATOM/HETATM records after RNA detection. Returns:
        (new_line, old_atomname_if_changed)
    """
    if not _is_atom_record(line):
        return line, None

    atom_field = _safe_slice(line, 12, 16)
    atom = atom_field.strip()

    new_atom = SUGAR_ATOM_MAP.get(atom)
    if new_atom is None:
        return line, None

    atom_formatted = new_atom.rjust(4)
    line = line[:12] + atom_formatted + line[16:]
    return line, atom


def replace_backbone_op(line: str) -> Tuple[str, Optional[str]]:
    """
    Normalize backbone phosphate atom names:

        O1P -> OP1
        O2P -> OP2

    Applied to all ATOM/HETATM records.

    Returns:
        (new_line, old_atomname_if_changed)
    """
    if not _is_atom_record(line):
        return line, None

    atom_field = _safe_slice(line, 12, 16)
    atom = atom_field.strip()

    new_atom = BACKBONE_OP_MAP.get(atom)
    if new_atom is None:
        return line, None

    atom_formatted = new_atom.rjust(4)
    line = line[:12] + atom_formatted + line[16:]
    return line, atom


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def make_outname(inp: str) -> str:
    """
    Make default output filename.

    If input is xxx.pdb, output is xxx_Dout.pdb.
    Otherwise output is xxx_Dout.
    """
    root, ext = os.path.splitext(inp)
    if ext.lower() == ".pdb":
        return root + "_Dout.pdb"
    return inp + "_Dout"


def normalize_pdb_naming(
    inp: str,
    outp: Optional[str] = None,
    deleteH: bool = False,
    return_report: bool = False,
) -> str:
    """
    Normalize nucleotide residue and atom naming in a PDB file.

    Parameters
    ----------
    inp:
        Input PDB path.

    outp:
        Output PDB path. If None, uses make_outname(inp).

    deleteH:
        If True, delete hydrogen ATOM/HETATM records from output.

    return_report:
        If False, return only output path as str.
        If True, return a dict-like report is NOT returned here for backward
        compatibility. Use normalize_pdb_naming_with_report(...) instead.

    Returns
    -------
    str
        Output PDB path.

    Notes
    -----
    The return_report argument is accepted for forward compatibility but this
    function intentionally returns str. For detailed counts, call
    normalize_pdb_naming_with_report(...).
    """
    report = normalize_pdb_naming_with_report(
        inp=inp,
        outp=outp,
        deleteH=deleteH,
    )
    return str(report["output"])


def normalize_pdb_naming_with_report(
    inp: str,
    outp: Optional[str] = None,
    deleteH: bool = False,
) -> Dict[str, Any]:
    """
    Normalize nucleotide residue and atom naming in a PDB file and return counts.

    Returns a dictionary with:
        input
        output
        total_lines
        deleted_hydrogen_atoms
        residue_changes
        atom_changes
        rna_residue_count
    """
    if outp is None:
        outp = make_outname(inp)

    try:
        with open(inp, "r", encoding="utf-8", errors="ignore") as fin:
            lines = fin.readlines()
    except OSError as e:
        raise OSError(f"Error reading {inp}: {e}") from e

    rna_residues = find_rna_residues(lines)

    res_counts: Counter = Counter()
    atom_counts: Counter = Counter()
    deleted_h_count = 0
    total = 0

    try:
        with open(outp, "w", encoding="utf-8") as fout:
            for line in lines:
                total += 1

                if deleteH and _is_hydrogen_atom(line):
                    deleted_h_count += 1
                    continue

                # Step 1: residue renaming / canonicalization
                line, changed_res = replace_resname(line, rna_residues)
                if changed_res is not None:
                    res_counts[changed_res] += 1

                # Step 2: DT methyl atom renaming
                line, changed_atom = replace_dt_methyl(line)
                if changed_atom is not None:
                    atom_counts[changed_atom] += 1

                # Step 3: Old sugar atom-name renaming, e.g. O2* -> O2'
                line, changed_atom2 = replace_sugar_old_names(line)
                if changed_atom2 is not None:
                    atom_counts[changed_atom2] += 1

                # Step 4: Backbone phosphate O1P/O2P renaming
                line, changed_atom3 = replace_backbone_op(line)
                if changed_atom3 is not None:
                    atom_counts[changed_atom3] += 1

                fout.write(line)
    except OSError as e:
        raise OSError(f"Error writing {outp}: {e}") from e

    return {
        "input": inp,
        "output": outp,
        "total_lines": total,
        "deleted_hydrogen_atoms": deleted_h_count,
        "residue_changes": dict(res_counts),
        "atom_changes": dict(atom_counts),
        "rna_residue_count": len(rna_residues),
    }


def format_report(report: Dict[str, Any], deleteH: bool = False) -> str:
    """
    Format a normalization report as human-readable text.
    """
    lines: List[str] = []

    lines.append(f"Processed {report.get('total_lines', 0)} lines.")

    if deleteH:
        lines.append(f"  Deleted hydrogen atoms : {report.get('deleted_hydrogen_atoms', 0)}")

    res_counts = Counter(report.get("residue_changes", {}))
    atom_counts = Counter(report.get("atom_changes", {}))

    if res_counts:
        report_order = [
            ("A", "DA"),
            ("T", "DT"),
            ("C", "DC"),
            ("G", "DG"),
            ("ADE", "A"),
            ("THY", "T"),
            ("CYT", "C"),
            ("GUA", "G"),
            ("URI", "U"),
            ("URA", "U"),
            ("ADE", "DA"),
            ("THY", "DT"),
            ("CYT", "DC"),
            ("GUA", "DG"),
        ]

        reported = set()
        for old_res, new_res in report_order:
            key = (old_res, new_res)
            if res_counts.get(key):
                lines.append(f"  {old_res} -> {new_res} : {res_counts[key]}")
                reported.add(key)

        for old_res, new_res in sorted(res_counts):
            key = (old_res, new_res)
            if key not in reported:
                lines.append(f"  {old_res} -> {new_res} : {res_counts[key]}")

    if atom_counts:
        for k in ("C5M", "H5M1", "H5M2", "H5M3"):
            if atom_counts.get(k):
                lines.append(f"  DT:{k} -> {DT_ATOM_MAP[k]} : {atom_counts[k]}")

        for k in ("O2*",):
            if atom_counts.get(k):
                lines.append(f"  {k} -> {SUGAR_ATOM_MAP[k]} : {atom_counts[k]}")

        for k in ("O1P", "O2P"):
            if atom_counts.get(k):
                lines.append(f"  {k} -> {BACKBONE_OP_MAP[k]} : {atom_counts[k]}")

    lines.append(f"Output written to: {report.get('output')}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize nucleotide residue and atom naming in a PDB file.",
    )
    parser.add_argument(
        "input_pdb",
        help="Input PDB file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_pdb",
        default=None,
        help="Output PDB file. Default: <input>_Dout.pdb or <input>_Dout.",
    )
    parser.add_argument(
        "--deleteH",
        action="store_true",
        help="Delete hydrogen ATOM/HETATM records from the output.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    try:
        report = normalize_pdb_naming_with_report(
            inp=args.input_pdb,
            outp=args.output_pdb,
            deleteH=args.deleteH,
        )
    except OSError as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

    print(format_report(report, deleteH=args.deleteH))


if __name__ == "__main__":
    main()
