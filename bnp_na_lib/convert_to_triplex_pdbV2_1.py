#!/usr/bin/env python3
"""
convert_to_triplex_pdbV2_1.py

Bundled with bnp_na V13.7 from the standalone convert_to_triplex_pdbV2_1.py script.

Convert a DNA duplex PDB into a DNA triplex by aligning a reference base-triple
onto a residue range of strand I (the purine strand) and building strand III
(the Hoogsteen strand) residue-by-residue.

Conventions used here
---------------------
Strand I   : purine strand in the duplex (X in Z·X-Y)
Strand II  : pyrimidine strand in the duplex (Y in Z·X-Y)
Strand III : Hoogsteen / reverse-Hoogsteen strand (Z in Z·X-Y)

Supported triplets
------------------
antiparallel : G·G-C  (embedded template)
parallel     : T·A-T  (embedded template)

Notes
-----
* Chain IDs A, B, C in the templates correspond to X, Y, Z respectively.
* The third strand is written in 5' -> 3' order in the final PDB.
* If no command-line arguments are given, or --gui is supplied, the script
  starts a Tk GUI.
* If -o/--out is omitted, the output filename is the input basename with
  _2TH added before the extension.

Author: updated from the original script supplied by Di Liu.
"""

from __future__ import annotations

import argparse
import copy
import io
import os
import string
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from edit_pdb_atom import file2rec, rec2file, pdb_ter_record

# ---------------------------------------------------------------------------
# Embedded templates supplied by the user.
# Chain IDs A, B, C correspond to X, Y, Z in Z·X-Y.
# ---------------------------------------------------------------------------

TEMPLATE_G_G_C = """\
ATOM      1  C1'  DC B   9       1.016  -7.216  34.261  1.00  0.00           C
ATOM      2  C2   DC B   9       1.193  -4.765  34.011  1.00  0.00           C
ATOM      3  C2'  DC B   9       0.505  -7.965  35.472  1.00  0.00           C
ATOM      4  C3'  DC B   9       0.805  -9.401  35.091  1.00  0.00           C
ATOM      5  C4   DC B   9      -0.711  -3.456  33.799  1.00  0.00           C
ATOM      6  C4'  DC B   9       0.607  -9.417  33.573  1.00  0.00           C
ATOM      7  C5   DC B   9      -1.559  -4.600  33.828  1.00  0.00           C
ATOM      8  C5'  DC B   9      -0.706  -9.991  33.093  1.00  0.00           C
ATOM      9  C6   DC B   9      -0.972  -5.795  33.959  1.00  0.00           C
ATOM     10  N1   DC B   9       0.386  -5.905  34.065  1.00  0.00           N
ATOM     11  N3   DC B   9       0.612  -3.557  33.881  1.00  0.00           N
ATOM     12  N4   DC B   9      -1.227  -2.238  33.692  1.00  0.00           N
ATOM     13  OP1  DC B   9      -4.189  -9.032  34.354  1.00  0.00           O
ATOM     14  O2   DC B   9       2.428  -4.888  34.081  1.00  0.00           O
ATOM     15  OP2  DC B   9      -3.343 -11.219  33.309  1.00  0.00           O
ATOM     16  O3'  DC B   9       2.157  -9.724  35.429  1.00  0.00           O
ATOM     17  O4'  DC B   9       0.710  -8.034  33.145  1.00  0.00           O
ATOM     18  O5'  DC B   9      -1.806  -9.282  33.670  1.00  0.00           O
ATOM     19  P    DC B   9      -3.308  -9.736  33.386  1.00  0.00           P
TER      20       DC B   9
ATOM     21  C1'  DG A  22       6.550   1.834  33.756  1.00  0.00           C
ATOM     22  C2   DG A  22       4.208  -1.793  33.667  1.00  0.00           C
ATOM     23  C2'  DG A  22       6.858   2.696  32.549  1.00  0.00           C
ATOM     24  C3'  DG A  22       8.248   3.206  32.877  1.00  0.00           C
ATOM     25  C4   DG A  22       4.488   0.390  33.916  1.00  0.00           C
ATOM     26  C4'  DG A  22       8.240   3.323  34.403  1.00  0.00           C
ATOM     27  C5   DG A  22       3.151   0.660  34.121  1.00  0.00           C
ATOM     28  C5'  DG A  22       8.201   4.719  34.976  1.00  0.00           C
ATOM     29  C6   DG A  22       2.225  -0.421  34.085  1.00  0.00           C
ATOM     30  C8   DG A  22       4.151   2.535  34.268  1.00  0.00           C
ATOM     31  N1   DG A  22       2.864  -1.631  33.852  1.00  0.00           N
ATOM     32  N2   DG A  22       4.618  -3.050  33.460  1.00  0.00           N
ATOM     33  N3   DG A  22       5.081  -0.802  33.684  1.00  0.00           N
ATOM     34  N7   DG A  22       2.951   2.017  34.341  1.00  0.00           N
ATOM     35  N9   DG A  22       5.126   1.606  33.996  1.00  0.00           N
ATOM     36  OP1  DG A  22       6.043   7.648  33.781  1.00  0.00           O
ATOM     37  OP2  DG A  22       8.552   7.603  34.325  1.00  0.00           O
ATOM     38  O3'  DG A  22       9.248   2.280  32.453  1.00  0.00           O
ATOM     39  O4'  DG A  22       7.110   2.528  34.859  1.00  0.00           O
ATOM     40  O5'  DG A  22       7.195   5.505  34.347  1.00  0.00           O
ATOM     41  O6   DG A  22       0.992  -0.397  34.222  1.00  0.00           O
ATOM     42  P    DG A  22       7.179   7.084  34.554  1.00  0.00           P
TER      43       DG A  22
ATOM     44  C1'  DG C  39      -4.154   5.394  33.675  1.00  0.00           C
ATOM     45  C2   DG C  39      -1.235   2.212  34.328  1.00  0.00           C
ATOM     46  C2'  DG C  39      -4.820   6.135  34.822  1.00  0.00           C
ATOM     47  C3'  DG C  39      -6.061   6.719  34.169  1.00  0.00           C
ATOM     48  C4   DG C  39      -1.875   4.298  33.999  1.00  0.00           C
ATOM     49  C4'  DG C  39      -5.615   6.978  32.728  1.00  0.00           C
ATOM     50  C5   DG C  39      -0.595   4.815  34.010  1.00  0.00           C
ATOM     51  C5'  DG C  39      -5.194   8.389  32.388  1.00  0.00           C
ATOM     52  C6   DG C  39       0.499   3.936  34.194  1.00  0.00           C
ATOM     53  C8   DG C  39      -1.882   6.473  33.702  1.00  0.00           C
ATOM     54  N1   DG C  39       0.077   2.622  34.348  1.00  0.00           N
ATOM     55  N2   DG C  39      -1.449   0.918  34.513  1.00  0.00           N
ATOM     56  N3   DG C  39      -2.264   3.018  34.147  1.00  0.00           N
ATOM     57  N7   DG C  39      -0.617   6.189  33.826  1.00  0.00           N
ATOM     58  N9   DG C  39      -2.699   5.375  33.787  1.00  0.00           N
ATOM     59  OP1  DG C  39      -1.955  10.156  33.399  1.00  0.00           O
ATOM     60  OP2  DG C  39      -3.997  11.028  32.112  1.00  0.00           O
ATOM     61  O3'  DG C  39      -7.194   5.830  34.136  1.00  0.00           O
ATOM     62  O4'  DG C  39      -4.454   6.146  32.515  1.00  0.00           O
ATOM     63  O5'  DG C  39      -3.907   8.677  32.953  1.00  0.00           O
ATOM     64  O6   DG C  39       1.703   4.203  34.225  1.00  0.00           O
ATOM     65  P    DG C  39      -3.058   9.922  32.431  1.00  0.00           P
TER      66       DG C  39
END
"""

TEMPLATE_T_A_T = """\
ATOM      1  C1'  DT B   3       1.565   6.774  14.023  1.00  0.00           C
ATOM      2  C2   DT B   3       0.417   4.635  13.814  1.00  0.00           C
ATOM      3  C2'  DT B   3       2.266   7.330  15.243  1.00  0.00           C
ATOM      4  C3'  DT B   3       2.558   8.750  14.803  1.00  0.00           C
ATOM      5  C4   DT B   3       1.679   2.530  13.733  1.00  0.00           C
ATOM      6  C4'  DT B   3       2.823   8.613  13.300  1.00  0.00           C
ATOM      7  C5   DT B   3       2.906   3.292  13.827  1.00  0.00           C
ATOM      8  C5'  DT B   3       4.273   8.647  12.878  1.00  0.00           C
ATOM      9  C6   DT B   3       2.818   4.629  13.897  1.00  0.00           C
ATOM     10  C7   DT B   3       4.208   2.558  13.915  1.00  0.00           C
ATOM     11  N1   DT B   3       1.615   5.304  13.902  1.00  0.00           N
ATOM     12  N3   DT B   3       0.524   3.277  13.733  1.00  0.00           N
ATOM     13  OP1  DT B   3       7.052   6.280  14.048  1.00  0.00           O
ATOM     14  O2   DT B   3      -0.663   5.200  13.809  1.00  0.00           O
ATOM     15  OP2  DT B   3       7.231   8.682  13.157  1.00  0.00           O
ATOM     16  O3'  DT B   3       1.417   9.581  15.039  1.00  0.00           O
ATOM     17  O4   DT B   3       1.611   1.308  13.646  1.00  0.00           O
ATOM     18  O4'  DT B   3       2.228   7.348  12.909  1.00  0.00           O
ATOM     19  O5'  DT B   3       5.008   7.597  13.510  1.00  0.00           O
ATOM     20  P    DT B   3       6.547   7.363  13.165  1.00  0.00           P
TER      21       DT B   3
ATOM     22  C1'  DA A  10      -7.068   0.476  13.406  1.00  0.00           C
ATOM     23  C2   DA A  10      -3.369   2.741  13.388  1.00  0.00           C
ATOM     24  C2'  DA A  10      -7.733  -0.107  12.177  1.00  0.00           C
ATOM     25  C3'  DA A  10      -9.199  -0.060  12.562  1.00  0.00           C
ATOM     26  C4   DA A  10      -4.589   0.917  13.529  1.00  0.00           C
ATOM     27  C4'  DA A  10      -9.180  -0.291  14.076  1.00  0.00           C
ATOM     28  C5   DA A  10      -3.474   0.112  13.674  1.00  0.00           C
ATOM     29  C5'  DA A  10      -9.602  -1.661  14.552  1.00  0.00           C
ATOM     30  C6   DA A  10      -2.218   0.754  13.659  1.00  0.00           C
ATOM     31  C8   DA A  10      -5.157  -1.186  13.757  1.00  0.00           C
ATOM     32  N1   DA A  10      -2.204   2.096  13.515  1.00  0.00           N
ATOM     33  N3   DA A  10      -4.606   2.252  13.381  1.00  0.00           N
ATOM     34  N6   DA A  10      -1.049   0.119  13.777  1.00  0.00           N
ATOM     35  N7   DA A  10      -3.846  -1.218  13.813  1.00  0.00           N
ATOM     36  N9   DA A  10      -5.672   0.073  13.577  1.00  0.00           N
ATOM     37  OP1  DA A  10      -8.454  -5.040  13.149  1.00  0.00           O
ATOM     38  OP2  DA A  10     -10.764  -4.346  14.025  1.00  0.00           O
ATOM     39  O3'  DA A  10      -9.761   1.216  12.244  1.00  0.00           O
ATOM     40  O4'  DA A  10      -7.822  -0.011  14.503  1.00  0.00           O
ATOM     41  O5'  DA A  10      -8.912  -2.682  13.834  1.00  0.00           O
ATOM     42  P    DA A  10      -9.284  -4.215  14.061  1.00  0.00           P
TER      43       DA A  10
ATOM     44  C1'  DT C  16      -2.958  -6.524  13.253  1.00  0.00           C
ATOM     45  C2   DT C  16      -2.567  -4.143  13.601  1.00  0.00           C
ATOM     46  C2'  DT C  16      -2.714  -7.309  11.983  1.00  0.00           C
ATOM     47  C3'  DT C  16      -3.331  -8.649  12.327  1.00  0.00           C
ATOM     48  C4   DT C  16      -0.281  -3.267  13.812  1.00  0.00           C
ATOM     49  C4'  DT C  16      -3.063  -8.790  13.830  1.00  0.00           C
ATOM     50  C5   DT C  16       0.223  -4.614  13.660  1.00  0.00           C
ATOM     51  C5'  DT C  16      -1.921  -9.702  14.219  1.00  0.00           C
ATOM     52  C6   DT C  16      -0.664  -5.607  13.505  1.00  0.00           C
ATOM     53  C7   DT C  16       1.703  -4.830  13.605  1.00  0.00           C
ATOM     54  N1   DT C  16      -2.028  -5.400  13.463  1.00  0.00           N
ATOM     55  N3   DT C  16      -1.651  -3.143  13.771  1.00  0.00           N
ATOM     56  OP1  DT C  16       1.689  -9.504  12.940  1.00  0.00           O
ATOM     57  O2   DT C  16      -3.765  -3.927  13.571  1.00  0.00           O
ATOM     58  OP2  DT C  16       0.327 -11.557  13.661  1.00  0.00           O
ATOM     59  O3'  DT C  16      -4.741  -8.642  12.064  1.00  0.00           O
ATOM     60  O4   DT C  16       0.414  -2.272  13.989  1.00  0.00           O
ATOM     61  O4'  DT C  16      -2.795  -7.448  14.313  1.00  0.00           O
ATOM     62  O5'  DT C  16      -0.701  -9.279  13.603  1.00  0.00           O
ATOM     63  P    DT C  16       0.645 -10.113  13.804  1.00  0.00           P
TER      64      THY C  16
END
"""

# ---------------------------------------------------------------------------
# Chemistry helpers
# ---------------------------------------------------------------------------

BASE_ALIASES = {
    "A": "A",
    "DA": "A",
    "ADE": "A",
    "G": "G",
    "DG": "G",
    "GUA": "G",
    "C": "C",
    "DC": "C",
    "CYT": "C",
    "T": "T",
    "DT": "T",
    "THY": "T",
    "U": "U",
    "DU": "U",
    "URA": "U",
}

PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T", "U"}

ALIGN_ATOMS = {
    "A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
    "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
    "C": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "T": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
}

MODE_ALIASES = {
    "antiparallel": "antiparallel",
    "anti": "antiparallel",
    "ap": "antiparallel",
    "parallel": "parallel",
    "para": "parallel",
    "p": "parallel",
}


@dataclass
class Residue:
    chain_id: str
    res_seq: int
    res_name: str
    index_in_chain: int
    records: List[object] = field(default_factory=list)

    @property
    def base(self) -> str:
        return normalize_base(self.res_name)

    @property
    def label(self) -> str:
        return f"{self.chain_id}:{self.res_seq}({self.res_name.strip()})"

    def atom_map(self) -> Dict[str, object]:
        return {rec.name.strip(): rec for rec in self.records if rec.recordName in ("ATOM", "HETATM")}


@dataclass
class Motif:
    mode: str
    label: str
    x_base: str
    y_base: str
    z_base: str
    template_text: str
    x_residue: Residue = field(init=False)
    y_residue: Residue = field(init=False)
    z_residue: Residue = field(init=False)
    ref_pair_coords: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        records = parse_pdb_text(self.template_text)
        chains = build_chain_residues(records)
        try:
            self.x_residue = chains["A"][0]
            self.y_residue = chains["B"][0]
            self.z_residue = chains["C"][0]
        except KeyError as exc:
            raise ValueError(
                f"Template for {self.label} must contain chains A/B/C as X/Y/Z."
            ) from exc
        self.ref_pair_coords = concat_base_coords(
            self.x_residue, self.x_base, self.y_residue, self.y_base
        )





# ---------------------------------------------------------------------------
# Low-level utilities
# ---------------------------------------------------------------------------


def normalize_mode(mode: str) -> str:
    key = mode.strip().lower()
    if key not in MODE_ALIASES:
        choices = ", ".join(sorted(set(MODE_ALIASES.values())))
        raise ValueError(f"Unsupported mode '{mode}'. Choose one of: {choices}.")
    return MODE_ALIASES[key]


def normalize_base(res_name: str) -> str:
    key = res_name.strip().upper()
    if key not in BASE_ALIASES:
        return key
    return BASE_ALIASES[key]


def parse_residue_range(text: str) -> Tuple[int, int]:
    cleaned = text.strip()
    for sep in (":", "-", "..", ","):
        if sep in cleaned:
            left, right = cleaned.split(sep, 1)
            return int(left.strip()), int(right.strip())
    value = int(cleaned)
    return value, value


def parse_pdb_text(text: str) -> List[object]:
    records: List[object] = []
    lines = text.splitlines(True)
    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    file2rec(lines, records)
    return records


def read_pdb_file(path: str) -> List[object]:
    records: List[object] = []
    with open(path, "r", encoding="utf-8") as handle:
        file2rec(handle, records)
    return records

def default_output_path(input_path: str) -> str:
    """Return the default output path: input basename + _2TH before the extension."""
    path = Path(input_path)
    if path.suffix:
        return str(path.with_name(f"{path.stem}_2TH{path.suffix}"))
    if path.name:
        return str(path.with_name(f"{path.name}_2TH.pdb"))
    return "triplex_2TH.pdb"


def display_chain_id(chain_id: str) -> str:
    return chain_id if chain_id.strip() else "<blank>"


def residue_sequence(residues: Sequence[Residue]) -> str:
    if not residues:
        return ""
    return "".join(res.base if res.base in {"A", "C", "G", "T", "U"} else "?" for res in residues)


def format_residue_number_list(residues: Sequence[Residue], max_items: int = 28) -> str:
    numbers = [str(res.res_seq) for res in residues]
    if len(numbers) <= max_items:
        return ", ".join(numbers)
    head = ", ".join(numbers[: max_items // 2])
    tail = ", ".join(numbers[-(max_items // 2) :])
    return f"{head}, ..., {tail}"


def format_chain_sequence_info(pdb_path: str) -> str:
    """Summarize chains and nucleotide sequences in PDB residue order."""
    records = read_pdb_file(pdb_path)
    chains = build_chain_residues(records)
    lines = [
        f"Input duplex PDB: {pdb_path}",
        "",
        "Detected strands/chains and sequences:",
        "  Sequence is reported in the residue order found in the PDB file; for a standard DNA strand this should be 5' -> 3'.",
    ]

    if not chains:
        lines.append("  No ATOM/HETATM nucleotide residues were detected.")
        return "\n".join(lines)

    for chain_id, residues in chains.items():
        if not residues:
            continue
        seq = residue_sequence(residues)
        first = residues[0].res_seq
        last = residues[-1].res_seq
        purines = sum(1 for res in residues if res.base in PURINES)
        pyrimidines = sum(1 for res in residues if res.base in PYRIMIDINES)
        lines.append(
            f"  Chain {display_chain_id(chain_id)}: {len(residues)} residue(s), "
            f"resSeq {first}..{last}, sequence {seq}"
        )
        lines.append(
            f"    residue numbers: {format_residue_number_list(residues)}"
        )
        lines.append(
            f"    base counts: purines={purines}, pyrimidines={pyrimidines}"
        )

    lines.extend([
        "",
        "Triplex convention used by this script:",
        "  Triples are written as Z·X-Y, where X-Y is the Watson-Crick duplex base pair and Z is the third-strand base.",
        "  Strand I   = purine strand in the duplex, i.e. X in Z·X-Y.",
        "  Strand II  = pyrimidine Watson-Crick partner in the duplex, i.e. Y in Z·X-Y.",
        "  Strand III = Hoogsteen or reverse-Hoogsteen strand added by this script, i.e. Z in Z·X-Y.",
        "",
        "Mode-specific sequence requirement:",
        "  antiparallel G·G-C: only the selected range on strand I must be G; strand II must provide the paired C segment.",
        "  parallel T·A-T: only the selected range on strand I must be A; strand II must provide the paired T segment.",
    ])
    return "\n".join(lines)


def format_selection_preview(
    pdb_path: str,
    strand_i_chain: str,
    residue_range: Optional[Tuple[int, int]],
    mode: str,
) -> str:
    """Return a short preview of the selected strand-I range for GUI assistance."""
    if not strand_i_chain or residue_range is None:
        return ""

    records = read_pdb_file(pdb_path)
    chains = build_chain_residues(records)
    if strand_i_chain not in chains:
        available = ", ".join(display_chain_id(cid) for cid in chains) or "<none>"
        return (
            "\nSelected range preview:\n"
            f"  Chain {display_chain_id(strand_i_chain)} was not found. Available chains: {available}"
        )

    start_res, end_res = residue_range
    low, high = sorted((start_res, end_res))
    selected = [res for res in chains[strand_i_chain] if low <= res.res_seq <= high]
    if not selected:
        return (
            "\nSelected range preview:\n"
            f"  No residues found on chain {display_chain_id(strand_i_chain)} in range {start_res}:{end_res}."
        )

    lines = [
        "",
        "Selected range preview:",
        f"  Strand I chain: {display_chain_id(strand_i_chain)}",
        f"  Residues: {selected[0].res_seq}..{selected[-1].res_seq}",
        f"  Sequence in selected range: {residue_sequence(selected)}",
    ]

    try:
        normalized_mode = normalize_mode(mode)
        expected = MOTIFS[normalized_mode].x_base
        bad = [res.label for res in selected if res.base != expected]
        if bad:
            lines.append(
                f"  WARNING: {normalized_mode} mode requires only the selected range on strand I to be all {expected}; "
                f"nonmatching residue(s): {', '.join(bad)}"
            )
        else:
            lines.append(
                f"  OK for {normalized_mode} mode: all selected strand-I residues are {expected}."
            )
    except Exception as exc:
        lines.append(f"  Mode check not available: {exc}")

    return "\n".join(lines)



def build_chain_residues(records: Sequence[object]) -> Dict[str, List[Residue]]:
    chains: Dict[str, List[Residue]] = {}
    residue_lookup: Dict[Tuple[str, int], Residue] = {}
    for rec in records:
        if rec.recordName not in ("ATOM", "HETATM"):
            continue
        key = (rec.chainID, rec.resSeq)
        if key not in residue_lookup:
            chain_list = chains.setdefault(rec.chainID, [])
            residue = Residue(
                chain_id=rec.chainID,
                res_seq=rec.resSeq,
                res_name=rec.resName.strip(),
                index_in_chain=len(chain_list),
            )
            chain_list.append(residue)
            residue_lookup[key] = residue
        residue_lookup[key].records.append(rec)
    return chains


def residue_coords(residue: Residue, atom_names: Sequence[str]) -> np.ndarray:
    amap = residue.atom_map()
    missing = [name for name in atom_names if name not in amap]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Residue {residue.label} is missing required atoms: {missing_str}")
    coords = [[amap[name].x, amap[name].y, amap[name].z] for name in atom_names]
    return np.asarray(coords, dtype=float)


def concat_base_coords(
    residue_x: Residue,
    base_x: str,
    residue_y: Residue,
    base_y: str,
) -> np.ndarray:
    atoms_x = ALIGN_ATOMS[base_x]
    atoms_y = ALIGN_ATOMS[base_y]
    coords_x = residue_coords(residue_x, atoms_x)
    coords_y = residue_coords(residue_y, atoms_y)
    return np.vstack([coords_x, coords_y])


def kabsch_transform(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return rotation R, translation t, and RMSD for mapping mobile -> target."""
    if mobile.shape != target.shape:
        raise ValueError(
            f"Point clouds must have the same shape, got {mobile.shape} and {target.shape}."
        )
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Point clouds must be shaped as (N, 3).")
    if mobile.shape[0] < 3:
        raise ValueError("At least three points are required for superposition.")

    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    mobile0 = mobile - mobile_center
    target0 = target - target_center

    covariance = mobile0.T @ target0
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T

    translation = target_center - rotation @ mobile_center
    moved = (rotation @ mobile.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((moved - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def transform_record(record: object, rotation: np.ndarray, translation: np.ndarray) -> object:
    new_record = copy.deepcopy(record)
    xyz = np.asarray([record.x, record.y, record.z], dtype=float)
    new_xyz = rotation @ xyz + translation
    new_record.update_xyz(*new_xyz.tolist())
    return new_record


def choose_unused_chain_id(records: Sequence[object], preferred: Optional[str] = None) -> str:
    used = {rec.chainID for rec in records if hasattr(rec, "chainID") and rec.chainID.strip()}
    if preferred and preferred not in used:
        return preferred
    pool = string.ascii_uppercase + string.ascii_lowercase + string.digits
    for candidate in pool:
        if candidate not in used:
            return candidate
    raise ValueError("Could not find an unused chain ID for strand III.")


def make_ter_for_last_residue(last_residue_records: Sequence[object]) -> object:
    if not last_residue_records:
        raise ValueError("Cannot create TER for an empty residue.")
    last = last_residue_records[-1]
    ter = pdb_ter_record("TER\n")
    ter.update_resName(last.resName)
    ter.update_chainID(last.chainID)
    ter.update_resSeq(last.resSeq)
    return ter


def get_link_distance(prev_residue_records: Sequence[object], next_residue_records: Sequence[object]) -> Optional[float]:
    prev = {rec.name.strip(): rec for rec in prev_residue_records}
    nxt = {rec.name.strip(): rec for rec in next_residue_records}
    if "O3'" not in prev or "P" not in nxt:
        return None
    o3 = np.asarray([prev["O3'"].x, prev["O3'"].y, prev["O3'"].z], dtype=float)
    p = np.asarray([nxt["P"].x, nxt["P"].y, nxt["P"].z], dtype=float)
    return float(np.linalg.norm(o3 - p))


MOTIFS: Dict[str, Motif] = {
    "antiparallel": Motif(
        mode="antiparallel",
        label="G·G-C",
        x_base="G",
        y_base="C",
        z_base="G",
        template_text=TEMPLATE_G_G_C,
    ),
    "parallel": Motif(
        mode="parallel",
        label="T·A-T",
        x_base="A",
        y_base="T",
        z_base="T",
        template_text=TEMPLATE_T_A_T,
    ),
}

# ---------------------------------------------------------------------------
# Pairing / chain detection
# ---------------------------------------------------------------------------


@dataclass
class PairingCandidate:
    chain_ii: str
    window_start_index: int
    order: str
    score: float
    strand_ii_residues: List[Residue]


@dataclass
class ConversionResult:
    output_path: str
    mode: str
    motif_label: str
    strand_i: str
    strand_ii: str
    strand_iii: str
    strand_i_residues: List[Residue]
    strand_ii_residues: List[Residue]
    strand_iii_resseqs: List[int]
    per_step_rmsd: List[float]
    link_distances: List[Optional[float]]


def validate_strand_i_residues(residues: Sequence[Residue], expected_base: str) -> None:
    bad = [res.label for res in residues if res.base != expected_base]
    if bad:
        joined = ", ".join(bad)
        raise ValueError(
            f"The selected strand I range must be all {expected_base} for this mode, but found: {joined}. "
            "The full strand-I sequence may be mixed; only the requested conversion range is checked."
        )


def score_pairing_window(
    strand_i_residues: Sequence[Residue],
    strand_ii_window: Sequence[Residue],
    motif: Motif,
    order: str,
) -> Optional[PairingCandidate]:
    if order == "same":
        paired_ii = list(strand_ii_window)
    elif order == "reverse":
        paired_ii = list(reversed(strand_ii_window))
    else:
        raise ValueError(f"Unknown order '{order}'.")

    rmsds: List[float] = []
    try:
        for res_i, res_ii in zip(strand_i_residues, paired_ii):
            target = concat_base_coords(res_i, motif.x_base, res_ii, motif.y_base)
            _, _, rmsd = kabsch_transform(motif.ref_pair_coords, target)
            rmsds.append(rmsd)
    except ValueError:
        return None

    return PairingCandidate(
        chain_ii=paired_ii[0].chain_id,
        window_start_index=strand_ii_window[0].index_in_chain,
        order=order,
        score=float(np.mean(rmsds)),
        strand_ii_residues=paired_ii,
    )


def find_best_strand_ii_mapping(
    chains: Dict[str, List[Residue]],
    strand_i_chain: str,
    strand_i_residues: Sequence[Residue],
    motif: Motif,
    requested_strand_ii: Optional[str] = None,
) -> PairingCandidate:
    if requested_strand_ii:
        candidate_chain_ids = [requested_strand_ii]
    else:
        candidate_chain_ids = [cid for cid in chains if cid != strand_i_chain]

    if not candidate_chain_ids:
        raise ValueError("Could not find any candidate chain for strand II.")

    window_size = len(strand_i_residues)
    candidates: List[PairingCandidate] = []

    for chain_id in candidate_chain_ids:
        residues = chains.get(chain_id, [])
        if len(residues) < window_size:
            continue
        for start in range(0, len(residues) - window_size + 1):
            window = residues[start : start + window_size]
            if any(res.base != motif.y_base for res in window):
                continue
            for order in ("reverse", "same"):
                candidate = score_pairing_window(strand_i_residues, window, motif, order)
                if candidate is not None:
                    candidates.append(candidate)

    if not candidates:
        base = motif.y_base
        if requested_strand_ii:
            raise ValueError(
                f"Could not map strand II on chain {requested_strand_ii}. "
                f"Expected a contiguous {base}-only segment of length {window_size}."
            )
        raise ValueError(
            f"Could not auto-detect strand II. Expected a contiguous {base}-only segment of length {window_size}."
        )

    best = min(candidates, key=lambda item: (item.score, item.chain_ii, item.window_start_index, item.order))
    if best.score > 3.5:
        raise ValueError(
            f"The best strand II match ({best.chain_ii}, {best.order} order) fits poorly "
            f"to the reference motif (mean RMSD {best.score:.2f} Å)."
        )
    return best


# ---------------------------------------------------------------------------
# Conversion core
# ---------------------------------------------------------------------------


def convert_duplex_to_triplex(
    duplex_path: str,
    output_path: str,
    strand_i_chain: str,
    residue_range: Tuple[int, int],
    mode: str,
    strand_ii_chain: Optional[str] = None,
    strand_iii_chain: Optional[str] = None,
    strand_iii_start_resseq: int = 1,
) -> ConversionResult:
    mode = normalize_mode(mode)
    motif = MOTIFS[mode]

    duplex_records = read_pdb_file(duplex_path)
    chains = build_chain_residues(duplex_records)

    if strand_i_chain not in chains:
        available = ", ".join(sorted(chains)) or "<none>"
        raise ValueError(f"Chain {strand_i_chain} was not found. Available chains: {available}")

    start_res, end_res = residue_range
    low, high = sorted((start_res, end_res))
    strand_i_residues = [res for res in chains[strand_i_chain] if low <= res.res_seq <= high]
    if not strand_i_residues:
        raise ValueError(
            f"No residues from chain {strand_i_chain} were found in the requested range {start_res}:{end_res}."
        )

    validate_strand_i_residues(strand_i_residues, motif.x_base)

    pairing = find_best_strand_ii_mapping(
        chains=chains,
        strand_i_chain=strand_i_chain,
        strand_i_residues=strand_i_residues,
        motif=motif,
        requested_strand_ii=strand_ii_chain,
    )
    strand_ii_residues = pairing.strand_ii_residues

    if strand_iii_chain:
        used = {rec.chainID for rec in duplex_records if hasattr(rec, "chainID")}
        if strand_iii_chain in used:
            raise ValueError(
                f"Requested strand III chain ID '{strand_iii_chain}' already exists in the input PDB. "
                f"Choose a different chain ID or leave it blank to auto-select one."
            )
    else:
        strand_iii_chain = choose_unused_chain_id(duplex_records, preferred="C")

    build_pairs: List[Tuple[Residue, Residue]] = list(zip(strand_i_residues, strand_ii_residues))
    if mode == "antiparallel":
        build_pairs = list(reversed(build_pairs))

    third_strand_records: List[object] = []
    third_strand_resseqs: List[int] = []
    per_step_rmsd: List[float] = []
    residue_blocks: List[List[object]] = []

    for offset, (res_i, res_ii) in enumerate(build_pairs):
        target_coords = concat_base_coords(res_i, motif.x_base, res_ii, motif.y_base)
        rotation, translation, rmsd = kabsch_transform(motif.ref_pair_coords, target_coords)
        new_resseq = strand_iii_start_resseq + offset

        transformed_block: List[object] = []
        for template_record in motif.z_residue.records:
            new_record = transform_record(template_record, rotation, translation)
            new_record.update_chainID(strand_iii_chain)
            new_record.update_resSeq(new_resseq)
            transformed_block.append(new_record)

        third_strand_records.extend(transformed_block)
        residue_blocks.append(transformed_block)
        third_strand_resseqs.append(new_resseq)
        per_step_rmsd.append(rmsd)

    if third_strand_records:
        third_strand_records.append(make_ter_for_last_residue(residue_blocks[-1]))

    all_records = list(duplex_records) + third_strand_records
    with open(output_path, "w", encoding="utf-8") as handle:
        rec2file(all_records, handle, reorder_serial=True)
        handle.write("END\n")

    link_distances: List[Optional[float]] = []
    for prev_block, next_block in zip(residue_blocks[:-1], residue_blocks[1:]):
        link_distances.append(get_link_distance(prev_block, next_block))

    return ConversionResult(
        output_path=output_path,
        mode=mode,
        motif_label=motif.label,
        strand_i=strand_i_chain,
        strand_ii=pairing.chain_ii,
        strand_iii=strand_iii_chain,
        strand_i_residues=list(strand_i_residues),
        strand_ii_residues=list(strand_ii_residues),
        strand_iii_resseqs=third_strand_resseqs,
        per_step_rmsd=per_step_rmsd,
        link_distances=link_distances,
    )


def format_result_summary(result: ConversionResult) -> str:
    lines = [
        f"Triplex written to: {result.output_path}",
        f"Mode: {result.mode} ({result.motif_label})",
        f"Strand I  (purine; X in Z·X-Y): {result.strand_i}",
        f"Strand II (WC partner; Y):      {result.strand_ii}",
        f"Strand III (Hoogsteen; Z):      {result.strand_iii}",
        f"Converted residues: {result.strand_i_residues[0].res_seq} to {result.strand_i_residues[-1].res_seq} on strand I",
        f"Strand I selected sequence: {residue_sequence(result.strand_i_residues)}",
        f"Strand II paired sequence:   {residue_sequence(result.strand_ii_residues)}",
        f"Added {len(result.strand_iii_resseqs)} residue(s) to strand III with resSeq {result.strand_iii_resseqs[0]}..{result.strand_iii_resseqs[-1]}",
        f"Mean template-fit RMSD: {np.mean(result.per_step_rmsd):.3f} Å",
    ]
    if result.link_distances:
        numeric = [d for d in result.link_distances if d is not None]
        if numeric:
            lines.append(
                f"Mean O3'-P linkage distance in strand III: {np.mean(numeric):.3f} Å"
            )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI and GUI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a DNA duplex PDB into a triplex over a residue range on strand I. "
            "Supported motifs: antiparallel G·G-C and parallel T·A-T."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python convert_to_triplex_pdbV2_1.py duplex.pdb --strand-I A --range 10:18 --mode antiparallel\n"
            "  python convert_to_triplex_pdbV2_1.py duplex.pdb --strand-I A --range 7-12 --mode parallel --strand-III D\n"
            "  python convert_to_triplex_pdbV2_1.py --gui\n"
        ),
    )
    parser.add_argument("duplex", help="input duplex PDB file")
    parser.add_argument(
        "-I",
        "--strand-I",
        dest="strand_i",
        required=True,
        help="chain ID of strand I (the purine strand)",
    )
    parser.add_argument(
        "-r",
        "--range",
        required=True,
        help="residue range on strand I, e.g. 10:18 or 10-18",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        help="triplex mode: antiparallel or parallel",
    )
    parser.add_argument(
        "-II",
        "--strand-II",
        dest="strand_ii",
        default=None,
        help="optional chain ID of strand II (auto-detected if omitted)",
    )
    parser.add_argument(
        "-III",
        "--strand-III",
        dest="strand_iii",
        default=None,
        help="optional chain ID for strand III (auto-selected if omitted)",
    )
    parser.add_argument(
        "--strand-III-start",
        dest="strand_iii_start",
        type=int,
        default=1,
        help="starting residue number for strand III (default: 1)",
    )
    parser.add_argument(
        "-o",
        "--out",
        default=None,
        help="output PDB file (default: input filename with _2TH before the extension)",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="launch the GUI",
    )
    return parser


def run_gui() -> None:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except Exception as exc:  # pragma: no cover - GUI availability depends on environment
        raise RuntimeError(
            "GUI mode requested, but tkinter is not available in this Python environment."
        ) from exc

    root = tk.Tk()
    icon_path = Path(__file__).resolve().parents[1] / "assets" / "bnp_na_icon.png"
    if icon_path.exists():
        try:
            root._bnp_na_window_icon = tk.PhotoImage(file=str(icon_path))
            root.iconphoto(True, root._bnp_na_window_icon)
        except Exception:
            pass
    root.title("convert_to_triplex_pdbV2_1")
    root.geometry("1000x760")
    root.minsize(940, 680)

    duplex_var = tk.StringVar()
    out_var = tk.StringVar(value="")
    strand_i_var = tk.StringVar(value="A")
    strand_ii_var = tk.StringVar(value="")
    strand_iii_var = tk.StringVar(value="")
    range_start_var = tk.StringVar(value="1")
    range_end_var = tk.StringVar(value="1")
    mode_var = tk.StringVar(value="antiparallel")
    strand_iii_start_var = tk.StringVar(value="1")
    last_auto_output = {"path": ""}

    main = ttk.Frame(root, padding=14)
    main.pack(fill="both", expand=True)
    main.columnconfigure(0, weight=0)
    main.columnconfigure(1, weight=1)
    main.columnconfigure(2, weight=0)
    main.rowconfigure(11, weight=1)

    def set_log(text: str) -> None:
        output_box.configure(state="normal")
        output_box.delete("1.0", "end")
        output_box.insert("1.0", text)
        output_box.configure(state="disabled")

    def append_log(text: str) -> None:
        output_box.configure(state="normal")
        output_box.insert("end", text)
        output_box.see("end")
        output_box.configure(state="disabled")

    def update_default_output(force: bool = False) -> None:
        duplex = duplex_var.get().strip()
        if not duplex:
            return
        proposed = default_output_path(duplex)
        current = out_var.get().strip()
        if force or not current or current == last_auto_output["path"] or current == "triplex_V2.pdb":
            out_var.set(proposed)
            last_auto_output["path"] = proposed

    def current_range_or_none() -> Optional[Tuple[int, int]]:
        try:
            return (int(range_start_var.get().strip()), int(range_end_var.get().strip()))
        except Exception:
            return None

    def load_sequence_info(event: Optional[object] = None) -> None:
        duplex = duplex_var.get().strip()
        if not duplex:
            set_log(
                "Choose an input duplex PDB file to list the detected strands/chains and sequences.\n"
                "The output filename will default to the input filename with _2TH inserted before the extension."
            )
            return

        update_default_output(force=False)
        if not os.path.isfile(duplex):
            set_log(f"Input file not found: {duplex}\n")
            return

        try:
            text = format_chain_sequence_info(duplex)
            preview = format_selection_preview(
                pdb_path=duplex,
                strand_i_chain=strand_i_var.get().strip(),
                residue_range=current_range_or_none(),
                mode=mode_var.get().strip(),
            )
            set_log(text + preview)
        except Exception as exc:
            set_log(f"Could not read strand/sequence information from {duplex}:\n{exc}\n")

    def choose_input() -> None:
        selected = filedialog.askopenfilename(
            title="Select duplex PDB",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if selected:
            duplex_var.set(selected)
            update_default_output(force=True)
            load_sequence_info()

    def choose_output() -> None:
        initial = out_var.get().strip() or default_output_path(duplex_var.get().strip() or "triplex.pdb")
        selected = filedialog.asksaveasfilename(
            title="Save output PDB",
            initialfile=os.path.basename(initial),
            initialdir=os.path.dirname(initial) or None,
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if selected:
            out_var.set(selected)

    ttk.Label(main, text="Input duplex PDB").grid(row=0, column=0, sticky="w", pady=4)
    input_entry = ttk.Entry(main, textvariable=duplex_var, width=82)
    input_entry.grid(row=0, column=1, sticky="ew", pady=4)
    input_entry.bind("<Return>", load_sequence_info)
    input_entry.bind("<FocusOut>", load_sequence_info)
    ttk.Button(main, text="Browse/load", command=choose_input).grid(
        row=0, column=2, sticky="ew", padx=(8, 0), pady=4
    )

    ttk.Label(main, text="Output PDB").grid(row=1, column=0, sticky="w", pady=4)
    ttk.Entry(main, textvariable=out_var, width=82).grid(row=1, column=1, sticky="ew", pady=4)
    ttk.Button(main, text="Browse", command=choose_output).grid(
        row=1, column=2, sticky="ew", padx=(8, 0), pady=4
    )

    ttk.Label(main, text="Strand I (purine chain)").grid(row=2, column=0, sticky="w", pady=4)
    strand_i_entry = ttk.Entry(main, textvariable=strand_i_var, width=12)
    strand_i_entry.grid(row=2, column=1, sticky="w", pady=4)
    strand_i_entry.bind("<Return>", load_sequence_info)
    strand_i_entry.bind("<FocusOut>", load_sequence_info)

    ttk.Label(main, text="Residue range on strand I").grid(row=3, column=0, sticky="w", pady=4)
    range_box = ttk.Frame(main)
    range_box.grid(row=3, column=1, sticky="w")
    start_entry = ttk.Entry(range_box, textvariable=range_start_var, width=10)
    start_entry.pack(side="left")
    ttk.Label(range_box, text=" to ").pack(side="left")
    end_entry = ttk.Entry(range_box, textvariable=range_end_var, width=10)
    end_entry.pack(side="left")
    for entry in (start_entry, end_entry):
        entry.bind("<Return>", load_sequence_info)
        entry.bind("<FocusOut>", load_sequence_info)

    ttk.Label(main, text="Triplex mode").grid(row=4, column=0, sticky="w", pady=4)
    mode_box = ttk.Combobox(
        main,
        textvariable=mode_var,
        values=["antiparallel", "parallel"],
        state="readonly",
        width=16,
    )
    mode_box.grid(row=4, column=1, sticky="w", pady=4)
    mode_box.bind("<<ComboboxSelected>>", load_sequence_info)

    ttk.Label(main, text="Strand II chain (optional)").grid(row=5, column=0, sticky="w", pady=4)
    ttk.Entry(main, textvariable=strand_ii_var, width=12).grid(row=5, column=1, sticky="w", pady=4)

    ttk.Label(main, text="Strand III chain (optional)").grid(row=6, column=0, sticky="w", pady=4)
    ttk.Entry(main, textvariable=strand_iii_var, width=12).grid(row=6, column=1, sticky="w", pady=4)

    ttk.Label(main, text="Strand III first resSeq").grid(row=7, column=0, sticky="w", pady=4)
    ttk.Entry(main, textvariable=strand_iii_start_var, width=12).grid(row=7, column=1, sticky="w", pady=4)

    explanation = (
        "Convention: triples are written as Z·X-Y. X-Y is the Watson-Crick base pair already in the duplex; "
        "Z is the third-strand Hoogsteen/reverse-Hoogsteen base. Strand I is the purine duplex strand (X), "
        "strand II is the pyrimidine Watson-Crick partner (Y), and strand III is the added Hoogsteen strand (Z).\n"
        "Supported templates: antiparallel G·G-C and parallel T·A-T. Only the specified residue range on strand I "
        "must be G for antiparallel or A for parallel; the rest of strand I may contain other bases. Leave strand II "
        "blank to auto-detect the paired segment, and leave strand III blank to auto-assign an unused chain ID."
    )
    ttk.Label(main, text=explanation, justify="left", wraplength=930).grid(
        row=8, column=0, columnspan=3, sticky="ew", pady=(12, 8)
    )

    log_label = ttk.Label(main, text="GUI log / strand and sequence information")
    log_label.grid(row=9, column=0, columnspan=3, sticky="w", pady=(4, 2))

    log_frame = ttk.Frame(main)
    log_frame.grid(row=11, column=0, columnspan=3, sticky="nsew", pady=(0, 8))
    log_frame.columnconfigure(0, weight=1)
    log_frame.rowconfigure(0, weight=1)
    output_box = tk.Text(log_frame, height=18, width=112, wrap="word")
    output_box.grid(row=0, column=0, sticky="nsew")
    yscroll = ttk.Scrollbar(log_frame, orient="vertical", command=output_box.yview)
    yscroll.grid(row=0, column=1, sticky="ns")
    output_box.configure(yscrollcommand=yscroll.set)

    def run_from_gui() -> None:
        try:
            duplex = duplex_var.get().strip()
            update_default_output(force=False)
            out = out_var.get().strip()
            if not duplex:
                raise ValueError("Please choose an input duplex PDB file.")
            if not os.path.isfile(duplex):
                raise ValueError(f"Input file not found: {duplex}")
            if not out:
                out = default_output_path(duplex)
                out_var.set(out)

            residue_range = (int(range_start_var.get().strip()), int(range_end_var.get().strip()))
            set_log(format_chain_sequence_info(duplex) + format_selection_preview(
                pdb_path=duplex,
                strand_i_chain=strand_i_var.get().strip(),
                residue_range=residue_range,
                mode=mode_var.get().strip(),
            ) + "\n\nRunning conversion...\n")

            result = convert_duplex_to_triplex(
                duplex_path=duplex,
                output_path=out,
                strand_i_chain=strand_i_var.get().strip(),
                residue_range=residue_range,
                mode=mode_var.get().strip(),
                strand_ii_chain=(strand_ii_var.get().strip() or None),
                strand_iii_chain=(strand_iii_var.get().strip() or None),
                strand_iii_start_resseq=int(strand_iii_start_var.get().strip()),
            )
            summary = format_result_summary(result)
            append_log("\n" + summary + "\n")
            messagebox.showinfo("Success", summary)
        except Exception as exc:
            append_log(f"\nError: {exc}\n")
            messagebox.showerror("Error", str(exc))

    button_row = ttk.Frame(main)
    button_row.grid(row=12, column=0, columnspan=3, sticky="ew")
    ttk.Button(button_row, text="Refresh strand/sequence info", command=load_sequence_info).pack(side="left")
    ttk.Button(button_row, text="Run conversion", command=run_from_gui).pack(side="left", padx=(8, 0))
    ttk.Button(button_row, text="Quit", command=root.destroy).pack(side="left", padx=(8, 0))

    set_log(
        "Choose an input duplex PDB file to list chains and sequences here.\n"
        "Default output naming: input_filename_2TH.pdb."
    )
    root.mainloop()

def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if not argv or "--gui" in argv:
        try:
            run_gui()
            return 0
        except Exception as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

    parser = build_parser()
    args = parser.parse_args(argv)
    residue_range = parse_residue_range(args.range)

    try:
        result = convert_duplex_to_triplex(
            duplex_path=args.duplex,
            output_path=args.out or default_output_path(args.duplex),
            strand_i_chain=args.strand_i,
            residue_range=residue_range,
            mode=args.mode,
            strand_ii_chain=args.strand_ii,
            strand_iii_chain=args.strand_iii,
            strand_iii_start_resseq=args.strand_iii_start,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(format_result_summary(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
