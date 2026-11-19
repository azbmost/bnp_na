"""
core_PDNA.py

P-DNA template level for the bnp_na P-DNA builder (bnp_na V13.22).

P-DNA here is the Pauling-like DNA of the Pf1 filamentous phage (PDB 1PFI):
both strands wind around the helix axis with their phosphates at the center
(P at r = 2.50 A) and every base pointing outward (C1' at r = 4.28 A), so no
base pairs are formed.

Template analysis of 1PFI-DNA_only_ligated.pdb (chains C and D, each poly(dC),
35 nt, written 5' -> 3'):

  * The model is exactly helical. One screw step maps nucleotide i onto i+1
    of chain C (and i+1 onto i of chain D): twist 131.840 deg and rise
    6.100 A per nucleotide about the z axis through the origin (fit RMSD
    0.0006 A over all 1292 mapped atoms, i.e. coordinate rounding), right
    handed, 2.731 nt per turn, pitch 16.66 A.
  * A dyad on the x axis (180 deg about x) maps chain C onto chain D exactly.
  * Registering the strands by nucleotide level (C1' height), C_i sits at the
    level of D_(34-i), 0.31 A apart in z. Chain C therefore pairs by level with
    D33..D1 over C1..C33, and each strand carries a 2-nt 3' overhang: C34-C35
    above D1 at the top, and D34-D35 below C1 at the bottom. In the original
    deposition order (C1 with D35, C2 with D34, ...) the strands are offset by
    exactly these two levels.
  * The two nucleotides of one level sit at C1' azimuths -114.69 deg (strand
    C, running 5'->3' along +z) and +114.69 deg (strand D), the same 130.6 deg
    C1'-axis-C1' angle as a Watson-Crick pair in bnp_na's B-DNA.

The level C17/D17 is kept here. Its dyad is the x axis and its two C1' atoms
are at z = -0.157 and +0.157, so the level center is the origin. The builder
places level k at Screw^(k-1)(C17/D17): strand I (chain A) takes chain C's
geometry and runs 5'->3' along +z, strand II (chain B) takes chain D's
geometry and runs 5'->3' along -z, and the bases are replaced by the
requested ones in the standard base reference frame (standard_bases.py).
"""

from __future__ import annotations

import math

TEMPLATE_NAME = "1PFI-DNA_only_ligated.pdb (Pf1 phage DNA, PDB 1PFI)"
TEMPLATE_BASE = "C"

# Exact helical symmetry of the template, per nucleotide along one strand.
TWIST_DEG = 131.84
RISE = 6.10
NT_PER_TURN = 360.0 / TWIST_DEG
PITCH = NT_PER_TURN * RISE
P_RADIUS = 2.50
C1_RADIUS = 4.28
# Angle about the axis between the two C1' atoms of one level, through the
# narrow side (as for the minor-groove side of a Watson-Crick pair).
LEVEL_C1_ANGLE_DEG = 130.62
# Height difference between the two C1' atoms of one level.
LEVEL_C1_DZ = 0.314

STRAND_I_CHAIN = "C"
STRAND_II_CHAIN = "D"

PDB_LINES = [
    "ATOM    305  P   DC  C  17      -2.470   0.386  -4.575  1.00  0.00           P\n",
    "ATOM    306  OP1 DC  C  17      -3.206  -0.707  -5.366  1.00  0.00           O\n",
    "ATOM    307  OP2 DC  C  17      -3.556   1.399  -4.178  1.00  0.00           O\n",
    "ATOM    308  O5' DC  C  17      -1.784  -0.243  -3.217  1.00  0.00           O\n",
    "ATOM    309  C5' DC  C  17      -2.026  -1.578  -2.696  1.00  0.00           C\n",
    "ATOM    310  C4' DC  C  17      -1.673  -1.742  -1.181  1.00  0.00           C\n",
    "ATOM    311  O4' DC  C  17      -2.557  -2.763  -0.634  1.00  0.00           O\n",
    "ATOM    312  C3' DC  C  17      -0.237  -2.188  -0.771  1.00  0.00           C\n",
    "ATOM    313  O3' DC  C  17       0.061  -1.737   0.589  1.00  0.00           O\n",
    "ATOM    314  C2' DC  C  17      -0.403  -3.711  -0.790  1.00  0.00           C\n",
    "ATOM    315  C1' DC  C  17      -1.789  -3.892  -0.157  1.00  0.00           C\n",
    "ATOM    316  N1  DC  C  17      -2.448  -5.182  -0.521  1.00  0.00           N\n",
    "ATOM    317  C2  DC  C  17      -2.461  -6.290   0.302  1.00  0.00           C\n",
    "ATOM    318  O2  DC  C  17      -1.856  -6.365   1.368  1.00  0.00           O\n",
    "ATOM    319  N3  DC  C  17      -3.262  -7.495  -0.131  1.00  0.00           N\n",
    "ATOM    320  C4  DC  C  17      -3.869  -7.505  -1.271  1.00  0.00           C\n",
    "ATOM    321  N4  DC  C  17      -4.611  -8.709  -1.604  1.00  0.00           N\n",
    "ATOM    322  C5  DC  C  17      -3.845  -6.456  -2.147  1.00  0.00           C\n",
    "ATOM    323  C6  DC  C  17      -3.168  -5.311  -1.855  1.00  0.00           C\n",
    "ATOM    971  P   DC  D  17      -2.470  -0.386   4.575  1.00  0.00           P\n",
    "ATOM    972  OP1 DC  D  17      -3.206   0.707   5.366  1.00  0.00           O\n",
    "ATOM    973  OP2 DC  D  17      -3.556  -1.399   4.178  1.00  0.00           O\n",
    "ATOM    974  O5' DC  D  17      -1.784   0.243   3.217  1.00  0.00           O\n",
    "ATOM    975  C5' DC  D  17      -2.026   1.578   2.696  1.00  0.00           C\n",
    "ATOM    976  C4' DC  D  17      -1.673   1.742   1.181  1.00  0.00           C\n",
    "ATOM    977  O4' DC  D  17      -2.557   2.763   0.634  1.00  0.00           O\n",
    "ATOM    978  C3' DC  D  17      -0.237   2.188   0.771  1.00  0.00           C\n",
    "ATOM    979  O3' DC  D  17       0.061   1.737  -0.589  1.00  0.00           O\n",
    "ATOM    980  C2' DC  D  17      -0.403   3.711   0.790  1.00  0.00           C\n",
    "ATOM    981  C1' DC  D  17      -1.789   3.892   0.157  1.00  0.00           C\n",
    "ATOM    982  N1  DC  D  17      -2.448   5.182   0.521  1.00  0.00           N\n",
    "ATOM    983  C2  DC  D  17      -2.461   6.290  -0.302  1.00  0.00           C\n",
    "ATOM    984  O2  DC  D  17      -1.856   6.365  -1.368  1.00  0.00           O\n",
    "ATOM    985  N3  DC  D  17      -3.262   7.495   0.131  1.00  0.00           N\n",
    "ATOM    986  C4  DC  D  17      -3.869   7.505   1.271  1.00  0.00           C\n",
    "ATOM    987  N4  DC  D  17      -4.611   8.709   1.604  1.00  0.00           N\n",
    "ATOM    988  C5  DC  D  17      -3.845   6.456   2.147  1.00  0.00           C\n",
    "ATOM    989  C6  DC  D  17      -3.168   5.311   1.855  1.00  0.00           C\n",
]


def level_atoms():
    """Return {'I': [(name, (x, y, z)), ...], 'II': [...]} for the template level."""
    out = {"I": [], "II": []}
    for line in PDB_LINES:
        chain = line[21]
        key = "I" if chain == STRAND_I_CHAIN else "II"
        out[key].append((line[12:16].strip(), (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return out


def helical_parameters():
    """The fixed P-DNA helical parameters shown in the GUI (not adjustable)."""
    return {
        "h-Rise (Å/nt)": RISE,
        "h-Twist (°/nt)": TWIST_DEG,
        "nt per turn": NT_PER_TURN,
        "Pitch (Å)": PITCH,
        "P radius (Å)": P_RADIUS,
        "C1' radius (Å)": C1_RADIUS,
    }


def screw_matrix(steps: float):
    """Rotation (3x3 nested lists) and z shift for ``steps`` screw steps."""
    angle = math.radians(TWIST_DEG * steps)
    c, s = math.cos(angle), math.sin(angle)
    return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], RISE * steps
