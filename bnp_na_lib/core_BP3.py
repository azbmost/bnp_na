"""
core_BP3.py

B-P junction core with a three-nucleotide transition, for make_BP.py
--transition 3 (bnp_na V13.23): two B-DNA base pairs, three transition
levels and two P-DNA levels.

There is no experimental structure of a B-DNA/P-DNA junction, so this core is a
model. It was built once, offline, and is bundled here so that a junction build
is only rigid fitting, trimming, base grafting and ligation, which is immediate.

Why a second core: in the two-nucleotide core (core_BP.py) strand I swings back
about the helix axis as it enters the P-DNA (one backbone step of -134 deg), so
a junction adds only about 25 deg of twist. With two nucleotides this cannot be
avoided: the two strands' windings through a junction are tied together by the
B-DNA and P-DNA geometries (they differ by about 200 deg), strand II's step out
of the P-DNA is fixed near 120 deg by the P nucleotide it hangs from, and a
step leaving the B-DNA radius can turn at most about 45 deg. A one-nucleotide
transition cannot even be ligated on both strands. Three nucleotides per strand
give enough freedom for every backbone step to turn forward.

  * B part: the top two base pairs of a bnp_na default B-DNA (DSSR rebuild at
    the default parameters, phenix.geometry_minimization with min_P_C5.params,
    phosphate regularization), aligned to +Z.
  * P part: bnp_na P-DNA (1PFI template, core_PDNA.py) on the same axis. Its
    first P-form level (the P interior level, E6/F2) is centered 14.0 A above
    the top B base pair and turned 0 deg about the axis relative to the
    bnp_na P-DNA frame.
  * Transition: on each strand the three nucleotides between them, which take
    the place of the P-DNA's first three levels, were built by torsion-space
    loop closure. Each is the ideal DSSR C2'-endo sugar-base unit from B-DNA
    (both helices are C2'-endo: B-DNA P = 154 deg, the 1PFI P-DNA P = 176 deg)
    joined through phosphate linkages with Phenix's ideal bond lengths and
    angles. Backbone torsions and glycosidic rotations were optimized so each
    strand closes onto the next helix, every phosphate-to-phosphate step turns
    forward about the axis between the B-DNA twist (34 deg) and the P-DNA
    twist (132 deg), and every glycosidic bond stays anti (chi within 50 deg
    of the unit's -108 deg). Contacts were scored with the union of A, C, G
    and T at every base near the junction, including their ring, amino and
    methyl hydrogens, and with riding sugar hydrogens, so the core stays
    clash-free whatever bases are grafted onto it.
  * Two junctions back to back: a P-DNA may keep no P-form level between two
    junctions, in which case the upper junction is the dyad image of the
    lower one (the P-DNA dyad between its third and fourth levels) and the
    two transitions bond to each other directly, E5 O3' to the other core's
    F3 P. The transition was therefore refined against its own dyad image as
    well: the image of the transition (both strands, union bases with
    hydrogens) and of the B pairs was added to the scoring context, with the
    cross bond treated as a covalent neighbour, and the refinement was
    repeated three times, rebuilding the image from each result. The ideal
    P-DNA neighbour stayed in the context throughout, so the same core serves
    P-DNAs that keep one or more P-form levels.
  * Result: phosphate steps, in the +Z direction, strand I (E) 35, 47, 35 and
    50 deg from the last B pair to the first P-form level; strand II (F) 33,
    98, 105 and 132 deg. Strand I's phosphates move inward 9.5 -> 9.3 -> 6.7
    -> 6.6 -> 2.5 A, strand II's 9.5 -> 6.1 -> 3.8 -> 2.6 -> 2.5 A. Puckers P
    = 152-155 deg; chi -60 to -158 deg.
  * Finish: phenix.geometry_minimization moved only the transition residues,
    with every atom restrained to its start (sigma 0.08 A; largest shift
    0.12 A), to bring the closing O3'-P-O5' angles in line.
  * Checks on the transition residues, against Phenix restraints: bonds RMSD
    0.005 A (largest 0.9 sigma), angles RMSD 0.8 deg (largest 2.5 sigma),
    chirality as built; the eight O3'-P linkages into, within and out of the
    transition are 1.600-1.612 A. MolProbity clashscore is 0 for 44 base
    assignments (all A, all C, all G, all T and 40 random ones) of the
    transition, P and B residues around the junction, and for random-sequence
    B-P-B builds whose P-DNA keeps no, one or two P-form levels between the
    two junctions.

Residues (chain E runs 5' B-DNA -> 3' P-DNA, chain F runs 5' P-DNA -> 3' B-DNA):

  E1/F7  B outer pair           E2/F6  B interior pair
  E3/F5  transition level 1     E4/F4  transition level 2
  E5/F3  transition level 3 (P side)
  E6/F2  P interior level       E7/F1  P outer level

At a junction, make_BP.py fits E1/F7 + E2/F6 onto the B-DNA's two terminal
base pairs and E6/F2 + E7/F1 onto the P-DNA's fourth and fifth levels from
that end. The core's B pairs and transition levels then replace the B-DNA's
two terminal pairs and the P-DNA's three terminal levels and take over their
bases; the core's P levels are dropped, as the P-DNA keeps its own
nucleotides there. Where a 6-level P-DNA joins two B-DNAs, all six levels
become transition and the two cores bond to each other. The G/C bases stored
below are placeholders.

The offline generator is not part of the app; the scripts are kept with the
maintainer's working files.
"""

from __future__ import annotations

# Sugar atoms used for fitting (order matters), as in core_BZ.py.
SUGAR_ATOMS = ["C1'", "C2'", "C3'", "C4'", "O5'", "C5'"]

STRAND_E = "E"  # runs 5' B-DNA -> 3' P-DNA
STRAND_F = "F"  # runs 5' P-DNA -> 3' B-DNA

# Pairs are ((chain, resSeq), (chain, resSeq)) with the E residue first.
# B end, interior -> outer (as B_END_BASEPAIRS in core_BZ.py).
B_END_BASEPAIRS = [(("E", 2), ("F", 6)), (("E", 1), ("F", 7))]
# Transition levels, B side -> P side. They take the place of the P-DNA's
# terminal levels and take over their bases.
TRANSITION_LEVELS = [(("E", 3), ("F", 5)), (("E", 4), ("F", 4)), (("E", 5), ("F", 3))]
# P end, interior -> outer. Used only to position the P-DNA helix; the helix
# keeps its own nucleotides at these two levels.
P_END_LEVELS = [(("E", 6), ("F", 2)), (("E", 7), ("F", 1))]
# P-form levels that must separate two junctions on one P-DNA. None: the
# transition was refined against its own dyad image, which is what the other
# junction is when a P-DNA keeps no level between them, so two cores bond to
# each other directly (O3'-P 1.58 A) and random B-P-B builds with 0, 1 or 2
# P-form levels kept between the junctions have a MolProbity clashscore of 0.
MIN_KEPT_P_LEVELS_BETWEEN_JUNCTIONS = 0

PDB_LINES = [
    'ATOM      1  P    DG E   1       9.452   0.989  25.628  1.00  0.00           P\n',
    'ATOM      2  OP1  DG E   1      10.779   1.288  25.046  1.00  0.00           O\n',
    'ATOM      3  OP2  DG E   1       9.310  -0.096  26.625  1.00  0.00           O\n',
    "ATOM      4  O5'  DG E   1       8.860   2.306  26.255  1.00  0.00           O\n",
    "ATOM      5  C5'  DG E   1       8.186   3.153  25.397  1.00  0.00           C\n",
    "ATOM      6  C4'  DG E   1       6.951   3.739  26.038  1.00  0.00           C\n",
    "ATOM      7  O4'  DG E   1       5.779   2.880  25.927  1.00  0.00           O\n",
    "ATOM      8  C3'  DG E   1       7.019   4.079  27.527  1.00  0.00           C\n",
    "ATOM      9  O3'  DG E   1       6.557   5.338  27.823  1.00  0.00           O\n",
    "ATOM     10  C2'  DG E   1       6.255   2.948  28.218  1.00  0.00           C\n",
    "ATOM     11  C1'  DG E   1       5.134   2.841  27.191  1.00  0.00           C\n",
    'ATOM     12  N9   DG E   1       4.380   1.593  27.275  1.00  0.00           N\n',
    'ATOM     13  C8   DG E   1       4.890   0.325  27.422  1.00  0.00           C\n',
    'ATOM     14  N7   DG E   1       3.966  -0.596  27.466  1.00  0.00           N\n',
    'ATOM     15  C5   DG E   1       2.773   0.104  27.341  1.00  0.00           C\n',
    'ATOM     16  C6   DG E   1       1.433  -0.356  27.320  1.00  0.00           C\n',
    'ATOM     17  O6   DG E   1       1.017  -1.518  27.412  1.00  0.00           O\n',
    'ATOM     18  N1   DG E   1       0.531   0.692  27.176  1.00  0.00           N\n',
    'ATOM     19  C2   DG E   1       0.876   2.017  27.064  1.00  0.00           C\n',
    'ATOM     20  N2   DG E   1      -0.141   2.882  26.931  1.00  0.00           N\n',
    'ATOM     21  N3   DG E   1       2.121   2.461  27.082  1.00  0.00           N\n',
    'ATOM     22  C4   DG E   1       3.014   1.458  27.222  1.00  0.00           C\n',
    'ATOM     23  P    DC E   2       7.253   6.142  29.028  1.00  0.00           P\n',
    'ATOM     24  OP1  DC E   2       8.181   7.137  28.445  1.00  0.00           O\n',
    'ATOM     25  OP2  DC E   2       7.747   5.166  30.025  1.00  0.00           O\n',
    "ATOM     26  O5'  DC E   2       6.022   6.897  29.655  1.00  0.00           O\n",
    "ATOM     27  C5'  DC E   2       4.988   7.217  28.797  1.00  0.00           C\n",
    "ATOM     28  C4'  DC E   2       3.648   7.002  29.439  1.00  0.00           C\n",
    "ATOM     29  O4'  DC E   2       3.156   5.635  29.327  1.00  0.00           O\n",
    "ATOM     30  C3'  DC E   2       3.515   7.322  30.929  1.00  0.00           C\n",
    "ATOM     31  O3'  DC E   2       2.411   8.106  31.224  1.00  0.00           O\n",
    "ATOM     32  C2'  DC E   2       3.513   5.958  31.618  1.00  0.00           C\n",
    "ATOM     33  C1'  DC E   2       2.644   5.241  30.591  1.00  0.00           C\n",
    'ATOM     34  N1   DC E   2       2.717   3.776  30.675  1.00  0.00           N\n',
    'ATOM     35  C2   DC E   2       1.525   3.052  30.600  1.00  0.00           C\n',
    'ATOM     36  O2   DC E   2       0.456   3.667  30.468  1.00  0.00           O\n',
    'ATOM     37  N3   DC E   2       1.565   1.701  30.675  1.00  0.00           N\n',
    'ATOM     38  C4   DC E   2       2.736   1.077  30.818  1.00  0.00           C\n',
    'ATOM     39  N4   DC E   2       2.728  -0.258  30.888  1.00  0.00           N\n',
    'ATOM     40  C5   DC E   2       3.966   1.790  30.897  1.00  0.00           C\n',
    'ATOM     41  C6   DC E   2       3.911   3.126  30.822  1.00  0.00           C\n',
    'ATOM     42  P    DG E   3       2.367   8.967  32.582  1.00  0.00           P\n',
    'ATOM     43  OP1  DG E   3       1.107   9.742  32.585  1.00  0.00           O\n',
    'ATOM     44  OP2  DG E   3       3.664   9.670  32.705  1.00  0.00           O\n',
    "ATOM     45  O5'  DG E   3       2.302   7.865  33.743  1.00  0.00           O\n",
    "ATOM     46  C5'  DG E   3       1.083   7.625  34.445  1.00  0.00           C\n",
    "ATOM     47  C4'  DG E   3       0.203   6.640  33.694  1.00  0.00           C\n",
    "ATOM     48  O4'  DG E   3       0.966   5.453  33.337  1.00  0.00           O\n",
    "ATOM     49  C3'  DG E   3      -1.004   6.119  34.477  1.00  0.00           C\n",
    "ATOM     50  O3'  DG E   3      -2.095   5.966  33.591  1.00  0.00           O\n",
    "ATOM     51  C2'  DG E   3      -0.510   4.763  34.965  1.00  0.00           C\n",
    "ATOM     52  C1'  DG E   3       0.228   4.317  33.727  1.00  0.00           C\n",
    'ATOM     53  N9   DG E   3       1.137   3.202  33.945  1.00  0.00           N\n',
    'ATOM     54  C8   DG E   3       2.482   3.270  34.208  1.00  0.00           C\n',
    'ATOM     55  N7   DG E   3       3.041   2.101  34.352  1.00  0.00           N\n',
    'ATOM     56  C5   DG E   3       1.999   1.200  34.171  1.00  0.00           C\n',
    'ATOM     57  C6   DG E   3       1.999  -0.215  34.211  1.00  0.00           C\n',
    'ATOM     58  O6   DG E   3       2.952  -0.977  34.421  1.00  0.00           O\n',
    'ATOM     59  N1   DG E   3       0.730  -0.737  33.974  1.00  0.00           N\n',
    'ATOM     60  C2   DG E   3      -0.396   0.013  33.729  1.00  0.00           C\n',
    'ATOM     61  N2   DG E   3      -1.532  -0.671  33.523  1.00  0.00           N\n',
    'ATOM     62  N3   DG E   3      -0.409   1.339  33.689  1.00  0.00           N\n',
    'ATOM     63  C4   DG E   3       0.820   1.862  33.918  1.00  0.00           C\n',
    'ATOM     64  P    DG E   4      -3.557   5.631  34.157  1.00  0.00           P\n',
    'ATOM     65  OP1  DG E   4      -3.394   5.088  35.526  1.00  0.00           O\n',
    'ATOM     66  OP2  DG E   4      -4.266   4.832  33.132  1.00  0.00           O\n',
    "ATOM     67  O5'  DG E   4      -4.273   7.062  34.240  1.00  0.00           O\n",
    "ATOM     68  C5'  DG E   4      -5.687   7.152  34.460  1.00  0.00           C\n",
    "ATOM     69  C4'  DG E   4      -6.091   6.403  35.720  1.00  0.00           C\n",
    "ATOM     70  O4'  DG E   4      -7.225   7.047  36.352  1.00  0.00           O\n",
    "ATOM     71  C3'  DG E   4      -6.511   4.941  35.519  1.00  0.00           C\n",
    "ATOM     72  O3'  DG E   4      -6.008   4.170  36.594  1.00  0.00           O\n",
    "ATOM     73  C2'  DG E   4      -8.026   5.023  35.571  1.00  0.00           C\n",
    "ATOM     74  C1'  DG E   4      -8.173   6.053  36.667  1.00  0.00           C\n",
    'ATOM     75  N9   DG E   4      -9.485   6.665  36.738  1.00  0.00           N\n',
    'ATOM     76  C8   DG E   4     -10.593   6.352  35.992  1.00  0.00           C\n',
    'ATOM     77  N7   DG E   4     -11.635   7.077  36.288  1.00  0.00           N\n',
    'ATOM     78  C5   DG E   4     -11.184   7.921  37.297  1.00  0.00           C\n',
    'ATOM     79  C6   DG E   4     -11.864   8.930  38.016  1.00  0.00           C\n',
    'ATOM     80  O6   DG E   4     -13.043   9.293  37.902  1.00  0.00           O\n',
    'ATOM     81  N1   DG E   4     -11.034   9.545  38.952  1.00  0.00           N\n',
    'ATOM     82  C2   DG E   4      -9.715   9.223  39.163  1.00  0.00           C\n',
    'ATOM     83  N2   DG E   4      -9.075   9.926  40.109  1.00  0.00           N\n',
    'ATOM     84  N3   DG E   4      -9.068   8.281  38.496  1.00  0.00           N\n',
    'ATOM     85  C4   DG E   4      -9.863   7.675  37.583  1.00  0.00           C\n',
    'ATOM     86  P    DG E   5      -6.047   2.567  36.533  1.00  0.00           P\n',
    'ATOM     87  OP1  DG E   5      -5.991   2.182  35.105  1.00  0.00           O\n',
    'ATOM     88  OP2  DG E   5      -7.174   2.110  37.377  1.00  0.00           O\n',
    "ATOM     89  O5'  DG E   5      -4.683   2.123  37.245  1.00  0.00           O\n",
    "ATOM     90  C5'  DG E   5      -3.509   1.922  36.464  1.00  0.00           C\n",
    "ATOM     91  C4'  DG E   5      -2.314   1.566  37.340  1.00  0.00           C\n",
    "ATOM     92  O4'  DG E   5      -1.442   2.719  37.490  1.00  0.00           O\n",
    "ATOM     93  C3'  DG E   5      -2.642   1.125  38.770  1.00  0.00           C\n",
    "ATOM     94  O3'  DG E   5      -1.715   0.138  39.185  1.00  0.00           O\n",
    "ATOM     95  C2'  DG E   5      -2.420   2.403  39.557  1.00  0.00           C\n",
    "ATOM     96  C1'  DG E   5      -1.165   2.881  38.863  1.00  0.00           C\n",
    'ATOM     97  N9   DG E   5      -0.829   4.269  39.121  1.00  0.00           N\n',
    'ATOM     98  C8   DG E   5      -1.558   5.181  39.842  1.00  0.00           C\n',
    'ATOM     99  N7   DG E   5      -0.995   6.355  39.908  1.00  0.00           N\n',
    'ATOM    100  C5   DG E   5       0.184   6.207  39.187  1.00  0.00           C\n',
    'ATOM    101  C6   DG E   5       1.209   7.141  38.909  1.00  0.00           C\n',
    'ATOM    102  O6   DG E   5       1.281   8.327  39.258  1.00  0.00           O\n',
    'ATOM    103  N1   DG E   5       2.227   6.577  38.143  1.00  0.00           N\n',
    'ATOM    104  C2   DG E   5       2.250   5.275  37.702  1.00  0.00           C\n',
    'ATOM    105  N2   DG E   5       3.316   4.913  36.974  1.00  0.00           N\n',
    'ATOM    106  N3   DG E   5       1.296   4.391  37.955  1.00  0.00           N\n',
    'ATOM    107  C4   DG E   5       0.299   4.926  38.700  1.00  0.00           C\n',
    'ATOM    108  P    DG E   6      -2.235  -1.121  40.025  1.00  0.00           P\n',
    'ATOM    109  OP1  DG E   6      -2.199  -2.438  39.234  1.00  0.00           O\n',
    'ATOM    110  OP2  DG E   6      -3.707  -0.927  40.422  1.00  0.00           O\n',
    "ATOM    111  O5'  DG E   6      -1.311  -1.234  41.383  1.00  0.00           O\n",
    "ATOM    112  C5'  DG E   6      -0.732  -2.461  41.904  1.00  0.00           C\n",
    "ATOM    113  C4'  DG E   6      -0.350  -2.390  43.419  1.00  0.00           C\n",
    "ATOM    114  O4'  DG E   6      -0.476  -3.734  43.966  1.00  0.00           O\n",
    "ATOM    115  C3'  DG E   6       1.078  -1.919  43.829  1.00  0.00           C\n",
    "ATOM    116  O3'  DG E   6       1.059  -1.378  45.189  1.00  0.00           O\n",
    "ATOM    117  C2'  DG E   6       1.828  -3.255  43.810  1.00  0.00           C\n",
    "ATOM    118  C1'  DG E   6       0.805  -4.207  44.444  1.00  0.00           C\n",
    'ATOM    119  N9   DG E   6       0.938  -5.600  44.028  1.00  0.00           N\n',
    'ATOM    120  C8   DG E   6       0.375  -6.193  42.922  1.00  0.00           C\n',
    'ATOM    121  N7   DG E   6       0.674  -7.459  42.812  1.00  0.00           N\n',
    'ATOM    122  C5   DG E   6       1.483  -7.720  43.910  1.00  0.00           C\n',
    'ATOM    123  C6   DG E   6       2.114  -8.920  44.325  1.00  0.00           C\n',
    'ATOM    124  O6   DG E   6       2.081 -10.033  43.785  1.00  0.00           O\n',
    'ATOM    125  N1   DG E   6       2.842  -8.740  45.496  1.00  0.00           N\n',
    'ATOM    126  C2   DG E   6       2.950  -7.554  46.181  1.00  0.00           C\n',
    'ATOM    127  N2   DG E   6       3.700  -7.582  47.294  1.00  0.00           N\n',
    'ATOM    128  N3   DG E   6       2.369  -6.428  45.805  1.00  0.00           N\n',
    'ATOM    129  C4   DG E   6       1.655  -6.580  44.670  1.00  0.00           C\n',
    'ATOM    130  P    DG E   7       2.326  -0.917  46.126  1.00  0.00           P\n',
    'ATOM    131  OP1  DG E   7       3.283  -0.012  45.334  1.00  0.00           O\n',
    'ATOM    132  OP2  DG E   7       3.164  -2.143  46.522  1.00  0.00           O\n',
    "ATOM    133  O5'  DG E   7       1.794  -0.153  47.483  1.00  0.00           O\n",
    "ATOM    134  C5'  DG E   7       2.322   1.096  48.004  1.00  0.00           C\n",
    "ATOM    135  C4'  DG E   7       2.014   1.334  49.519  1.00  0.00           C\n",
    "ATOM    136  O4'  DG E   7       3.100   2.137  50.066  1.00  0.00           O\n",
    "ATOM    137  C3'  DG E   7       0.710   2.083  49.929  1.00  0.00           C\n",
    "ATOM    138  O3'  DG E   7       0.321   1.708  51.289  1.00  0.00           O\n",
    "ATOM    139  C2'  DG E   7       1.205   3.533  49.910  1.00  0.00           C\n",
    "ATOM    140  C1'  DG E   7       2.597   3.406  50.543  1.00  0.00           C\n",
    'ATOM    141  N9   DG E   7       3.547   4.434  50.128  1.00  0.00           N\n',
    'ATOM    142  C8   DG E   7       4.364   4.410  49.022  1.00  0.00           C\n',
    'ATOM    143  N7   DG E   7       5.107   5.477  48.912  1.00  0.00           N\n',
    'ATOM    144  C5   DG E   7       4.762   6.255  50.010  1.00  0.00           C\n',
    'ATOM    145  C6   DG E   7       5.236   7.525  50.425  1.00  0.00           C\n',
    'ATOM    146  O6   DG E   7       6.086   8.243  49.885  1.00  0.00           O\n',
    'ATOM    147  N1   DG E   7       4.615   7.947  51.596  1.00  0.00           N\n',
    'ATOM    148  C2   DG E   7       3.660   7.237  52.281  1.00  0.00           C\n',
    'ATOM    149  N2   DG E   7       3.181   7.814  53.394  1.00  0.00           N\n',
    'ATOM    150  N3   DG E   7       3.208   6.052  51.905  1.00  0.00           N\n',
    'ATOM    151  C4   DG E   7       3.798   5.623  50.770  1.00  0.00           C\n',
    'TER     152       DG E   7\n',
    'ATOM    153  P    DG F   1       2.495  -0.164  55.275  1.00  0.00           P\n',
    'ATOM    154  OP1  DG F   1       2.974  -1.391  56.066  1.00  0.00           O\n',
    'ATOM    155  OP2  DG F   1       3.776   0.587  54.879  1.00  0.00           O\n',
    "ATOM    156  O5'  DG F   1       1.688  -0.627  53.917  1.00  0.00           O\n",
    "ATOM    157  C5'  DG F   1       1.632  -1.983  53.397  1.00  0.00           C\n",
    "ATOM    158  C4'  DG F   1       1.251  -2.066  51.882  1.00  0.00           C\n",
    "ATOM    159  O4'  DG F   1       1.890  -3.256  51.334  1.00  0.00           O\n",
    "ATOM    160  C3'  DG F   1      -0.248  -2.187  51.471  1.00  0.00           C\n",
    "ATOM    161  O3'  DG F   1      -0.440  -1.682  50.111  1.00  0.00           O\n",
    "ATOM    162  C2'  DG F   1      -0.419  -3.709  51.490  1.00  0.00           C\n",
    "ATOM    163  C1'  DG F   1       0.894  -4.189  50.858  1.00  0.00           C\n",
    'ATOM    164  N9   DG F   1       1.313  -5.524  51.273  1.00  0.00           N\n',
    'ATOM    165  C8   DG F   1       2.063  -5.852  52.379  1.00  0.00           C\n',
    'ATOM    166  N7   DG F   1       2.279  -7.134  52.489  1.00  0.00           N\n',
    'ATOM    167  C5   DG F   1       1.635  -7.689  51.391  1.00  0.00           C\n',
    'ATOM    168  C6   DG F   1       1.521  -9.040  50.976  1.00  0.00           C\n',
    'ATOM    169  O6   DG F   1       1.983 -10.053  51.516  1.00  0.00           O\n',
    'ATOM    170  N1   DG F   1       0.779  -9.157  49.805  1.00  0.00           N\n',
    'ATOM    171  C2   DG F   1       0.219  -8.107  49.120  1.00  0.00           C\n',
    'ATOM    172  N2   DG F   1      -0.461  -8.424  48.007  1.00  0.00           N\n',
    'ATOM    173  N3   DG F   1       0.317  -6.843  49.496  1.00  0.00           N\n',
    'ATOM    174  C4   DG F   1       1.033  -6.706  50.631  1.00  0.00           C\n',
    'ATOM    175  P    DG F   2      -1.786  -1.749  49.175  1.00  0.00           P\n',
    'ATOM    176  OP1  DG F   2      -3.020  -1.287  49.966  1.00  0.00           O\n',
    'ATOM    177  OP2  DG F   2      -2.081  -3.205  48.778  1.00  0.00           O\n',
    "ATOM    178  O5'  DG F   2      -1.593  -0.839  47.817  1.00  0.00           O\n",
    "ATOM    179  C5'  DG F   2      -2.566   0.107  47.296  1.00  0.00           C\n",
    "ATOM    180  C4'  DG F   2      -2.374   0.446  45.781  1.00  0.00           C\n",
    "ATOM    181  O4'  DG F   2      -3.686   0.763  45.234  1.00  0.00           O\n",
    "ATOM    182  C3'  DG F   2      -1.464   1.643  45.371  1.00  0.00           C\n",
    "ATOM    183  O3'  DG F   2      -0.960   1.449  44.011  1.00  0.00           O\n",
    "ATOM    184  C2'  DG F   2      -2.484   2.786  45.390  1.00  0.00           C\n",
    "ATOM    185  C1'  DG F   2      -3.717   2.128  44.757  1.00  0.00           C\n",
    'ATOM    186  N9   DG F   2      -4.992   2.707  45.173  1.00  0.00           N\n',
    'ATOM    187  C8   DG F   2      -5.735   2.367  46.279  1.00  0.00           C\n',
    'ATOM    188  N7   DG F   2      -6.835   3.061  46.389  1.00  0.00           N\n',
    'ATOM    189  C5   DG F   2      -6.819   3.911  45.291  1.00  0.00           C\n',
    'ATOM    190  C6   DG F   2      -7.749   4.897  44.876  1.00  0.00           C\n',
    'ATOM    191  O6   DG F   2      -8.813   5.228  45.416  1.00  0.00           O\n',
    'ATOM    192  N1   DG F   2      -7.342   5.528  43.705  1.00  0.00           N\n',
    'ATOM    193  C2   DG F   2      -6.186   5.245  43.020  1.00  0.00           C\n',
    'ATOM    194  N2   DG F   2      -5.969   5.963  41.907  1.00  0.00           N\n',
    'ATOM    195  N3   DG F   2      -5.309   4.329  43.396  1.00  0.00           N\n',
    'ATOM    196  C4   DG F   2      -5.686   3.704  44.531  1.00  0.00           C\n',
    'ATOM    197  P    DG F   3      -0.127   2.618  43.285  1.00  0.00           P\n',
    'ATOM    198  OP1  DG F   3       0.618   3.358  44.326  1.00  0.00           O\n',
    'ATOM    199  OP2  DG F   3      -1.063   3.342  42.396  1.00  0.00           O\n',
    "ATOM    200  O5'  DG F   3       0.928   1.836  42.367  1.00  0.00           O\n",
    "ATOM    201  C5'  DG F   3       2.182   2.445  42.043  1.00  0.00           C\n",
    "ATOM    202  C4'  DG F   3       2.519   2.242  40.575  1.00  0.00           C\n",
    "ATOM    203  O4'  DG F   3       3.301   3.364  40.079  1.00  0.00           O\n",
    "ATOM    204  C3'  DG F   3       3.350   0.994  40.256  1.00  0.00           C\n",
    "ATOM    205  O3'  DG F   3       2.898   0.432  39.037  1.00  0.00           O\n",
    "ATOM    206  C2'  DG F   3       4.753   1.565  40.098  1.00  0.00           C\n",
    "ATOM    207  C1'  DG F   3       4.413   2.854  39.380  1.00  0.00           C\n",
    'ATOM    208  N9   DG F   3       5.480   3.845  39.386  1.00  0.00           N\n',
    'ATOM    209  C8   DG F   3       5.638   4.888  40.264  1.00  0.00           C\n',
    'ATOM    210  N7   DG F   3       6.688   5.620  40.015  1.00  0.00           N\n',
    'ATOM    211  C5   DG F   3       7.260   5.026  38.897  1.00  0.00           C\n',
    'ATOM    212  C6   DG F   3       8.422   5.379  38.169  1.00  0.00           C\n',
    'ATOM    213  O6   DG F   3       9.202   6.319  38.376  1.00  0.00           O\n',
    'ATOM    214  N1   DG F   3       8.646   4.511  37.103  1.00  0.00           N\n',
    'ATOM    215  C2   DG F   3       7.849   3.439  36.781  1.00  0.00           C\n',
    'ATOM    216  N2   DG F   3       8.226   2.716  35.716  1.00  0.00           N\n',
    'ATOM    217  N3   DG F   3       6.758   3.097  37.454  1.00  0.00           N\n',
    'ATOM    218  C4   DG F   3       6.526   3.933  38.496  1.00  0.00           C\n',
    'ATOM    219  P    DG F   4       3.706  -0.775  38.352  1.00  0.00           P\n',
    'ATOM    220  OP1  DG F   4       4.158  -1.667  39.442  1.00  0.00           O\n',
    'ATOM    221  OP2  DG F   4       4.692  -0.207  37.405  1.00  0.00           O\n',
    "ATOM    222  O5'  DG F   4       2.595  -1.542  37.497  1.00  0.00           O\n",
    "ATOM    223  C5'  DG F   4       1.810  -2.559  38.107  1.00  0.00           C\n",
    "ATOM    224  C4'  DG F   4       0.928  -3.251  37.082  1.00  0.00           C\n",
    "ATOM    225  O4'  DG F   4       1.735  -3.720  35.965  1.00  0.00           O\n",
    "ATOM    226  C3'  DG F   4       0.182  -4.488  37.591  1.00  0.00           C\n",
    "ATOM    227  O3'  DG F   4      -1.102  -4.535  36.996  1.00  0.00           O\n",
    "ATOM    228  C2'  DG F   4       1.054  -5.627  37.080  1.00  0.00           C\n",
    "ATOM    229  C1'  DG F   4       1.409  -5.069  35.720  1.00  0.00           C\n",
    'ATOM    230  N9   DG F   4       2.538  -5.725  35.075  1.00  0.00           N\n',
    'ATOM    231  C8   DG F   4       3.707  -6.141  35.663  1.00  0.00           C\n',
    'ATOM    232  N7   DG F   4       4.537  -6.698  34.826  1.00  0.00           N\n',
    'ATOM    233  C5   DG F   4       3.874  -6.648  33.605  1.00  0.00           C\n',
    'ATOM    234  C6   DG F   4       4.278  -7.100  32.326  1.00  0.00           C\n',
    'ATOM    235  O6   DG F   4       5.338  -7.655  32.007  1.00  0.00           O\n',
    'ATOM    236  N1   DG F   4       3.304  -6.853  31.361  1.00  0.00           N\n',
    'ATOM    237  C2   DG F   4       2.095  -6.246  31.601  1.00  0.00           C\n',
    'ATOM    238  N2   DG F   4       1.285  -6.093  30.544  1.00  0.00           N\n',
    'ATOM    239  N3   DG F   4       1.705  -5.818  32.794  1.00  0.00           N\n',
    'ATOM    240  C4   DG F   4       2.643  -6.051  33.744  1.00  0.00           C\n',
    'ATOM    241  P    DG F   5      -2.089  -5.767  37.289  1.00  0.00           P\n',
    'ATOM    242  OP1  DG F   5      -2.035  -6.027  38.744  1.00  0.00           O\n',
    'ATOM    243  OP2  DG F   5      -1.781  -6.854  36.332  1.00  0.00           O\n',
    "ATOM    244  O5'  DG F   5      -3.534  -5.191  36.921  1.00  0.00           O\n",
    "ATOM    245  C5'  DG F   5      -4.535  -6.055  36.390  1.00  0.00           C\n",
    "ATOM    246  C4'  DG F   5      -4.756  -5.775  34.912  1.00  0.00           C\n",
    "ATOM    247  O4'  DG F   5      -3.731  -6.444  34.121  1.00  0.00           O\n",
    "ATOM    248  C3'  DG F   5      -6.095  -6.256  34.339  1.00  0.00           C\n",
    "ATOM    249  O3'  DG F   5      -6.581  -5.279  33.387  1.00  0.00           O\n",
    "ATOM    250  C2'  DG F   5      -5.701  -7.567  33.665  1.00  0.00           C\n",
    "ATOM    251  C1'  DG F   5      -4.367  -7.148  33.080  1.00  0.00           C\n",
    'ATOM    252  N9   DG F   5      -3.518  -8.251  32.650  1.00  0.00           N\n',
    'ATOM    253  C8   DG F   5      -2.779  -9.091  33.445  1.00  0.00           C\n',
    'ATOM    254  N7   DG F   5      -2.106  -9.984  32.772  1.00  0.00           N\n',
    'ATOM    255  C5   DG F   5      -2.417  -9.715  31.444  1.00  0.00           C\n',
    'ATOM    256  C6   DG F   5      -1.985 -10.355  30.257  1.00  0.00           C\n',
    'ATOM    257  O6   DG F   5      -1.216 -11.319  30.141  1.00  0.00           O\n',
    'ATOM    258  N1   DG F   5      -2.541  -9.767  29.123  1.00  0.00           N\n',
    'ATOM    259  C2   DG F   5      -3.403  -8.697  29.135  1.00  0.00           C\n',
    'ATOM    260  N2   DG F   5      -3.834  -8.268  27.939  1.00  0.00           N\n',
    'ATOM    261  N3   DG F   5      -3.816  -8.089  30.239  1.00  0.00           N\n',
    'ATOM    262  C4   DG F   5      -3.284  -8.650  31.353  1.00  0.00           C\n',
    'ATOM    263  P    DG F   6      -7.576  -5.685  32.185  1.00  0.00           P\n',
    'ATOM    264  OP1  DG F   6      -8.748  -6.381  32.763  1.00  0.00           O\n',
    'ATOM    265  OP2  DG F   6      -6.730  -6.370  31.181  1.00  0.00           O\n',
    "ATOM    266  O5'  DG F   6      -8.048  -4.315  31.570  1.00  0.00           O\n",
    "ATOM    267  C5'  DG F   6      -8.097  -3.236  32.431  1.00  0.00           C\n",
    "ATOM    268  C4'  DG F   6      -7.596  -1.973  31.779  1.00  0.00           C\n",
    "ATOM    269  O4'  DG F   6      -6.153  -1.801  31.879  1.00  0.00           O\n",
    "ATOM    270  C3'  DG F   6      -7.888  -1.779  30.290  1.00  0.00           C\n",
    "ATOM    271  O3'  DG F   6      -8.409  -0.543  29.993  1.00  0.00           O\n",
    "ATOM    272  C2'  DG F   6      -6.563  -2.088  29.593  1.00  0.00           C\n",
    "ATOM    273  C1'  DG F   6      -5.662  -1.397  30.610  1.00  0.00           C\n",
    'ATOM    274  N9   DG F   6      -4.261  -1.796  30.518  1.00  0.00           N\n',
    'ATOM    275  C8   DG F   6      -3.770  -3.073  30.375  1.00  0.00           C\n',
    'ATOM    276  N7   DG F   6      -2.467  -3.117  30.322  1.00  0.00           N\n',
    'ATOM    277  C5   DG F   6      -2.071  -1.790  30.438  1.00  0.00           C\n',
    'ATOM    278  C6   DG F   6      -0.776  -1.214  30.446  1.00  0.00           C\n',
    'ATOM    279  O6   DG F   6       0.320  -1.781  30.350  1.00  0.00           O\n',
    'ATOM    280  N1   DG F   6      -0.831   0.169  30.584  1.00  0.00           N\n',
    'ATOM    281  C2   DG F   6      -1.985   0.904  30.699  1.00  0.00           C\n',
    'ATOM    282  N2   DG F   6      -1.831   2.231  30.824  1.00  0.00           N\n',
    'ATOM    283  N3   DG F   6      -3.199   0.379  30.692  1.00  0.00           N\n',
    'ATOM    284  C4   DG F   6      -3.169  -0.963  30.559  1.00  0.00           C\n',
    'ATOM    285  P    DC F   7      -9.462  -0.430  28.785  1.00  0.00           P\n',
    'ATOM    286  OP1  DC F   7     -10.823  -0.345  29.362  1.00  0.00           O\n',
    'ATOM    287  OP2  DC F   7      -9.149  -1.473  27.781  1.00  0.00           O\n',
    "ATOM    288  O5'  DC F   7      -9.081   0.968  28.170  1.00  0.00           O\n",
    "ATOM    289  C5'  DC F   7      -8.514   1.887  29.032  1.00  0.00           C\n",
    "ATOM    290  C4'  DC F   7      -7.394   2.638  28.378  1.00  0.00           C\n",
    "ATOM    291  O4'  DC F   7      -6.101   1.974  28.479  1.00  0.00           O\n",
    "ATOM    292  C3'  DC F   7      -7.527   2.962  26.889  1.00  0.00           C\n",
    "ATOM    293  O3'  DC F   7      -7.255   4.288  26.592  1.00  0.00           O\n",
    "ATOM    294  C2'  DC F   7      -6.602   1.965  26.192  1.00  0.00           C\n",
    "ATOM    295  C1'  DC F   7      -5.468   2.035  27.209  1.00  0.00           C\n",
    'ATOM    296  N1   DC F   7      -4.523   0.912  27.118  1.00  0.00           N\n',
    'ATOM    297  C2   DC F   7      -3.157   1.195  27.180  1.00  0.00           C\n',
    'ATOM    298  O2   DC F   7      -2.795   2.374  27.308  1.00  0.00           O\n',
    'ATOM    299  N3   DC F   7      -2.266   0.180  27.099  1.00  0.00           N\n',
    'ATOM    300  C4   DC F   7      -2.699  -1.075  26.960  1.00  0.00           C\n',
    'ATOM    301  N4   DC F   7      -1.784  -2.046  26.883  1.00  0.00           N\n',
    'ATOM    302  C5   DC F   7      -4.085  -1.392  26.894  1.00  0.00           C\n',
    'ATOM    303  C6   DC F   7      -4.955  -0.377  26.975  1.00  0.00           C\n',
    'TER     304       DC F   7\n',
    'END\n',
]


def core_records():
    """Return a *new* list of edit_pdb_atom records for the B-P junction core."""
    from .edit_pdb_atom import pdb_atom_record, pdb_ter_record

    recs = []
    for line in PDB_LINES:
        if line.startswith(("ATOM", "HETATM")):
            recs.append(pdb_atom_record(line))
        elif line.startswith("TER"):
            recs.append(pdb_ter_record(line))
    return recs
