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
  * Result: phosphate steps, in the +Z direction, strand I (E) 35, 45, 36 and
    51 deg from the last B pair to the first P-form level; strand II (F) 33,
    71, 131 and 133 deg. Strand I's phosphates move inward 9.5 -> 9.3 -> 7.0
    -> 6.4 -> 2.5 A, strand II's 9.5 -> 6.9 -> 4.2 -> 2.6 -> 2.5 A. Puckers P
    = 153-155 deg; chi -58 to -158 deg.
  * Finish: phenix.geometry_minimization moved only the transition residues,
    with every atom restrained to its start (sigma 0.08 A; largest shift
    0.12 A), to bring the closing O3'-P-O5' angles in line.
  * Checks on the transition residues, against Phenix restraints: bonds RMSD
    0.005 A (largest 0.8 sigma), angles RMSD 0.8 deg (largest 2.1 sigma),
    chirality as built; the eight O3'-P linkages into, within and out of the
    transition are 1.599-1.610 A. MolProbity clashscore is 0 for 44 base
    assignments (all A, all C, all G, all T and 40 random ones) of the
    transition, P and B residues around the junction.

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
nucleotides there. The G/C bases stored below are placeholders.

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
# P-form levels that must separate two junctions on one P-DNA. This core's
# transitions wind forward and stay clear of each other with a single level
# between them (closest contact 2.7 A, a Watson-Crick hydrogen bond in the
# B-DNA, over random builds with 1, 2 and 3 levels kept).
MIN_KEPT_P_LEVELS_BETWEEN_JUNCTIONS = 1

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
    'ATOM     42  P    DG E   3       2.386   9.000  32.561  1.00  0.00           P\n',
    'ATOM     43  OP1  DG E   3       1.104   9.739  32.583  1.00  0.00           O\n',
    'ATOM     44  OP2  DG E   3       3.664   9.743  32.624  1.00  0.00           O\n',
    "ATOM     45  O5'  DG E   3       2.383   7.930  33.754  1.00  0.00           O\n",
    "ATOM     46  C5'  DG E   3       1.159   7.565  34.396  1.00  0.00           C\n",
    "ATOM     47  C4'  DG E   3       0.369   6.581  33.549  1.00  0.00           C\n",
    "ATOM     48  O4'  DG E   3       1.039   5.289  33.525  1.00  0.00           O\n",
    "ATOM     49  C3'  DG E   3      -1.053   6.289  34.038  1.00  0.00           C\n",
    "ATOM     50  O3'  DG E   3      -1.907   6.142  32.917  1.00  0.00           O\n",
    "ATOM     51  C2'  DG E   3      -0.883   4.961  34.764  1.00  0.00           C\n",
    "ATOM     52  C1'  DG E   3       0.088   4.291  33.818  1.00  0.00           C\n",
    'ATOM     53  N9   DG E   3       0.766   3.128  34.375  1.00  0.00           N\n',
    'ATOM     54  C8   DG E   3       2.002   3.097  34.971  1.00  0.00           C\n',
    'ATOM     55  N7   DG E   3       2.354   1.906  35.369  1.00  0.00           N\n',
    'ATOM     56  C5   DG E   3       1.284   1.096  35.011  1.00  0.00           C\n',
    'ATOM     57  C6   DG E   3       1.093  -0.296  35.187  1.00  0.00           C\n',
    'ATOM     58  O6   DG E   3       1.859  -1.115  35.712  1.00  0.00           O\n',
    'ATOM     59  N1   DG E   3      -0.134  -0.716  34.678  1.00  0.00           N\n',
    'ATOM     60  C2   DG E   3      -1.058   0.104  34.077  1.00  0.00           C\n',
    'ATOM     61  N2   DG E   3      -2.185  -0.485  33.649  1.00  0.00           N\n',
    'ATOM     62  N3   DG E   3      -0.892   1.409  33.906  1.00  0.00           N\n',
    'ATOM     63  C4   DG E   3       0.299   1.834  34.396  1.00  0.00           C\n',
    'ATOM     64  P    DG E   4      -3.498   6.078  33.113  1.00  0.00           P\n',
    'ATOM     65  OP1  DG E   4      -3.822   4.751  33.683  1.00  0.00           O\n',
    'ATOM     66  OP2  DG E   4      -4.112   6.484  31.829  1.00  0.00           O\n',
    "ATOM     67  O5'  DG E   4      -3.796   7.226  34.188  1.00  0.00           O\n",
    "ATOM     68  C5'  DG E   4      -5.100   7.367  34.752  1.00  0.00           C\n",
    "ATOM     69  C4'  DG E   4      -5.246   6.492  35.986  1.00  0.00           C\n",
    "ATOM     70  O4'  DG E   4      -5.991   7.194  37.017  1.00  0.00           O\n",
    "ATOM     71  C3'  DG E   4      -5.993   5.172  35.774  1.00  0.00           C\n",
    "ATOM     72  O3'  DG E   4      -5.378   4.169  36.555  1.00  0.00           O\n",
    "ATOM     73  C2'  DG E   4      -7.382   5.491  36.304  1.00  0.00           C\n",
    "ATOM     74  C1'  DG E   4      -6.994   6.329  37.501  1.00  0.00           C\n",
    'ATOM     75  N9   DG E   4      -8.076   7.125  38.053  1.00  0.00           N\n',
    'ATOM     76  C8   DG E   4      -9.396   7.116  37.674  1.00  0.00           C\n',
    'ATOM     77  N7   DG E   4     -10.136   7.940  38.362  1.00  0.00           N\n',
    'ATOM     78  C5   DG E   4      -9.249   8.527  39.256  1.00  0.00           C\n',
    'ATOM     79  C6   DG E   4      -9.470   9.502  40.257  1.00  0.00           C\n',
    'ATOM     80  O6   DG E   4     -10.532  10.060  40.564  1.00  0.00           O\n',
    'ATOM     81  N1   DG E   4      -8.296   9.818  40.937  1.00  0.00           N\n',
    'ATOM     82  C2   DG E   4      -7.065   9.262  40.682  1.00  0.00           C\n',
    'ATOM     83  N2   DG E   4      -6.048   9.693  41.442  1.00  0.00           N\n',
    'ATOM     84  N3   DG E   4      -6.846   8.350  39.748  1.00  0.00           N\n',
    'ATOM     85  C4   DG E   4      -7.979   8.032  39.077  1.00  0.00           C\n',
    'ATOM     86  P    DG E   5      -5.824   2.634  36.411  1.00  0.00           P\n',
    'ATOM     87  OP1  DG E   5      -5.925   2.344  34.964  1.00  0.00           O\n',
    'ATOM     88  OP2  DG E   5      -6.994   2.412  37.292  1.00  0.00           O\n',
    "ATOM     89  O5'  DG E   5      -4.583   1.826  37.016  1.00  0.00           O\n",
    "ATOM     90  C5'  DG E   5      -3.322   1.875  36.356  1.00  0.00           C\n",
    "ATOM     91  C4'  DG E   5      -2.184   1.543  37.311  1.00  0.00           C\n",
    "ATOM     92  O4'  DG E   5      -1.351   2.715  37.523  1.00  0.00           O\n",
    "ATOM     93  C3'  DG E   5      -2.604   1.094  38.713  1.00  0.00           C\n",
    "ATOM     94  O3'  DG E   5      -1.691   0.126  39.185  1.00  0.00           O\n",
    "ATOM     95  C2'  DG E   5      -2.462   2.373  39.517  1.00  0.00           C\n",
    "ATOM     96  C1'  DG E   5      -1.174   2.880  38.912  1.00  0.00           C\n",
    'ATOM     97  N9   DG E   5      -0.887   4.275  39.195  1.00  0.00           N\n',
    'ATOM     98  C8   DG E   5      -1.784   5.257  39.533  1.00  0.00           C\n',
    'ATOM     99  N7   DG E   5      -1.226   6.419  39.732  1.00  0.00           N\n',
    'ATOM    100  C5   DG E   5       0.127   6.188  39.514  1.00  0.00           C\n',
    'ATOM    101  C6   DG E   5       1.228   7.073  39.582  1.00  0.00           C\n',
    'ATOM    102  O6   DG E   5       1.226   8.280  39.860  1.00  0.00           O\n',
    'ATOM    103  N1   DG E   5       2.427   6.429  39.287  1.00  0.00           N\n',
    'ATOM    104  C2   DG E   5       2.546   5.098  38.967  1.00  0.00           C\n',
    'ATOM    105  N2   DG E   5       3.786   4.655  38.713  1.00  0.00           N\n',
    'ATOM    106  N3   DG E   5       1.523   4.258  38.899  1.00  0.00           N\n',
    'ATOM    107  C4   DG E   5       0.349   4.872  39.184  1.00  0.00           C\n',
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
    'ATOM    197  P    DG F   3      -0.081   2.594  43.300  1.00  0.00           P\n',
    'ATOM    198  OP1  DG F   3       0.985   2.996  44.244  1.00  0.00           O\n',
    'ATOM    199  OP2  DG F   3      -1.011   3.612  42.761  1.00  0.00           O\n',
    "ATOM    200  O5'  DG F   3       0.592   1.851  42.053  1.00  0.00           O\n",
    "ATOM    201  C5'  DG F   3       1.994   1.569  42.054  1.00  0.00           C\n",
    "ATOM    202  C4'  DG F   3       2.426   1.008  40.708  1.00  0.00           C\n",
    "ATOM    203  O4'  DG F   3       3.383   1.903  40.074  1.00  0.00           O\n",
    "ATOM    204  C3'  DG F   3       3.115  -0.359  40.752  1.00  0.00           C\n",
    "ATOM    205  O3'  DG F   3       2.706  -1.119  39.629  1.00  0.00           O\n",
    "ATOM    206  C2'  DG F   3       4.587   0.014  40.639  1.00  0.00           C\n",
    "ATOM    207  C1'  DG F   3       4.476   1.134  39.628  1.00  0.00           C\n",
    'ATOM    208  N9   DG F   3       5.656   1.982  39.539  1.00  0.00           N\n',
    'ATOM    209  C8   DG F   3       6.005   3.009  40.379  1.00  0.00           C\n',
    'ATOM    210  N7   DG F   3       7.120   3.596  40.044  1.00  0.00           N\n',
    'ATOM    211  C5   DG F   3       7.535   2.915  38.905  1.00  0.00           C\n',
    'ATOM    212  C6   DG F   3       8.679   3.103  38.094  1.00  0.00           C\n',
    'ATOM    213  O6   DG F   3       9.584   3.938  38.226  1.00  0.00           O\n',
    'ATOM    214  N1   DG F   3       8.717   2.196  37.038  1.00  0.00           N\n',
    'ATOM    215  C2   DG F   3       7.770   1.230  36.796  1.00  0.00           C\n',
    'ATOM    216  N2   DG F   3       7.979   0.447  35.728  1.00  0.00           N\n',
    'ATOM    217  N3   DG F   3       6.694   1.043  37.548  1.00  0.00           N\n',
    'ATOM    218  C4   DG F   3       6.642   1.919  38.582  1.00  0.00           C\n',
    'ATOM    219  P    DG F   4       3.249  -2.615  39.414  1.00  0.00           P\n',
    'ATOM    220  OP1  DG F   4       3.088  -3.328  40.700  1.00  0.00           O\n',
    'ATOM    221  OP2  DG F   4       4.580  -2.540  38.771  1.00  0.00           O\n',
    "ATOM    222  O5'  DG F   4       2.230  -3.233  38.347  1.00  0.00           O\n",
    "ATOM    223  C5'  DG F   4       2.344  -4.597  37.949  1.00  0.00           C\n",
    "ATOM    224  C4'  DG F   4       1.223  -4.962  36.990  1.00  0.00           C\n",
    "ATOM    225  O4'  DG F   4       1.769  -5.554  35.780  1.00  0.00           O\n",
    "ATOM    226  C3'  DG F   4       0.206  -5.979  37.519  1.00  0.00           C\n",
    "ATOM    227  O3'  DG F   4      -1.084  -5.622  37.056  1.00  0.00           O\n",
    "ATOM    228  C2'  DG F   4       0.673  -7.273  36.866  1.00  0.00           C\n",
    "ATOM    229  C1'  DG F   4       1.051  -6.734  35.504  1.00  0.00           C\n",
    'ATOM    230  N9   DG F   4       1.887  -7.624  34.713  1.00  0.00           N\n',
    'ATOM    231  C8   DG F   4       2.741  -8.595  35.174  1.00  0.00           C\n',
    'ATOM    232  N7   DG F   4       3.359  -9.235  34.220  1.00  0.00           N\n',
    'ATOM    233  C5   DG F   4       2.882  -8.650  33.053  1.00  0.00           C\n',
    'ATOM    234  C6   DG F   4       3.189  -8.926  31.700  1.00  0.00           C\n',
    'ATOM    235  O6   DG F   4       3.972  -9.773  31.248  1.00  0.00           O\n',
    'ATOM    236  N1   DG F   4       2.481  -8.099  30.831  1.00  0.00           N\n',
    'ATOM    237  C2   DG F   4       1.590  -7.129  31.221  1.00  0.00           C\n',
    'ATOM    238  N2   DG F   4       1.003  -6.430  30.238  1.00  0.00           N\n',
    'ATOM    239  N3   DG F   4       1.293  -6.860  32.484  1.00  0.00           N\n',
    'ATOM    240  C4   DG F   4       1.975  -7.658  33.342  1.00  0.00           C\n',
    'ATOM    241  P    DG F   5      -2.375  -6.458  37.515  1.00  0.00           P\n',
    'ATOM    242  OP1  DG F   5      -2.335  -6.539  38.991  1.00  0.00           O\n',
    'ATOM    243  OP2  DG F   5      -2.448  -7.696  36.706  1.00  0.00           O\n',
    "ATOM    244  O5'  DG F   5      -3.602  -5.529  37.080  1.00  0.00           O\n",
    "ATOM    245  C5'  DG F   5      -4.779  -6.115  36.533  1.00  0.00           C\n",
    "ATOM    246  C4'  DG F   5      -4.884  -5.826  35.044  1.00  0.00           C\n",
    "ATOM    247  O4'  DG F   5      -3.822  -6.517  34.326  1.00  0.00           O\n",
    "ATOM    248  C3'  DG F   5      -6.189  -6.276  34.377  1.00  0.00           C\n",
    "ATOM    249  O3'  DG F   5      -6.587  -5.295  33.394  1.00  0.00           O\n",
    "ATOM    250  C2'  DG F   5      -5.778  -7.594  33.727  1.00  0.00           C\n",
    "ATOM    251  C1'  DG F   5      -4.397  -7.207  33.241  1.00  0.00           C\n",
    'ATOM    252  N9   DG F   5      -3.543  -8.332  32.875  1.00  0.00           N\n',
    'ATOM    253  C8   DG F   5      -2.928  -9.220  33.723  1.00  0.00           C\n',
    'ATOM    254  N7   DG F   5      -2.215 -10.120  33.102  1.00  0.00           N\n',
    'ATOM    255  C5   DG F   5      -2.365  -9.805  31.756  1.00  0.00           C\n',
    'ATOM    256  C6   DG F   5      -1.825 -10.429  30.606  1.00  0.00           C\n',
    'ATOM    257  O6   DG F   5      -1.081 -11.418  30.546  1.00  0.00           O\n',
    'ATOM    258  N1   DG F   5      -2.228  -9.791  29.435  1.00  0.00           N\n',
    'ATOM    259  C2   DG F   5      -3.049  -8.689  29.383  1.00  0.00           C\n',
    'ATOM    260  N2   DG F   5      -3.327  -8.213  28.161  1.00  0.00           N\n',
    'ATOM    261  N3   DG F   5      -3.561  -8.095  30.453  1.00  0.00           N\n',
    'ATOM    262  C4   DG F   5      -3.178  -8.706  31.601  1.00  0.00           C\n',
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
