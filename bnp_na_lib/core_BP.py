"""
core_BP.py

B-P junction core for make_BP.py (bnp_na V13.22): two B-DNA base pairs, two
transition levels and two P-DNA levels.

There is no experimental structure of a B-DNA/P-DNA junction, so this core is a
model. It was built once, offline, and is bundled here so that a junction build
is only rigid fitting, trimming, base grafting and ligation, which is immediate.

  * B part: the top two base pairs of a bnp_na default B-DNA (DSSR rebuild at
    the default parameters, phenix.geometry_minimization with min_P_C5.params,
    phosphate regularization), aligned to +Z.
  * P part: bnp_na P-DNA (1PFI template, core_PDNA.py) on the same axis. Its
    first P-form level (the P interior level, E5/F2) is centered 10.0 A above
    the top B base pair and turned 250 deg about the axis relative to the
    bnp_na P-DNA frame.
  * Transition: on each strand the two nucleotides between them, which take
    the place of the P-DNA's first two levels, were built by torsion-space
    loop closure. Each is an ideal DSSR sugar-base unit (E3 C3'-endo from
    A-DNA; E4, F3 and F4 C2'-endo from B-DNA) joined through phosphate
    linkages with Phenix's ideal bond lengths and angles. Backbone torsions and
    glycosidic rotations were optimized so each strand closes onto the next
    helix, scoring contacts with riding sugar hydrogens and with the union of
    A, C, G and T at every base near the junction, so the core stays clash-free
    whatever bases are grafted onto it.
  * Finish: phenix.geometry_minimization moved only the transition residues,
    with every atom restrained to its start (sigma 0.08 A; largest shift
    0.11 A), to bring the closing O3'-P-O5' angle in line. The F3 base was then
    turned 15 deg about its glycosidic bond, which clears a thymine methyl that
    otherwise touches a purine two levels up.
  * Checks on the transition residues, against Phenix restraints: bonds RMSD
    0.005 A (largest 0.8 sigma), angles RMSD 1.0 deg (largest 2.3 sigma),
    chirality as built; the six O3'-P linkages into, within and out of the
    transition are 1.603-1.611 A. MolProbity clashscore is 0 for 88 base
    assignments (all A, all C, all G, all T and 84 random ones) of the
    transition, P and B residues around the junction.

Residues (chain E runs 5' B-DNA -> 3' P-DNA, chain F runs 5' P-DNA -> 3' B-DNA):

  E1/F6  B outer pair           E2/F5  B interior pair
  E3/F4  transition level 1     E4/F3  transition level 2
         (B side)                      (P side)
  E5/F2  P interior level       E6/F1  P outer level

At a junction, make_BP.py fits E1/F6 + E2/F5 onto the B-DNA's two terminal
base pairs and E5/F2 + E6/F1 onto the P-DNA's third and fourth levels from that
end. The core's B pairs and transition levels then replace the B-DNA's two
terminal pairs and the P-DNA's two terminal levels and take over their bases;
the core's P levels are dropped, as the P-DNA keeps its own nucleotides there.
The G/C bases stored below are placeholders.

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
B_END_BASEPAIRS = [(("E", 2), ("F", 5)), (("E", 1), ("F", 6))]
# Transition levels, B side -> P side. They take the place of the P-DNA's two
# terminal levels and take over their bases.
TRANSITION_LEVELS = [(("E", 3), ("F", 4)), (("E", 4), ("F", 3))]
# P end, interior -> outer. Used only to position the P-DNA helix; the helix
# keeps its own nucleotides at these two levels.
P_END_LEVELS = [(("E", 5), ("F", 2)), (("E", 6), ("F", 1))]

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
    'ATOM     42  P    DG E   3       2.541   9.337  32.250  1.00  0.00           P\n',
    'ATOM     43  OP1  DG E   3       2.877  10.540  31.457  1.00  0.00           O\n',
    'ATOM     44  OP2  DG E   3       3.456   8.928  33.340  1.00  0.00           O\n',
    "ATOM     45  O5'  DG E   3       1.054   9.514  32.826  1.00  0.00           O\n",
    "ATOM     46  C5'  DG E   3       0.294   8.373  33.239  1.00  0.00           C\n",
    "ATOM     47  C4'  DG E   3       0.963   7.678  34.409  1.00  0.00           C\n",
    "ATOM     48  O4'  DG E   3       1.220   8.646  35.461  1.00  0.00           O\n",
    "ATOM     49  C3'  DG E   3       0.152   6.595  35.092  1.00  0.00           C\n",
    "ATOM     50  O3'  DG E   3       0.244   5.357  34.381  1.00  0.00           O\n",
    "ATOM     51  C2'  DG E   3       0.861   6.533  36.437  1.00  0.00           C\n",
    "ATOM     52  C1'  DG E   3       1.052   8.019  36.723  1.00  0.00           C\n",
    'ATOM     53  N9   DG E   3      -0.097   8.614  37.398  1.00  0.00           N\n',
    'ATOM     54  C8   DG E   3      -1.048   9.441  36.850  1.00  0.00           C\n',
    'ATOM     55  N7   DG E   3      -1.968   9.808  37.701  1.00  0.00           N\n',
    'ATOM     56  C5   DG E   3      -1.606   9.178  38.886  1.00  0.00           C\n',
    'ATOM     57  C6   DG E   3      -2.225   9.200  40.159  1.00  0.00           C\n',
    'ATOM     58  O6   DG E   3      -3.251   9.802  40.504  1.00  0.00           O\n',
    'ATOM     59  N1   DG E   3      -1.531   8.423  41.083  1.00  0.00           N\n',
    'ATOM     60  C2   DG E   3      -0.385   7.714  40.811  1.00  0.00           C\n',
    'ATOM     61  N2   DG E   3       0.143   7.022  41.831  1.00  0.00           N\n',
    'ATOM     62  N3   DG E   3       0.204   7.685  39.624  1.00  0.00           N\n',
    'ATOM     63  C4   DG E   3      -0.459   8.437  38.713  1.00  0.00           C\n',
    'ATOM     64  P    DG E   4      -0.602   4.081  34.881  1.00  0.00           P\n',
    'ATOM     65  OP1  DG E   4      -2.040   4.411  34.776  1.00  0.00           O\n',
    'ATOM     66  OP2  DG E   4      -0.040   3.658  36.183  1.00  0.00           O\n',
    "ATOM     67  O5'  DG E   4      -0.272   2.938  33.809  1.00  0.00           O\n",
    "ATOM     68  C5'  DG E   4       0.763   1.988  34.064  1.00  0.00           C\n",
    "ATOM     69  C4'  DG E   4       0.265   0.837  34.931  1.00  0.00           C\n",
    "ATOM     70  O4'  DG E   4      -1.126   1.049  35.305  1.00  0.00           O\n",
    "ATOM     71  C3'  DG E   4       0.303  -0.547  34.270  1.00  0.00           C\n",
    "ATOM     72  O3'  DG E   4       0.652  -1.526  35.231  1.00  0.00           O\n",
    "ATOM     73  C2'  DG E   4      -1.143  -0.736  33.840  1.00  0.00           C\n",
    "ATOM     74  C1'  DG E   4      -1.835  -0.140  35.047  1.00  0.00           C\n",
    'ATOM     75  N9   DG E   4      -3.240   0.172  34.838  1.00  0.00           N\n',
    'ATOM     76  C8   DG E   4      -3.761   1.217  34.116  1.00  0.00           C\n',
    'ATOM     77  N7   DG E   4      -5.065   1.242  34.109  1.00  0.00           N\n',
    'ATOM     78  C5   DG E   4      -5.430   0.144  34.878  1.00  0.00           C\n',
    'ATOM     79  C6   DG E   4      -6.713  -0.341  35.226  1.00  0.00           C\n',
    'ATOM     80  O6   DG E   4      -7.820   0.120  34.915  1.00  0.00           O\n',
    'ATOM     81  N1   DG E   4      -6.634  -1.481  36.020  1.00  0.00           N\n',
    'ATOM     82  C2   DG E   4      -5.464  -2.074  36.423  1.00  0.00           C\n',
    'ATOM     83  N2   DG E   4      -5.588  -3.165  37.181  1.00  0.00           N\n',
    'ATOM     84  N3   DG E   4      -4.257  -1.632  36.105  1.00  0.00           N\n',
    'ATOM     85  C4   DG E   4      -4.317  -0.523  35.333  1.00  0.00           C\n',
    'ATOM     86  P    DG E   5       2.043  -1.441  36.025  1.00  0.00           P\n',
    'ATOM     87  OP1  DG E   5       3.187  -0.788  35.234  1.00  0.00           O\n',
    'ATOM     88  OP2  DG E   5       2.567  -2.831  36.422  1.00  0.00           O\n',
    "ATOM     89  O5'  DG E   5       1.707  -0.573  37.383  1.00  0.00           O\n",
    "ATOM     90  C5'  DG E   5       2.516   0.516  37.904  1.00  0.00           C\n",
    "ATOM     91  C4'  DG E   5       2.272   0.820  39.419  1.00  0.00           C\n",
    "ATOM     92  O4'  DG E   5       3.517   1.343  39.966  1.00  0.00           O\n",
    "ATOM     93  C3'  DG E   5       1.183   1.856  39.829  1.00  0.00           C\n",
    "ATOM     94  O3'  DG E   5       0.715   1.584  41.189  1.00  0.00           O\n",
    "ATOM     95  C2'  DG E   5       2.007   3.148  39.810  1.00  0.00           C\n",
    "ATOM     96  C1'  DG E   5       3.329   2.695  40.444  1.00  0.00           C\n",
    'ATOM     97  N9   DG E   5       4.495   3.470  40.028  1.00  0.00           N\n',
    'ATOM     98  C8   DG E   5       5.283   3.253  38.922  1.00  0.00           C\n',
    'ATOM     99  N7   DG E   5       6.258   4.114  38.812  1.00  0.00           N\n',
    'ATOM    100  C5   DG E   5       6.106   4.951  39.910  1.00  0.00           C\n',
    'ATOM    101  C6   DG E   5       6.867   6.073  40.325  1.00  0.00           C\n',
    'ATOM    102  O6   DG E   5       7.863   6.570  39.785  1.00  0.00           O\n',
    'ATOM    103  N1   DG E   5       6.364   6.630  41.496  1.00  0.00           N\n',
    'ATOM    104  C2   DG E   5       5.268   6.166  42.181  1.00  0.00           C\n',
    'ATOM    105  N2   DG E   5       4.939   6.840  43.294  1.00  0.00           N\n',
    'ATOM    106  N3   DG E   5       4.549   5.122  41.805  1.00  0.00           N\n',
    'ATOM    107  C4   DG E   5       5.020   4.565  40.670  1.00  0.00           C\n',
    'ATOM    108  P    DG E   6      -0.289   2.483  42.125  1.00  0.00           P\n',
    'ATOM    109  OP1  DG E   6      -1.539   2.900  41.334  1.00  0.00           O\n',
    'ATOM    110  OP2  DG E   6       0.396   3.801  42.522  1.00  0.00           O\n',
    "ATOM    111  O5'  DG E   6      -0.711   1.654  43.483  1.00  0.00           O\n",
    "ATOM    112  C5'  DG E   6      -2.063   1.530  44.004  1.00  0.00           C\n",
    "ATOM    113  C4'  DG E   6      -2.126   1.146  45.519  1.00  0.00           C\n",
    "ATOM    114  O4'  DG E   6      -3.346   1.724  46.066  1.00  0.00           O\n",
    "ATOM    115  C3'  DG E   6      -2.172  -0.357  45.929  1.00  0.00           C\n",
    "ATOM    116  O3'  DG E   6      -1.657  -0.524  47.289  1.00  0.00           O\n",
    "ATOM    117  C2'  DG E   6      -3.684  -0.605  45.910  1.00  0.00           C\n",
    "ATOM    118  C1'  DG E   6      -4.229   0.682  46.543  1.00  0.00           C\n",
    'ATOM    119  N9   DG E   6      -5.583   1.034  46.128  1.00  0.00           N\n',
    'ATOM    120  C8   DG E   6      -5.948   1.766  45.022  1.00  0.00           C\n',
    'ATOM    121  N7   DG E   6      -7.239   1.918  44.912  1.00  0.00           N\n',
    'ATOM    122  C5   DG E   6      -7.761   1.246  46.010  1.00  0.00           C\n',
    'ATOM    123  C6   DG E   6      -9.105   1.065  46.425  1.00  0.00           C\n',
    'ATOM    124  O6   DG E   6     -10.140   1.476  45.885  1.00  0.00           O\n',
    'ATOM    125  N1   DG E   6      -9.185   0.318  47.596  1.00  0.00           N\n',
    'ATOM    126  C2   DG E   6      -8.108  -0.189  48.281  1.00  0.00           C\n',
    'ATOM    127  N2   DG E   6      -8.390  -0.883  49.394  1.00  0.00           N\n',
    'ATOM    128  N3   DG E   6      -6.850  -0.028  47.905  1.00  0.00           N\n',
    'ATOM    129  C4   DG E   6      -6.750   0.695  46.770  1.00  0.00           C\n',
    'TER     130       DG E   6\n',
    'ATOM    131  P    DG F   1      -1.033   2.277  51.275  1.00  0.00           P\n',
    'ATOM    132  OP1  DG F   1      -0.177   3.278  52.066  1.00  0.00           O\n',
    'ATOM    133  OP2  DG F   1      -2.299   3.052  50.878  1.00  0.00           O\n',
    "ATOM    134  O5'  DG F   1      -0.243   1.784  49.917  1.00  0.00           O\n",
    "ATOM    135  C5'  DG F   1       0.978   2.374  49.396  1.00  0.00           C\n",
    "ATOM    136  C4'  DG F   1       1.231   2.078  47.881  1.00  0.00           C\n",
    "ATOM    137  O4'  DG F   1       1.978   3.203  47.334  1.00  0.00           O\n",
    "ATOM    138  C3'  DG F   1       2.045   0.814  47.471  1.00  0.00           C\n",
    "ATOM    139  O3'  DG F   1       1.690   0.406  46.111  1.00  0.00           O\n",
    "ATOM    140  C2'  DG F   1       3.468   1.381  47.490  1.00  0.00           C\n",
    "ATOM    141  C1'  DG F   1       3.271   2.765  46.857  1.00  0.00           C\n",
    'ATOM    142  N9   DG F   1       4.251   3.765  47.273  1.00  0.00           N\n',
    'ATOM    143  C8   DG F   1       4.186   4.580  48.379  1.00  0.00           C\n',
    'ATOM    144  N7   DG F   1       5.214   5.376  48.489  1.00  0.00           N\n',
    'ATOM    145  C5   DG F   1       6.007   5.070  47.391  1.00  0.00           C\n',
    'ATOM    146  C6   DG F   1       7.252   5.607  46.976  1.00  0.00           C\n',
    'ATOM    147  O6   DG F   1       7.927   6.493  47.516  1.00  0.00           O\n',
    'ATOM    148  N1   DG F   1       7.705   5.009  45.805  1.00  0.00           N\n',
    'ATOM    149  C2   DG F   1       7.044   4.019  45.120  1.00  0.00           C\n',
    'ATOM    150  N2   DG F   1       7.644   3.569  44.007  1.00  0.00           N\n',
    'ATOM    151  N3   DG F   1       5.884   3.509  45.496  1.00  0.00           N\n',
    'ATOM    152  C4   DG F   1       5.425   4.076  46.631  1.00  0.00           C\n',
    'ATOM    153  P    DG F   2       2.385  -0.749  45.175  1.00  0.00           P\n',
    'ATOM    154  OP1  DG F   2       2.560  -2.055  45.966  1.00  0.00           O\n',
    'ATOM    155  OP2  DG F   2       3.808  -0.323  44.778  1.00  0.00           O\n',
    "ATOM    156  O5'  DG F   2       1.491  -1.009  43.817  1.00  0.00           O\n",
    "ATOM    157  C5'  DG F   2       1.116  -2.313  43.296  1.00  0.00           C\n",
    "ATOM    158  C4'  DG F   2       0.727  -2.303  41.781  1.00  0.00           C\n",
    "ATOM    159  O4'  DG F   2       1.067  -3.610  41.234  1.00  0.00           O\n",
    "ATOM    160  C3'  DG F   2      -0.758  -2.066  41.371  1.00  0.00           C\n",
    "ATOM    161  O3'  DG F   2      -0.825  -1.530  40.011  1.00  0.00           O\n",
    "ATOM    162  C2'  DG F   2      -1.284  -3.505  41.390  1.00  0.00           C\n",
    "ATOM    163  C1'  DG F   2      -0.122  -4.282  40.757  1.00  0.00           C\n",
    'ATOM    164  N9   DG F   2      -0.030  -5.678  41.173  1.00  0.00           N\n',
    'ATOM    165  C8   DG F   2       0.620  -6.173  42.279  1.00  0.00           C\n',
    'ATOM    166  N7   DG F   2       0.527  -7.471  42.389  1.00  0.00           N\n',
    'ATOM    167  C5   DG F   2      -0.230  -7.857  41.291  1.00  0.00           C\n',
    'ATOM    168  C6   DG F   2      -0.660  -9.143  40.876  1.00  0.00           C\n',
    'ATOM    169  O6   DG F   2      -0.450 -10.237  41.416  1.00  0.00           O\n',
    'ATOM    170  N1   DG F   2      -1.408  -9.082  39.705  1.00  0.00           N\n',
    'ATOM    171  C2   DG F   2      -1.705  -7.929  39.020  1.00  0.00           C\n',
    'ATOM    172  N2   DG F   2      -2.440  -8.076  37.907  1.00  0.00           N\n',
    'ATOM    173  N3   DG F   2      -1.311  -6.724  39.396  1.00  0.00           N\n',
    'ATOM    174  C4   DG F   2      -0.582  -6.760  40.531  1.00  0.00           C\n',
    'ATOM    175  P    DG F   3      -2.224  -1.537  39.213  1.00  0.00           P\n',
    'ATOM    176  OP1  DG F   3      -2.328  -2.827  38.496  1.00  0.00           O\n',
    'ATOM    177  OP2  DG F   3      -2.303  -0.273  38.448  1.00  0.00           O\n',
    "ATOM    178  O5'  DG F   3      -3.334  -1.477  40.372  1.00  0.00           O\n",
    "ATOM    179  C5'  DG F   3      -4.679  -1.899  40.115  1.00  0.00           C\n",
    "ATOM    180  C4'  DG F   3      -4.841  -3.394  40.358  1.00  0.00           C\n",
    "ATOM    181  O4'  DG F   3      -4.542  -3.717  41.742  1.00  0.00           O\n",
    "ATOM    182  C3'  DG F   3      -6.246  -3.957  40.108  1.00  0.00           C\n",
    "ATOM    183  O3'  DG F   3      -6.127  -5.246  39.538  1.00  0.00           O\n",
    "ATOM    184  C2'  DG F   3      -6.818  -4.045  41.513  1.00  0.00           C\n",
    "ATOM    185  C1'  DG F   3      -5.584  -4.525  42.241  1.00  0.00           C\n",
    'ATOM    186  N9   DG F   3      -5.647  -4.381  43.685  1.00  0.00           N\n',
    'ATOM    187  C8   DG F   3      -6.670  -3.841  44.424  1.00  0.00           C\n',
    'ATOM    188  N7   DG F   3      -6.438  -3.850  45.708  1.00  0.00           N\n',
    'ATOM    189  C5   DG F   3      -5.182  -4.435  45.822  1.00  0.00           C\n',
    'ATOM    190  C6   DG F   3      -4.403  -4.712  46.970  1.00  0.00           C\n',
    'ATOM    191  O6   DG F   3      -4.679  -4.488  48.157  1.00  0.00           O\n',
    'ATOM    192  N1   DG F   3      -3.191  -5.312  46.638  1.00  0.00           N\n',
    'ATOM    193  C2   DG F   3      -2.783  -5.606  45.359  1.00  0.00           C\n',
    'ATOM    194  N2   DG F   3      -1.580  -6.185  45.236  1.00  0.00           N\n',
    'ATOM    195  N3   DG F   3      -3.503  -5.351  44.277  1.00  0.00           N\n',
    'ATOM    196  C4   DG F   3      -4.686  -4.767  44.584  1.00  0.00           C\n',
    'ATOM    197  P    DG F   4      -7.378  -5.930  38.801  1.00  0.00           P\n',
    'ATOM    198  OP1  DG F   4      -8.547  -5.040  38.979  1.00  0.00           O\n',
    'ATOM    199  OP2  DG F   4      -7.451  -7.339  39.248  1.00  0.00           O\n',
    "ATOM    200  O5'  DG F   4      -6.955  -5.918  37.259  1.00  0.00           O\n",
    "ATOM    201  C5'  DG F   4      -5.648  -6.343  36.866  1.00  0.00           C\n",
    "ATOM    202  C4'  DG F   4      -5.374  -5.929  35.432  1.00  0.00           C\n",
    "ATOM    203  O4'  DG F   4      -4.125  -6.514  34.967  1.00  0.00           O\n",
    "ATOM    204  C3'  DG F   4      -6.431  -6.378  34.427  1.00  0.00           C\n",
    "ATOM    205  O3'  DG F   4      -6.639  -5.355  33.443  1.00  0.00           O\n",
    "ATOM    206  C2'  DG F   4      -5.791  -7.620  33.811  1.00  0.00           C\n",
    "ATOM    207  C1'  DG F   4      -4.361  -7.135  33.724  1.00  0.00           C\n",
    'ATOM    208  N9   DG F   4      -3.371  -8.185  33.521  1.00  0.00           N\n',
    'ATOM    209  C8   DG F   4      -3.581  -9.543  33.502  1.00  0.00           C\n',
    'ATOM    210  N7   DG F   4      -2.492 -10.232  33.294  1.00  0.00           N\n',
    'ATOM    211  C5   DG F   4      -1.500  -9.266  33.164  1.00  0.00           C\n',
    'ATOM    212  C6   DG F   4      -0.111  -9.407  32.926  1.00  0.00           C\n',
    'ATOM    213  O6   DG F   4       0.541 -10.449  32.778  1.00  0.00           O\n',
    'ATOM    214  N1   DG F   4       0.529  -8.171  32.865  1.00  0.00           N\n',
    'ATOM    215  C2   DG F   4      -0.093  -6.955  33.015  1.00  0.00           C\n',
    'ATOM    216  N2   DG F   4       0.690  -5.870  32.925  1.00  0.00           N\n',
    'ATOM    217  N3   DG F   4      -1.391  -6.810  33.238  1.00  0.00           N\n',
    'ATOM    218  C4   DG F   4      -2.028  -8.004  33.300  1.00  0.00           C\n',
    'ATOM    219  P    DG F   5      -7.576  -5.685  32.185  1.00  0.00           P\n',
    'ATOM    220  OP1  DG F   5      -8.748  -6.381  32.763  1.00  0.00           O\n',
    'ATOM    221  OP2  DG F   5      -6.730  -6.370  31.181  1.00  0.00           O\n',
    "ATOM    222  O5'  DG F   5      -8.048  -4.315  31.570  1.00  0.00           O\n",
    "ATOM    223  C5'  DG F   5      -8.097  -3.236  32.431  1.00  0.00           C\n",
    "ATOM    224  C4'  DG F   5      -7.596  -1.973  31.779  1.00  0.00           C\n",
    "ATOM    225  O4'  DG F   5      -6.153  -1.801  31.879  1.00  0.00           O\n",
    "ATOM    226  C3'  DG F   5      -7.888  -1.779  30.290  1.00  0.00           C\n",
    "ATOM    227  O3'  DG F   5      -8.409  -0.543  29.993  1.00  0.00           O\n",
    "ATOM    228  C2'  DG F   5      -6.563  -2.088  29.593  1.00  0.00           C\n",
    "ATOM    229  C1'  DG F   5      -5.662  -1.397  30.610  1.00  0.00           C\n",
    'ATOM    230  N9   DG F   5      -4.261  -1.796  30.518  1.00  0.00           N\n',
    'ATOM    231  C8   DG F   5      -3.770  -3.073  30.375  1.00  0.00           C\n',
    'ATOM    232  N7   DG F   5      -2.467  -3.117  30.322  1.00  0.00           N\n',
    'ATOM    233  C5   DG F   5      -2.071  -1.790  30.438  1.00  0.00           C\n',
    'ATOM    234  C6   DG F   5      -0.776  -1.214  30.446  1.00  0.00           C\n',
    'ATOM    235  O6   DG F   5       0.320  -1.781  30.350  1.00  0.00           O\n',
    'ATOM    236  N1   DG F   5      -0.831   0.169  30.584  1.00  0.00           N\n',
    'ATOM    237  C2   DG F   5      -1.985   0.904  30.699  1.00  0.00           C\n',
    'ATOM    238  N2   DG F   5      -1.831   2.231  30.824  1.00  0.00           N\n',
    'ATOM    239  N3   DG F   5      -3.199   0.379  30.692  1.00  0.00           N\n',
    'ATOM    240  C4   DG F   5      -3.169  -0.963  30.559  1.00  0.00           C\n',
    'ATOM    241  P    DC F   6      -9.462  -0.430  28.785  1.00  0.00           P\n',
    'ATOM    242  OP1  DC F   6     -10.823  -0.345  29.362  1.00  0.00           O\n',
    'ATOM    243  OP2  DC F   6      -9.149  -1.473  27.781  1.00  0.00           O\n',
    "ATOM    244  O5'  DC F   6      -9.081   0.968  28.170  1.00  0.00           O\n",
    "ATOM    245  C5'  DC F   6      -8.514   1.887  29.032  1.00  0.00           C\n",
    "ATOM    246  C4'  DC F   6      -7.394   2.638  28.378  1.00  0.00           C\n",
    "ATOM    247  O4'  DC F   6      -6.101   1.974  28.479  1.00  0.00           O\n",
    "ATOM    248  C3'  DC F   6      -7.527   2.962  26.889  1.00  0.00           C\n",
    "ATOM    249  O3'  DC F   6      -7.255   4.288  26.592  1.00  0.00           O\n",
    "ATOM    250  C2'  DC F   6      -6.602   1.965  26.192  1.00  0.00           C\n",
    "ATOM    251  C1'  DC F   6      -5.468   2.035  27.209  1.00  0.00           C\n",
    'ATOM    252  N1   DC F   6      -4.523   0.912  27.118  1.00  0.00           N\n',
    'ATOM    253  C2   DC F   6      -3.157   1.195  27.180  1.00  0.00           C\n',
    'ATOM    254  O2   DC F   6      -2.795   2.374  27.308  1.00  0.00           O\n',
    'ATOM    255  N3   DC F   6      -2.266   0.180  27.099  1.00  0.00           N\n',
    'ATOM    256  C4   DC F   6      -2.699  -1.075  26.960  1.00  0.00           C\n',
    'ATOM    257  N4   DC F   6      -1.784  -2.046  26.883  1.00  0.00           N\n',
    'ATOM    258  C5   DC F   6      -4.085  -1.392  26.894  1.00  0.00           C\n',
    'ATOM    259  C6   DC F   6      -4.955  -0.377  26.975  1.00  0.00           C\n',
    'TER     260       DC F   6\n',
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
