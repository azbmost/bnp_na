# Change Log

This file records the public GitHub-ready `bnp_na` version history from the repository preparation work onward.

## Unreleased

- Refined the main GUI typography: the sequence input label is no longer bold, while functional module titles use bold label-frame headings.
- Updated the mirror-image chirality module title to `Mirror-image L-form chirality (L-DNA)`.
- Replaced the one-line current helical-parameter display with a compact three-row horizontal table showing current values, defaults, and default/custom source status.
- Added default-value labels beside every helical-parameter entry in the customization dialog.
- Reduced vertical padding between input rows and switched explanatory/status text to smaller GUI fonts to free more vertical space.

## V13.5

- Renamed the bundled helical-axis angle tool from `bnp_na_lib/angle_helical_axisV2.py` to `bnp_na_lib/angle_helical_axisV2_1.py` so the filename indicates the incorporated script version.
- Incorporated the `angle_helical_axisV2_1.py` updates into the bundled helical-axis angle tool.
- Added `--axis-margin` to the helical-axis angle tool CLI.
- Added an `Axis drawing margin` GUI field to control how far the BILD axis arrow extends beyond the fitted/custom axis span.
- Added explanatory `.comment` lines before each drawn helical-axis BILD object.
- Updated README documentation for the new axis-margin behavior.
- Added this top-level change log.

## V13.4

- Added `bnp_na_lib/xyz_bild.py` for writing simple XYZ coordinate-axis BILD files.
- Added a `Write XYZ axes BILD` launcher to the main GUI.
- Moved the helical-axis angle tool and XYZ axes BILD tool launchers to the bottom `Analysis tools` section immediately above `Log output`.
- Documented the XYZ axes BILD GUI and CLI usage.

## V13.3

- Added `bnp_na_lib/angle_helical_axisV2.py` for measuring around-axis angles between two atom or XYZ points.
- Added GUI and CLI support for fitted PDB axes and custom axis point/vector definitions.
- Added Chimera/ChimeraX BILD output for the helical axis, radial vectors, points, and point projections.
- Added a main-GUI launcher for the helical-axis angle tool.

## V13.2

- Added optional mirror-image L-DNA/L-RNA modeling with `bnp_na_lib/pdb_inv_rotV2.py`.
- Added GUI controls for both inversion (`i`) mode and coordinate-plane reflection (`o`) mode operations.
- Applied the chirality operation after align-to-Z and before final placement/orientation.
- Added final-PDB REMARK records identifying `bnp_na`, DiLiuLab/AZBMOST provenance, L-form status, operation, and L-DNA/L-RNA residue ranges.
- Expanded README documentation for L-form modeling, default helical parameters, and the default Phenix minimization params file.

## V13.1

- Prepared the project as a public GitHub repository under the AZBMOST organization.
- Added MIT licensing and GitHub-oriented README setup instructions.
- Added clone, pull, dependency-install, GUI-run, and executable/PyInstaller guidance.
- Added app icon assets with a right-handed helix cue while keeping the script runnable without icon assets.
- Renamed the PDB naming helper to `pdb_name_standard.py` and documented that it was changed from the previous `pdb_make_dna_v3_2.py` script.
- Kept generated output, caches, and the `test/` folder out of version control.
