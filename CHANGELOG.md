# Change Log

This file records the public GitHub-ready `bnp_na` version history from the repository preparation work onward.

## V13.7

- Added triplex conversion from an input duplex PDB using bundled `convert_to_triplex_pdbV2_1.py`.
- Added `bnp_na_lib/build_triplex.py` so the main app can preview duplex chains/sequences, run the conversion, capture its log, and report the output PDB path.
- Added a `Triplex converter` launcher directly to the right of the `B-Z builder` button in the helix-type row.
- Added a native triplex converter dialog with input/output PDB selection, strand-I chain, residue range, antiparallel/parallel mode, optional strand-II/strand-III chain IDs, and strand-III residue-number start.
- Added triplex final-PDB `REMARK BNP_NA...` provenance records.
- Moved tool-button explanation text so each note immediately follows the button it describes.
- Updated README guidance for triplex conventions, supported `G·G-C` and `T·A-T` motifs, GUI usage, CLI usage, and troubleshooting.

## V13.6

- Added B-Z structure building from alternating B-DNA/Z-DNA PDB inputs using bundled `make_BZV2_3.py` and `core_BZ.py`.
- Added `bnp_na_lib/build_bz.py` so the main app can run the B-Z pipeline, capture its log, and report final/raw output paths.
- Added a `B-Z builder` launcher directly after the four helix-type choices, with a minimal same-row hint and a dedicated dialog for input PDBs, output path, axis correction mode, axis source, and Z-DNA terminal auto-trim.
- Added B-Z final-PDB `REMARK BNP_NA...` provenance records.
- Put the two analysis-tool launcher buttons on one row to save vertical GUI space.
- Compacted the current DSSR helical-parameter display by placing the table beside the Customize button and moving the status text below the button.
- Expanded README guidance for B-Z input order, axis correction, true collinear bond-length tradeoffs, axis source, Z-DNA terminal auto-trim, direct CLI usage, and troubleshooting.
- Refined the main GUI typography: the sequence input label is no longer bold, while functional module titles use bold label-frame headings.
- Updated the mirror-image chirality module title to `Mirror-image L-form chirality (L-DNA)`.
- Replaced the one-line current helical-parameter display with a compact three-row horizontal table showing current values, defaults, and default/custom source status.
- Added default-value labels beside every helical-parameter entry in the customization dialog.
- Reduced vertical padding between input rows and switched explanatory/status text to smaller GUI fonts to free more vertical space.
- Updated the first-line GUI title to identify `bnp_na` as Module #1 of the AZBMOST package.
- Added Tk app-name identity hints so task/menu labels can show `bnp_na` where the platform honors Tk application names.
- Added Å/° units to GUI helical-parameter, placement, XYZ axes BILD controls, and placement logs.
- Allowed simple arithmetic expressions such as `360/10.5` in numeric GUI fields for helical parameters, placement/orientation, Z-DNA length, and XYZ axes BILD settings.

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
