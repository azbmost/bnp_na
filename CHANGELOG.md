# Change Log

This file records the public GitHub-ready `bnp_na` version history from the repository preparation work onward.

## V13.16

- Added a `Regularize phosphates` pipeline checkbox directly below `Run phenix.geometry_minimization`; its per-helix-type default follows the minimization default.
- Added C1'-derived, per-chain helical regularization of every movable atom selected by `min_P_C5.params`: `P`, `OP1`, `OP2`, `O5'`, `C5'`, and `O3'`.
- Affected terminal atoms are excluded from their consensus-location calculation and are regularized afterward from internal positions, preserving periodic phosphate-linkage bond geometry.
- Added a standalone `Regularize phosphates` dialog under `Other tools` and a matching `bnp_na_lib/regularize_phosphates.py` CLI.
- Added regression coverage for internal periodicity, 5'-terminal exclusion, phosphate-only 3'-terminal propagation, and short-helix validation.

## V13.15

- Added `--axis_range` and `--axis-range` to `Measure angle around axis` for restricting a PDB-fitted axis to one or more residue ranges.
- Added an `Axis residue ranges` field to the angle-tool GUI.
- All listed ranges contribute selected axis atoms to the PCA fit; the written start-to-end order of the first range sets the positive axis direction.
- Added parsing, atom-selection, direction-orientation, and CLI-alias regression coverage.

## V13.14

- Fixed a native macOS Tk 8.6 crash when using the file dialogs in `Measure angle around axis`.
- The angle tool now omits native file-type restrictions on macOS, bypassing Tk's unsafe UTType conversion while retaining file filters on other platforms.
- The output dialog still supplies `.bild` as its default extension on macOS.
- Added regression coverage for macOS input/output dialog options and the cross-platform filter structure.

## V13.13

- `combine_PDB` now allows the same input PDB path to appear two or more times.
- Every repeated occurrence is processed independently and receives fresh consecutive chain IDs, atom/TER serials, connectivity references, and updated `LINK`/`REMARK`/`HET` metadata.
- Added regression coverage showing that a repeated two-chain PDB produces `A/B` for its first occurrence and `C/D` for its second occurrence.

## V13.12

- The `combine_PDB` tool now preserves `LINK` records and updates both fixed-column endpoint chain IDs.
- Source `REMARK` records are retained, with automatic chain updates for `re_helix` `CHAIN_RANGE`, `CHAIN_RESIDUES`, `JUNCTION`, and `SPECIAL` current-output fields, bnp_na `CHAIN` annotations, and common colon/DSSR-style residue labels.
- `re_helix` provenance fields such as `COMMAND`, `source=`, and `original_*=` remain unchanged so their original identities are not corrupted.
- Related `HET` and `HETNAM` records are retained, with `HET` chain IDs updated for linker residues.
- Added regression and real-example validation using `Inverted_TT-A38-Olson_rex.pdb`, including all 41 of its `LINK` records.

## V13.11

- Added a `combine_PDB` tool in `Other tools` for combining two or more PDB coordinate files.
- Added a dynamic input-count dropdown with scrollable file fields for 2 through 26 input PDB files.
- Source chains are reassigned consecutive uppercase IDs (`A` through `Z`) in input-file and first-appearance order.
- Combined outputs globally renumber `ATOM`, `HETATM`, and `TER` records and remap companion-record serials and `CONECT` references.
- Added `bnp_na_lib/combine_pdb.py` as a reusable helper and direct command-line tool.

## V13.10

- Added an optional `O3' before 5' phosphate` mode to the Add phosphates tool and its CLI.
- The tool reports the preceding `O3'` status separately for each chain and writes the atom as a one-atom residue `n-1` before the 5' phosphate on residue `n`.
- A new `O3'` can be added alongside a newly generated 5' phosphate or onto an existing 5' phosphate using neighboring backbone geometry.
- Add-phosphates provenance now reports phosphate-group and preceding-`O3'` additions separately while retaining file-order atom renumbering.

## V13.9

- Added an `Add phosphates` tool in `Other tools` for reporting and adding missing 5' and 3' terminal phosphates by chain.
- Added `bnp_na_lib/add_phosphates.py`, with a direct CLI and reusable helper API for terminal phosphate status reports and neighbor-geometry phosphate placement.
- The 5' addition fits the second residue sugar onto the first residue sugar, then transforms the second residue's `P/OP1/OP2` onto the first residue.
- The 3' addition fits the penultimate residue sugar onto the terminal residue sugar, then transforms the terminal residue's `P/OP1/OP2` into a phosphate-only `N+1` residue using the terminal residue name.
- Phosphate-added PDB outputs are renumbered in file order, with existing `CONECT` records remapped to the new atom serials.
- Corrected the A-RNA DSSR rebuild command to use `--backbone=RNA` instead of `--backbone=A-RNA`, preserving RNA `O2'` atoms in the rebuilt PDB.

## V13.8

- Renamed the main GUI's bottom `Analysis tools` area to `Other tools`.
- Added a `Get helical-axis info` tool in `Other tools`.
- Added `bnp_na_lib/helical_axis_info.py`, which filters a PDB to two selected chain IDs, runs `x3dna-dssr --more`, reports DSSR axis start/end points, axis vector, unit vector, and angle to a user-provided reference vector.
- Added optional Chimera/ChimeraX `.bild` output for the selected-chain DSSR helical axis.
- Made the helical-axis info BILD filename auto-update whenever the input PDB or chain IDs change.
- Added start-to-end distance reporting plus an optional helix-length-in-bp field for estimating the full helix length as `distance / (bp - 1) * bp`.
- Made chain-ID order define selected-chain helical-axis direction, so `A B` and `B A` report opposite vectors.
- Added controls for drawing the reference vector in the BILD output and for setting its drawing length.
- Added an `Align helix to z` tool in `Other tools` for applying the existing DSSR `align2z.py` workflow to any input helix PDB.
- Renamed `Open helical-axis angle tool` to `Measure angle around axis` and arranged all `Other tools` launchers on one row.
- Replaced visible `Other tools` descriptions with light-blue `?` help buttons.
- Changed `Other tools` log output to append instead of replacing previous run records.
- Updated README usage and repository-layout documentation for the new selected-chain helical-axis info tool.
- Renamed the bundled helical-axis angle tool to `bnp_na_lib/angle_helical_axisV2_2.py`.
- Added optional region-defined 2-fold symmetry-axis calculation to the helical-axis angle tool, including symmetry-axis points, +90-degree rotated points, RMSD reporting, GUI controls, CLI flags, and BILD output.
- Made 2-fold symmetry output opt-in with GUI radio buttons, disabled non-relevant symmetry fields by default, allowed blank symmetry regions to use the whole two-chain model, and added BILD arrow color/direction notes to the output report.
- Changed the helical-axis angle tool window from a transient child dialog to a normal top-level tool window to avoid disappearing when moved between displays.
- Set the default 2-fold symmetry point radius to `15.0 Å` in both the GUI and CLI.

## V13.7

- Added triplex conversion from an input duplex PDB using bundled `convert_to_triplex_pdbV2_1.py`.
- Added `bnp_na_lib/build_triplex.py` so the main app can preview duplex chains/sequences, run the conversion, capture its log, and report the output PDB path.
- Added a `Triplex converter` launcher directly to the right of the `B-Z builder` button in the helix-type row.
- Added a native triplex converter dialog with input/output PDB selection, strand-I chain, residue range, antiparallel/parallel mode, optional strand-II/strand-III chain IDs, and strand-III residue-number start.
- Added triplex final-PDB `REMARK BNP_NA...` provenance records.
- Moved tool-button explanation text so each note immediately follows the button it describes.
- Mirrored helical-axis angle and XYZ axes BILD tool CLI commands/results into the main `bnp_na` log box.
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
