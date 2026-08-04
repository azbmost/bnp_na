# bnp_na V13.13 updated folder

Top-level app folder contains only:

- `bnp_na.py` — main GUI/controller and `-v` / `--version` entry point.
- `CHANGELOG.md` — public version-by-version change log.
- `bnp_na_lib/` — helper modules and support files.

## Helper folder: `bnp_na_lib/`

- `build_common.py` — shared sequence expansion, DSSR rebuild/fiber wrappers, DSSR installation check, helical table writer, phenix.geometry_minimization wrapper, default 12-parameter rows.
- `build_bdna.py` — B-DNA builder using `x3dna-dssr rebuild --backbone=B-DNA --par-type=heli`; optional phenix.geometry_minimization; DSSR-based align-to-Z.
- `build_adna.py` — A-DNA builder using `x3dna-dssr rebuild --backbone=A-DNA --par-type=heli`; optional phenix.geometry_minimization; no twist warp.
- `build_arna.py` — A-RNA builder using `x3dna-dssr rebuild --backbone=RNA --par-type=heli`; optional phenix.geometry_minimization.
- `build_zdna.py` — Z-DNA builder using DSSR fiber; no phenix.geometry_minimization option.
- `build_bz.py` — bnp_na wrapper for running the B-Z structure builder from the main GUI and capturing the log/output paths.
- `build_triplex.py` — bnp_na wrapper for previewing duplex chains/sequences and running the triplex converter from the main GUI.
- `add_phosphates.py` — terminal phosphate reporter and neighbor-geometry phosphate-placement helper for the main GUI's Add phosphates tool.
- `combine_pdb.py` — combines multiple PDB coordinate files with consecutive A-Z chain IDs, global serial renumbering, remapped connectivity, and updated LINK/REMARK metadata.
- `align2z.py` — DSSR-dependent align-to-Z module using `x3dna-dssr --more` point-one/point-two endpoints.
- `angle_helical_axisV2_2.py` — helical-axis radial-angle and 2-fold symmetry-axis calculator with PDB-fit/custom-axis modes and Chimera/ChimeraX BILD output.
- `helical_axis_info.py` — DSSR selected-chain helical-axis reporter with unit-vector/angle reporting and optional Chimera/ChimeraX BILD output.
- `make_BZV2_3.py` — standalone B-Z structure builder incorporated into V13.6; combines alternating B-DNA/Z-DNA PDB inputs using B-Z junction cores.
- `core_BZ.py` — bundled B-Z junction core structure data used by `make_BZV2_3.py`.
- `convert_to_triplex_pdbV2_1.py` — standalone duplex-to-triplex converter incorporated into V13.7; supports antiparallel G·G-C and parallel T·A-T triplex motifs.
- `geometry_utils.py` — shared rotation/vector helpers.
- `na_placer.py` — final orient/place transformation after +Z alignment.
- `pdb_inv_rotV2.py` — optional inversion/reflection helper for mirror-image L-form models after align-to-Z and before final placement.
- `pdb_name_standard.py` — nucleotide residue/atom-name normalization; changed from the previous script name `pdb_make_dna_v3_2.py`.
- `xyz_bild.py` — coordinate-axis BILD writer with configurable arrow length and width.
- `edit_pdb_atom.py` — PDB parser/writer helper.
- `min_P_C5.params` — default params file shown in the GUI phenix.geometry_minimization field.
- `__init__.py` — helper package marker.

## V13.13 changes

1. Version is `bnp_na V13.13`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. `combine_PDB` allows the same input path to be selected repeatedly.
3. Every occurrence receives independent chain, serial, connectivity, LINK, REMARK, and HET remapping.
4. A repeated two-chain PDB is verified to produce chains A/B and then C/D.

## V13.12 changes

1. Version is `bnp_na V13.12`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. `combine_PDB` now retains and updates both chain endpoints in source `LINK` records.
3. It retains source `REMARK` records and updates current chain references in re_helix, bnp_na, and common DSSR-style formats.
4. Original/source provenance inside re_helix remarks remains unchanged.
5. Related `HET`/`HETNAM` records are retained, with linker-residue `HET` chains updated.

## V13.11 changes

1. Version is `bnp_na V13.11`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. A `combine_PDB` button is available from `Other tools`.
3. Its input-count dropdown dynamically displays 2 through 26 scrollable PDB file fields.
4. `combine_pdb.py` assigns chains consecutively as A-Z in input-file and first-appearance order.
5. Combined coordinate records are globally renumbered, and companion records plus `CONECT` references are updated.

## V13.10 changes

1. Version is `bnp_na V13.10`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. The `Add phosphates` dialog now has an optional `O3' before 5' phosphate` checkbox.
3. `add_phosphates.py` reports this preceding atom separately and can add it to an existing or newly generated 5' phosphate.
4. For a 5' phosphate on residue `n`, the new `O3'` is written in a one-atom residue `n-1` using the first nucleotide's residue name.
5. The CLI exposes the same behavior through `--add-5prime-o3`; use `--ends none` when adding only this atom.

## V13.9 changes

1. Version is `bnp_na V13.9`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. An `Add phosphates` button is available from `Other tools`.
3. `add_phosphates.py` reports whether each chain has a first-residue 5' phosphate and a terminal phosphate-only 3' residue.
4. The tool can add selected missing 5' and/or 3' phosphates to selected chains using neighboring residue sugar-atom fits.
5. Add-phosphates outputs are renumbered in file order, and existing `CONECT` records are remapped to the new serials.
6. A-RNA DSSR rebuild now uses `--backbone=RNA`, preserving RNA sugar `O2'` atoms.

## V13.8 changes

1. Version is `bnp_na V13.8`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
2. The bottom main-GUI tool area is now titled `Other tools`.
3. A `Get helical-axis info` button is available from `Other tools`.
4. `helical_axis_info.py` filters a PDB to two requested chain IDs, runs `x3dna-dssr --more`, and reports the DSSR axis start/end points, axis vector, unit vector, and angle to a reference vector.
5. The selected-chain DSSR `.out` filename ends with `.out`, and optional `.bild` output can draw the selected-chain helical axis.
6. The helical-axis info BILD filename auto-updates whenever the input PDB or chain IDs change.
7. The report includes start-to-end distance and can estimate full helix length from an optional bp count using `distance / (bp - 1) * bp`.
8. Chain-ID order controls the selected-chain helical-axis direction, so `A B` and `B A` report opposite vectors.
9. The helical-axis info BILD output can omit the reference vector or draw it with a user-provided length.
10. An `Align helix to z` button is available from `Other tools` and uses the existing DSSR `align2z.py` workflow on any input helix PDB.
11. The helical-axis angle launcher is now named `Measure angle around axis`, and all `Other tools` launchers are arranged on one row.
12. The `Other tools` descriptions are available through light-blue `?` help buttons, and tool logs append instead of replacing previous records.

## V13.7 changes

1. Log output is embedded in the main GUI.
2. Startup immediately checks `x3dna-dssr` and prints the executable/version/help output into the embedded log.
3. `Customize DSSR parameters` pre-fills defaults when no values were previously saved, and keeps saved values for the selected NA type.
4. Output folder is user-selectable. Final placed PDB is written to the selected folder; intermediate files are written to `<selected folder>/tmp_file/`.
5. B-DNA, A-DNA, and A-RNA have a `phenix.geometry_minimization` checkbox and params-file field. Z-DNA does not provide this option.
6. GUI text uses `GIDEON` in all capitals.
7. Version is `bnp_na V13.7`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
8. Sequence length is updated in the GUI after sequence input changes.
9. Optional mirror-image L-form generation is available after align-to-Z and before final orient/place, using `pdb_inv_rotV2.py` i-mode and o-mode operations.
10. A helical-axis angle tool is available from the main GUI and can also be run directly as `bnp_na_lib/angle_helical_axisV2_2.py`.
11. The helical-axis angle tool filename now indicates the incorporated V2.2 script version; earlier public versions used `angle_helical_axisV2.py` and `angle_helical_axisV2_1.py`.
12. The helical-axis angle tool includes an axis drawing margin control, opt-in region/whole-model 2-fold symmetry-axis output, and explanatory BILD `.comment` records.
13. An XYZ axes BILD writer is available from the main GUI and can also be run directly as `bnp_na_lib/xyz_bild.py`.
14. A B-Z structure builder is available from the main GUI and can also be run directly as `bnp_na_lib/make_BZV2_3.py`.
15. The B-Z builder uses bundled `core_BZ.py` junction data and writes final/raw PDB outputs with captured logs in the main app.
16. A triplex converter is available from the main GUI and can also be run directly as `bnp_na_lib/convert_to_triplex_pdbV2_1.py`.
17. The triplex converter uses bundled base-triple templates to add strand III to an input duplex PDB and writes final PDB outputs with captured logs in the main app.
