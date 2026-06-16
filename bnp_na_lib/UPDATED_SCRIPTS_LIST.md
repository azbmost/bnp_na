# bnp_na V13.2 updated folder

Top-level app folder contains only:

- `bnp_na.py` — main GUI/controller and `-v` / `--version` entry point.
- `bnp_na_lib/` — helper modules and support files.

## Helper folder: `bnp_na_lib/`

- `build_common.py` — shared sequence expansion, DSSR rebuild/fiber wrappers, DSSR installation check, helical table writer, phenix.geometry_minimization wrapper, default 12-parameter rows.
- `build_bdna.py` — B-DNA builder using `x3dna-dssr rebuild --backbone=B-DNA --par-type=heli`; optional phenix.geometry_minimization; DSSR-based align-to-Z.
- `build_adna.py` — A-DNA builder using `x3dna-dssr rebuild --backbone=A-DNA --par-type=heli`; optional phenix.geometry_minimization; no twist warp.
- `build_arna.py` — A-RNA builder using `x3dna-dssr rebuild --backbone=A-RNA --par-type=heli`; optional phenix.geometry_minimization.
- `build_zdna.py` — Z-DNA builder using DSSR fiber; no phenix.geometry_minimization option.
- `align2z.py` — DSSR-dependent align-to-Z module using `x3dna-dssr --more` point-one/point-two endpoints.
- `geometry_utils.py` — shared rotation/vector helpers.
- `na_placer.py` — final orient/place transformation after +Z alignment.
- `pdb_inv_rotV2.py` — optional inversion/reflection helper for mirror-image L-form models after align-to-Z and before final placement.
- `pdb_name_standard.py` — nucleotide residue/atom-name normalization; changed from the previous script name `pdb_make_dna_v3_2.py`.
- `edit_pdb_atom.py` — PDB parser/writer helper.
- `min_P_C5.params` — default params file shown in the GUI phenix.geometry_minimization field.
- `__init__.py` — helper package marker.

## V13.2 GUI changes

1. Log output is embedded in the main GUI.
2. Startup immediately checks `x3dna-dssr` and prints the executable/version/help output into the embedded log.
3. `Customize DSSR parameters` pre-fills defaults when no values were previously saved, and keeps saved values for the selected NA type.
4. Output folder is user-selectable. Final placed PDB is written to the selected folder; intermediate files are written to `<selected folder>/tmp_file/`.
5. B-DNA, A-DNA, and A-RNA have a `phenix.geometry_minimization` checkbox and params-file field. Z-DNA does not provide this option.
6. GUI text uses `GIDEON` in all capitals.
7. Version is `bnp_na V13.2`; run `python bnp_na.py -v` or `python bnp_na.py --version` to print it.
8. Sequence length is updated in the GUI after sequence input changes.
9. Optional mirror-image L-form generation is available after align-to-Z and before final orient/place, using `pdb_inv_rotV2.py` i-mode and o-mode operations.
