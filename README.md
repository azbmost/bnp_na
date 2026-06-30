# bnp_na

`bnp_na` is a Tkinter GUI for building and placing nucleic acid helices. It can generate B-DNA, A-DNA, A-RNA, and Z-DNA models, normalize PDB atom/residue names, align the helix to the +Z axis, and write a final oriented/placed PDB file. It also includes a B-Z structure builder for joining existing B-DNA and Z-DNA PDB segments through B-Z junction cores, plus a triplex converter for adding a third strand to an input duplex PDB.

The current app version is `V13.8`.

## What It Does

- Builds B-DNA, A-DNA, and A-RNA from a 5' to 3' sequence using DSSR helical-parameter tables.
- Builds Z-DNA from a positive even base-pair length using DSSR fiber generation.
- Lets you customize the 12 DSSR base-pair/helical parameters for B-DNA, A-DNA, and A-RNA.
- Optionally runs `phenix.geometry_minimization` for B-DNA, A-DNA, and A-RNA.
- Normalizes nucleotide residue and atom names in generated PDB files.
- Aligns the generated helix to +Z using DSSR axis information.
- Aligns any input helix PDB to +Z from `Other tools`.
- Optionally applies an inversion/reflection operation to make mirror-image L-form nucleic-acid models.
- Builds B-Z DNA constructs from alternating B-DNA/Z-DNA PDB files using bundled B-Z junction core data.
- Converts an input duplex DNA PDB into a triplex by adding strand III over a selected residue range.
- Measures around-axis angles between two atom or XYZ points and writes Chimera/ChimeraX BILD drawings.
- Reports DSSR helical-axis start/end/unit-vector information for two selected PDB chains and optionally writes an axis BILD drawing.
- Writes simple XYZ coordinate-axis BILD helpers with configurable arrow length and width.
- Applies roll, phi, theta, x, y, z, and delta_z placement values.
- Writes final PDB files to the selected output folder and intermediate files to `<output folder>/tmp_file/`.

## Requirements

- Python 3.9 or newer.
- Tkinter, normally included with Python from python.org and many system Python installs.
- NumPy, installed with `pip install -r requirements.txt`.
- `x3dna-dssr` available on `PATH`, or installed at `/usr/local/bin/x3dna-dssr`.
- Optional: `phenix.geometry_minimization` on `PATH`, or `PHENIX_ENV` pointing to a valid Phenix environment script.

`x3dna-dssr` is required for real model generation. The GUI can open without it, but build jobs will fail until DSSR is installed and discoverable.

## Clone The Repository

`git clone` downloads a complete local copy of the repository from GitHub:

```bash
git clone https://github.com/AZBMOST/bnp_na.git
cd bnp_na
```

Install the Python runtime dependency:

```bash
python3 -m pip install -r requirements.txt
```

Run the GUI:

```bash
python3 bnp_na.py
```

Print the version without opening the GUI:

```bash
python3 bnp_na.py --version
```

See `CHANGELOG.md` for the version-by-version change log.

## Pull Updates

After cloning once, use `git pull` inside the repository folder to fetch and merge the latest changes from GitHub into your local copy:

```bash
cd bnp_na
git pull
```

If you have local edits, commit or stash them before pulling so Git can merge cleanly. `git pull` does not reinstall Python packages, so rerun `python3 -m pip install -r requirements.txt` only when `requirements.txt` changes.

## Quick Start

1. Start the GUI with `python3 bnp_na.py`.
2. Check the startup status line for `x3dna-dssr`.
3. Choose `B-DNA`, `A-DNA`, `A-RNA`, or `Z-DNA`.
4. Enter a sequence for B-DNA, A-DNA, or A-RNA, or enter an even base-pair length for Z-DNA.
5. Choose an output folder.
6. Adjust DSSR parameters, minimization, hydrogens, mirror-image L-form conversion, and placement values if needed.
7. Click `Generate`.
8. Read the embedded log for the exact commands, intermediate files, and final placed PDB path.

For B-Z junction constructs, click the `B-Z builder` button on the same row as the `B-DNA`, `A-DNA`, `A-RNA`, and `Z-DNA` choices. That tool combines existing B-DNA and Z-DNA PDB files rather than generating a single helix from the sequence field.

For triplex construction from an existing duplex PDB, click the `Triplex converter` button directly to the right of `B-Z builder`.

## GUI Field Guide

Numeric GUI fields for DSSR helical parameters, Z-DNA helix length, placement/orientation, and the XYZ axes BILD helper accept either plain numbers or simple arithmetic expressions such as `360/10.5`, `6*4`, or `(20+10)/2`. Blank DSSR helical-parameter fields use the displayed default value.

### x3dna-dssr Status

At startup, the GUI checks whether `x3dna-dssr` can be found. If it reports `FOUND`, the log also records the executable path and version/help output. If it reports `NOT FOUND`, install DSSR or add it to `PATH`.

The app uses DSSR for four jobs:

- `rebuild` for B-DNA, A-DNA, and A-RNA.
- `fiber --model=Z-DNA` for Z-DNA.
- `--more` axis extraction before align-to-Z placement.
- `--more` axis extraction for selected-chain helical-axis information in `Other tools`.
- `--more` axis extraction for the standalone `Align helix to z` tool in `Other tools`.

### Sequence

For B-DNA and A-DNA, enter DNA letters using `A`, `T`, `C`, and `G`.

For A-RNA, enter RNA letters using `A`, `U`, `C`, and `G`.

Whitespace is ignored. Compact count syntax is supported:

```text
A10T5C2G
```

This expands to ten `A` bases, five `T` bases, two `C` bases, and one `G` base. Counts must be positive integers.

Two-strand syntax such as `lead:follow` is not supported in this version. The GUI expects one 5' to 3' sequence and chooses the complementary pair labels internally.

### Z-DNA Helix Length

When `Z-DNA` is selected, the sequence field is disabled. Z-DNA is generated through DSSR's built-in Z-DNA fiber model, so the GUI asks only for helix length.

The Z-DNA length must be a positive even integer. The app passes `length / 2` as the DSSR repeat count.

### Helix Name

The helix name is optional. If you leave it blank, the app creates a default name such as `B-DNA25` or `Z-DNA20`.

Names are sanitized before file creation. Spaces and unsafe filename characters are converted to underscores. The final placed file uses this name:

```text
<helix-name>_oriented_placed.pdb
```

### Output Folder

The output folder controls where the final placed PDB goes. The GUI also creates an intermediate folder below it:

```text
<output folder>/tmp_file/
```

The final placed PDB is meant to be the file you use downstream. The `tmp_file/` folder keeps DSSR tables, rebuilt PDB files, normalized files, minimized files, DSSR reports, and aligned-to-Z files so the workflow can be inspected later.

## Nucleic Acid Types

### B-DNA

B-DNA is built with:

```text
x3dna-dssr rebuild --backbone=B-DNA --par-type=heli
```

B-DNA defaults to running `phenix.geometry_minimization` after PDB name normalization.

### A-DNA

A-DNA is built with:

```text
x3dna-dssr rebuild --backbone=A-DNA --par-type=heli
```

A-DNA defaults to skipping `phenix.geometry_minimization`, but you can enable it.

### A-RNA

A-RNA is built with:

```text
x3dna-dssr rebuild --backbone=A-RNA --par-type=heli
```

A-RNA uses `U`, not `T`, in the input sequence. It defaults to skipping `phenix.geometry_minimization`, but you can enable it.

### Z-DNA

Z-DNA is built with:

```text
x3dna-dssr fiber --model=Z-DNA
```

The Z-DNA GUI path does not use the DSSR 12-parameter table and does not offer Phenix minimization.

## B-Z Structure Builder

The B-Z structure builder is based on the bundled `bnp_na_lib/make_BZV2_3.py` and `bnp_na_lib/core_BZ.py` scripts. This tool is different from the main `Generate` button: it does not read the sequence field and it does not create one isolated helix. Instead, it takes already-built B-DNA and Z-DNA PDB files and joins them through B-Z junction core structures.

In the main GUI, use the `B-Z builder` button directly after the four helix-type choices:

```text
B-Z builder
```

### Input Order

Add input PDB files in strict alternating order, starting with B-DNA:

```text
B1.pdb Z1.pdb B2.pdb Z2.pdb ...
```

The minimum input is two files:

```text
B1.pdb Z1.pdb
```

The builder assumes the type from the file order, not from the filename. A file named `my_z_model.pdb` in the first position is still treated as the B segment. The number of B models must be at least the number of Z models.

Each input helix should be a double-stranded DNA PDB with two main chains. The script picks the two chains with the largest residue counts.

### What The Builder Does

For every neighboring pair of input helices, the builder:

- Inserts one fresh copy of the bundled B-Z junction core.
- Fits the core to the previous helix using sugar atoms `C1'`, `C2'`, `C3'`, `C4'`, `O5'`, and `C5'`.
- Fits the next helix to the other side of the core.
- Trims overlapping base pairs at the helix ends while keeping the full junction core.
- Writes a final ligated PDB with two continuous DNA strands numbered 5' to 3'.
- Writes a raw aligned PDB with `_raw` before the file extension.

For example, if the output is:

```text
multi_BZ.pdb
```

the raw aligned file is:

```text
multi_BZ_raw.pdb
```

The final B-Z PDB gets `REMARK BNP_NA...` provenance records, including the `bnp_na` version and AZBMOST repository link. Since the B-Z builder does not itself make an L-form model, it writes `REMARK BNP_NA_L_FORM NO`.

### Axis Correction

The `Axis correction` option controls what happens after the local junction fit.

`codirectional` is the default. It rotates the moving helix so the fitted helix axes point in the same direction, but it does not laterally shift the helix to force both axes onto the exact same infinite line. This usually preserves the local junction geometry better.

`collinear` is stricter. It first makes the axes codirectional, then applies a lateral shift so the fitted axis lines coincide. This can make the global axis cleaner, but it may stretch bond lengths or local junction contacts. Because of that tradeoff, true collinear mode is useful to test but is not always the optimal final model compared with the default `codirectional` mode.

`none` skips post-fit axis correction and uses only the local sugar-atom junction fit.

### Axis Source

The `Axis source` option controls how the helix axis is estimated for the post-fit correction.

`auto` is the default. It tries DSSR `--more` axis extraction through `align2z.py`; if DSSR is not available for that axis calculation, it falls back to PCA/SVD on `C1'` atoms.

`dssr` requires DSSR axis extraction and stops if that fails.

`pca` uses only the internal `C1'` PCA axis estimate.

### Z-DNA Terminal Auto-Trim

For this builder, each selected Z-DNA chain should start with a pyrimidine (`C`, `T`, or `U`) and end with a purine (`A` or `G`) in residue-number order. This terminal phase keeps the Z-DNA end in register with the bundled B-Z core.

The GUI default is:

```text
Auto-trim terminal Z-DNA bp if needed: ON
```

When enabled, the builder can remove at most one base pair from each Z-DNA end to satisfy the terminal phase requirement. For Z-DNA generated by DSSR or by the `bnp_na` Z-DNA path, prepare the original Z-DNA input 2 bp longer than the target final Z segment so this correction can trim one base pair from each end if needed.

If you turn auto-trim off and a Z-DNA file has the wrong terminal phase, the build stops with an explicit error.

### Direct Command-Line Use

The bundled B-Z builder remains directly runnable:

```bash
python3 bnp_na_lib/make_BZV2_3.py --out multi_BZ.pdb B1.pdb Z1.pdb
```

Use codirectional axis correction explicitly:

```bash
python3 bnp_na_lib/make_BZV2_3.py \
  --axis-mode codirectional \
  --axis-source auto \
  --out multi_BZ.pdb \
  B1.pdb Z1.pdb B2.pdb
```

Use strict collinear correction:

```bash
python3 bnp_na_lib/make_BZV2_3.py \
  --axis-mode collinear \
  --out multi_BZ.pdb \
  B1.pdb Z1.pdb
```

Disable Z-DNA terminal auto-trim:

```bash
python3 bnp_na_lib/make_BZV2_3.py \
  --no-z-auto-trim \
  --out multi_BZ.pdb \
  B1.pdb Z1.pdb
```

## Triplex Converter

`bnp_na` V13.7 adds a triplex converter based on the bundled `bnp_na_lib/convert_to_triplex_pdbV2_1.py` script. This tool is different from the main `Generate` button: it reads an existing duplex DNA PDB and adds a third Hoogsteen or reverse-Hoogsteen strand over a selected residue range.

In the main GUI, use the `Triplex converter` button directly to the right of `B-Z builder`:

```text
Triplex converter
```

### Triplex Convention

The converter writes triplets as:

```text
Z·X-Y
```

`X-Y` is the Watson-Crick duplex base pair already present in the input PDB. `Z` is the third-strand base added by the converter.

The GUI and log use these strand names:

- Strand I: the purine duplex strand, `X` in `Z·X-Y`.
- Strand II: the Watson-Crick partner, `Y` in `Z·X-Y`.
- Strand III: the added Hoogsteen or reverse-Hoogsteen strand, `Z` in `Z·X-Y`.

### Supported Modes

Two embedded base-triple templates are bundled:

- `antiparallel`: `G·G-C`.
- `parallel`: `T·A-T`.

Only the selected residue range on strand I has to match the mode requirement. For `antiparallel`, the selected strand-I range must be all `G`. For `parallel`, the selected strand-I range must be all `A`. The rest of the input strand may contain other bases.

Strand II can be left blank. In that case, the converter searches other chains for a matching contiguous partner segment. Strand III can also be left blank; the converter then chooses an unused chain ID, preferring `C` when available.

### GUI Fields

`Input duplex PDB` is the existing duplex model to convert.

`Output PDB` defaults to the input filename with `_2TH` inserted before the extension. For example:

```text
duplex.pdb -> duplex_2TH.pdb
```

`Strand I purine chain` is the chain ID for the purine strand in the duplex.

`Strand I residue range` is the residue-number range on strand I to convert.

`Triplex mode` chooses the embedded base-triple template.

`Strand II chain` is optional. Use it when auto-detection chooses the wrong partner chain or you want to make the mapping explicit.

`Strand III chain` is optional. Use it when you need a specific chain ID for the added strand.

`Strand III first resSeq` controls the residue number assigned to the first residue of the added third strand.

`Refresh strand info` reads the input PDB, lists detected chains and sequences, and previews whether the selected strand-I range matches the chosen mode.

The final triplex PDB gets `REMARK BNP_NA...` provenance records, including the `bnp_na` version and AZBMOST repository link. Since the triplex converter does not itself make an L-form model, it writes `REMARK BNP_NA_L_FORM NO`.

### Direct Command-Line Use

The bundled triplex converter remains directly runnable:

```bash
python3 bnp_na_lib/convert_to_triplex_pdbV2_1.py duplex.pdb \
  --strand-I A \
  --range 10:18 \
  --mode antiparallel
```

Use parallel `T·A-T` mode and choose strand III chain/residue numbering:

```bash
python3 bnp_na_lib/convert_to_triplex_pdbV2_1.py duplex.pdb \
  --strand-I A \
  --range 7-12 \
  --mode parallel \
  --strand-III D \
  --strand-III-start 101 \
  --out duplex_triplex.pdb
```

## DSSR Parameter Customization

`Customize DSSR parameters` opens a table of 12 values used for B-DNA, A-DNA, and A-RNA DSSR rebuild jobs. The values are kept separately for each nucleic-acid type while the GUI is open.

The 12 columns are:

```text
Shear, Stretch, Stagger, Buckle, Propeller, Opening,
X-disp, Y-disp, h-Rise, Incl., Tip, h-Twist
```

The first six values are local base-pair parameters. The last six values are local helical-step parameters. In the generated DSSR table, the last row uses `999999` for the helical-step values because there is no next base pair after the final row.

Built-in default values:

| NA type | Shear | Stretch | Stagger | Buckle | Propeller | Opening | X-disp | Y-disp | h-Rise | Incl. | Tip | h-Twist |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| B-DNA | 0.0000 | -0.1500 | 0.0900 | 0.5000 | -11.4000 | 0.6000 | 0.0500 | 0.0200 | 3.4000 | 2.1000 | 0.0000 | 34.2857 |
| A-DNA | 0.0001 | -0.1448 | 0.0638 | 0.0003 | -10.5158 | -1.8170 | -4.4616 | 0.0001 | 2.5466 | 22.6460 | 0.0001 | 32.7273 |
| A-RNA | 0.0137 | -0.0848 | 0.0126 | -0.0044 | -2.0765 | -1.6676 | -4.0513 | 0.0678 | 2.8120 | 15.5148 | 0.7866 | 32.7273 |

Translational parameters are in Å. Angular parameters are in °. The app writes these values to the DSSR helical table with four digits after the decimal. If a field is left empty in the customization dialog, the corresponding default above is used.

Buttons in the parameter dialog:

- `Save` stores the current values for the selected nucleic-acid type.
- `Cancel` closes the dialog without applying the latest edits.
- `Restore defaults` fills the dialog with the built-in default values.
- `Clear fields` clears the text fields; empty fields mean "use the default value".

Only values that differ from the built-in defaults are reported as GUI overrides in the log.

## Phenix Geometry Minimization

The `phenix.geometry_minimization` section controls whether the normalized PDB is minimized before align-to-Z.

Defaults:

- B-DNA: minimization on.
- A-DNA: minimization off.
- A-RNA: minimization off.
- Z-DNA: minimization not available in the GUI.

The params file defaults to:

```text
bnp_na_lib/min_P_C5.params
```

The bundled default params file is intentionally small:

```text
pdb_interpretation {
    link_distance_cutoff = 7.0
}
selection = name " P " or name " OP1" or name " OP2" or name " O5'" or name " C5'" or name " O3'" or name " O1P" or name " O2P"
```

The `link_distance_cutoff = 7.0` line relaxes the bond/link-distance interpretation threshold used by Phenix while reading the generated nucleic-acid PDB. The `selection` line targets phosphate and nearby sugar-backbone atoms for the geometry-minimization run: `P`, `OP1`, `OP2`, `O5'`, `C5'`, and `O3'`. It also includes old phosphate atom names `O1P` and `O2P`, in case a generated or imported PDB still uses that convention.

When minimization is enabled, the app copies the params file into the intermediate output folder and runs Phenix there. The aligned input becomes the minimized PDB instead of the normalized PDB.

If Phenix is not found, either add `phenix.geometry_minimization` to `PATH` or set `PHENIX_ENV` to a Phenix environment script. The helper also checks the older default path:

```text
/Applications/phenix-1.21.2-5419/phenix_env.sh
```

## PDB Name Standardization

After DSSR creates the first PDB, the app runs `bnp_na_lib/pdb_name_standard.py`. This helper was changed from the previous script named `pdb_make_dna_v3_2.py`.

The standardization step does not build a new helix. It cleans the naming in the DSSR-generated PDB before optional Phenix minimization and final align-to-Z.

It can:

- Convert DNA residue names `A`, `T`, `C`, and `G` to `DA`, `DT`, `DC`, and `DG`.
- Convert standard three-letter nucleotide residue names such as `ADE`, `THY`, `CYT`, `GUA`, and `URI` to the expected canonical form.
- Preserve RNA residues when the residue contains `O2'` or old-style `O2*`.
- Rename old sugar atom `O2*` to `O2'`.
- Rename phosphate atoms `O1P` and `O2P` to `OP1` and `OP2`.
- Rename DT methyl atoms such as `C5M` and `H5M1` to the `C7`/`H7*` style.
- Delete hydrogen atoms when the GUI checkbox is enabled.

You can also run this helper directly:

```bash
python3 bnp_na_lib/pdb_name_standard.py input.pdb
python3 bnp_na_lib/pdb_name_standard.py input.pdb --deleteH
```

## Mirror-Image L-Form Modeling

The `Mirror-image L-form chirality` section can be used to generate mirror-image nucleic-acid models. This step runs after DSSR align-to-Z and before the final placement/orientation transform.

The order is:

```text
build -> normalize names -> optional minimization -> align to +Z -> optional inv/rot mirror -> orient/place
```

This position in the pipeline is important. Once the helix has been aligned to +Z, the coordinate axes are predictable. The mirror operation can change chirality first, and then the normal `roll`, `phi`, `theta`, `x`, `y`, and `z` controls can place the mirrored model.

Enable `Apply inv/rot after align-to-Z` to turn this on. When enabled, the operation dropdown becomes active.

### Why Inversion Plus Rotation Gives A Reflection

Changing chirality requires an improper transform, meaning a transform with determinant `-1`. A normal rotation has determinant `+1`, so it can rotate a model but cannot make its mirror image.

Point inversion maps:

```text
(x, y, z) -> (-x, -y, -z)
```

Point inversion has determinant `-1`, so it changes handedness. A 180-degree rotation has determinant `+1`. Combining point inversion with one or more 180-degree rotations still gives determinant `-1`, so the final operation is still a mirror/chirality-changing operation.

For example:

```text
inversion + 180-degree rotation around x
(-x, -y, -z) -> (-x, y, z)
```

That is equivalent to reflection across the `yz` plane.

### i Mode Operations

`i` mode means point inversion plus optional 180-degree rotations:

```text
i      inversion only
ix     inversion + 180-degree rotation around x
iy     inversion + 180-degree rotation around y
iz     inversion + 180-degree rotation around z
ixy    inversion + rotations around x and y
ixz    inversion + rotations around x and z
iyz    inversion + rotations around y and z
ixyz   inversion + rotations around x, y, and z
```

The original script also accepts instructions such as `x`, `xy`, or `xyz` as shorthand for `ix`, `ixy`, or `ixyz`. The GUI lists the explicit `i` forms.

### o Mode Operations

`o` mode names the reflection plane directly:

```text
oxy    reflection across the xy plane
oyz    reflection across the yz plane
oxz    reflection across the xz plane
```

Internally, these are implemented as inversion plus a 180-degree rotation around the perpendicular axis:

```text
oxy = inversion + Rz(180)
oyz = inversion + Rx(180)
oxz = inversion + Ry(180)
```

The GUI default is `oyz`. After align-to-Z, `oyz` and `oxz` keep the `z` coordinate sign, so the aligned +Z direction is preserved before placement. `oxy` and plain `i` change the `z` sign, so they reverse the aligned helix direction before placement.

When this option is enabled, the final placed PDB name includes the L-form operation label:

```text
<helix-name>_L_o_yz_oriented_placed.pdb
```

The intermediate mirrored PDB is written in:

```text
<output folder>/tmp_file/
```

The final placed PDB also contains machine-readable `REMARK` lines. These include provenance and the L-form residue annotations needed by future applications:

```text
REMARK BNP_NA bnp_na V13.8 from DiLiuLab's AZBMOST was used to create this file.
REMARK BNP_NA_REPOSITORY https://github.com/azbmost/bnp_na
REMARK BNP_NA_L_FORM YES
REMARK BNP_NA_L_FORM_KIND L-DNA
REMARK BNP_NA_L_RESIDUES_BEGIN COUNT <n>
REMARK BNP_NA_L_RESIDUE KIND L-DNA CHAIN <chain> RESSEQ <num> ICODE <icode> RESNAME <name>
REMARK BNP_NA_L_RESIDUES_END
```

For A-RNA mirror output, the residue kind is written as `L-RNA`. For ordinary non-mirrored output, the final PDB includes:

```text
REMARK BNP_NA_L_FORM NO
REMARK BNP_NA_L_RESIDUES NONE
```

You can also run the helper directly:

```bash
python3 bnp_na_lib/pdb_inv_rotV2.py model.pdb oyz
python3 bnp_na_lib/pdb_inv_rotV2.py model.pdb ix
```

## Helical-Axis Angle Tool

`bnp_na` includes `bnp_na_lib/angle_helical_axisV2_1.py`, an analysis tool for measuring how two points sit around a helical axis. This tool does not modify the model. It calculates radial vectors from a straight helical axis to two points, reports the angle between those radial directions, and writes a Chimera/ChimeraX `.bild` drawing.

The bundled filename is `angle_helical_axisV2_1.py` to indicate the V2.1 script update. Earlier public versions used `angle_helical_axisV2.py`. V13.5 added an adjustable axis drawing margin and more explanatory `.comment` records in the generated BILD file.

In the main `bnp_na` GUI, use the `Other tools` section near the bottom, immediately above `Log output`, and click:

```text
Measure angle around axis
```

This opens the angle tool in a separate window.

You can also launch it directly:

```bash
python3 bnp_na_lib/angle_helical_axisV2_1.py
python3 bnp_na_lib/angle_helical_axisV2_1.py --gui
```

### Axis Definition

The tool can define the helical axis in two ways.

`Fit from PDB` uses PCA/SVD on selected axis atoms from a PDB. The default axis atom name is `C1'`, with `C1*` accepted as an alias. You can provide one or more axis atom names in the GUI or with `--axis-atoms`.

`Custom axis` uses a point on the axis plus a direction vector. The direction vector is normalized automatically. In custom-axis mode, a PDB is not required unless one or both points are specified as atoms.

### Point Inputs

Point 1 and Point 2 each have a type selector:

- `Atom`: specify chain, residue number, and atom name, such as `A:5:C1'`.
- `XYZ`: specify explicit coordinates, such as `1 0 0` or `1,0,0`.

Atom point specifications require a PDB so the atom coordinates can be looked up. XYZ point specifications can be used without a PDB when the axis is also custom.

### Angle Definition

For each point, the tool projects the point onto the axis. The radial vector is:

```text
point - closest point on axis
```

The report includes:

- Distance from each point to the axis.
- The radial unit vector for each point.
- An unsigned angle from `0` to `180` degrees between the two radial unit vectors.
- A signed angle from `-180` to `180` degrees from Point 1 to Point 2 around the axis direction.
- The signed angle wrapped to `0` to `360` degrees.

The signed angle depends on axis direction. In PDB-fit mode, the tool tries to choose a consistent axis direction by correlating the fitted axis coordinate with residue numbering when enough selected atoms exist. In custom-axis mode, the direction is exactly the axis vector you provide after normalization.

This means custom-axis mode with XYZ points can also be used as a small geometry calculator.

### Command-Line Examples

Fit axis from PDB and use atom points:

```bash
python3 bnp_na_lib/angle_helical_axisV2_1.py -i helix.pdb \
  --point1 "A:5:C1*" \
  --point2 "B:18:C1*"
```

Use a custom axis and XYZ points without a PDB:

```bash
python3 bnp_na_lib/angle_helical_axisV2_1.py \
  --axis-point "0 0 0" --axis-vector "0 0 1" \
  --point1 "1 0 0" --point2 "0 1 0" \
  -o custom_axis_vectors.bild
```

Use a custom axis and atom points from a PDB:

```bash
python3 bnp_na_lib/angle_helical_axisV2_1.py -i helix.pdb \
  --axis-point "0 0 0" --axis-vector "0 0 1" \
  --point1 "A:5:C1*" --point2 "B:18:C1*" \
  -o helix_custom_axis_vectors.bild
```

Use a longer displayed axis arrow in the BILD drawing:

```bash
python3 bnp_na_lib/angle_helical_axisV2_1.py \
  --axis-point "0 0 0" --axis-vector "0 0 1" \
  --point1 "1 0 0" --point2 "0 1 0" \
  --axis-margin 12 \
  -o longer_axis_vectors.bild
```

### BILD Output

The tool writes Chimera/ChimeraX `.bild` files using `.comment`, `.color`, `.arrow`, and `.sphere` records. The BILD file draws:

- The fitted or custom helical axis as an arrow.
- The two radial vectors as arrows.
- Spheres at the two points.
- Spheres at the two point projections onto the axis.

The BILD file also includes `.comment` records immediately before the axis arrow, radial-vector arrows, point-marker spheres, and projection-marker spheres. These comments make the file easier to inspect or modify by hand.

Typical visualization workflow:

1. Open the PDB in Chimera or ChimeraX.
2. Open the generated `.bild` file.

Default drawing sizes are:

```text
--axis-margin 5.0
--axis-radius 1.0
--vector-radius 1.0
--sphere-radius 1.25
```

`--axis-margin` is in Angstrom. It controls how far the displayed axis arrow extends beyond the selected fitted range or the two custom-axis reference points. It changes only the BILD drawing length, not the fitted axis, custom axis, radial vectors, or reported angle. Increase it when the axis arrow looks too short in Chimera/ChimeraX; decrease it when the arrow visually dominates a small local measurement.

The PCA/SVD axis fit uses `numpy.linalg.svd`; NumPy's SVD reference is here: <https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html>.

## Align Helix To Z Tool

`bnp_na` V13.8 includes an `Align helix to z` button in the main GUI's `Other tools` section. This tool reads an existing helix PDB, runs `x3dna-dssr --more` to obtain DSSR `point-one` and `point-two` axis endpoints, translates `point-one` to the origin, and rotates the helix axis onto the coordinate-system +Z axis.

The dialog asks for:

- Input PDB file.
- Aligned output PDB file.
- DSSR `.out` file for the `--more` report.

The tool appends the equivalent `bnp_na_lib/align2z.py` CLI command and alignment report to the main log.

## Helical-Axis Info Tool

`bnp_na` V13.8 includes `bnp_na_lib/helical_axis_info.py`, a DSSR-based tool for reporting the helical axis formed by two selected chains in a larger PDB file. The tool filters the input PDB down to the requested chain IDs, runs `x3dna-dssr --more`, and parses DSSR's first `point-one` and `point-two` axis endpoints.

In the main `bnp_na` GUI, use the `Other tools` section near the bottom, immediately above `Log output`, and click:

```text
Get helical-axis info
```

The dialog asks for:

- Input PDB file.
- Two PDB chain IDs, such as `C D` or `C,D`.
- A reference vector for angle reporting, default `0 0 1`.
- Optional helix length in base pairs for estimating the full helix length from the DSSR endpoint distance.
- Optional Chimera/ChimeraX `.bild` output for the helical axis.
- Optional control for whether the reference vector is drawn in the `.bild` file, plus its drawing length.

The report includes:

- DSSR helical-axis start point and end point.
- Axis vector, start-to-end distance, and normalized unit vector.
- Reference-vector unit vector.
- Unsigned angle from `0` to `180` degrees between the DSSR axis direction and the reference vector.
- A note that the start-to-end distance spans the start and end base-pair centers, so it is approximately one base-pair step shorter than the full helix length.
- If helix length is provided as `n` bp, an estimated full helix length using `distance / (n - 1) * n`.
- Paths for the selected-chain PDB, DSSR `.out` file, and optional BILD output.

The order of the two chain IDs defines the reported axis direction. For example, `A B` and `B A` report opposite axis vectors and opposite unit vectors for the same duplex.

The selected-chain PDB and DSSR `.out` file are written to the selected output folder's `tmp_file/` folder. The DSSR output filename ends with `.out`, for example:

```text
model_chains_CD_dssr_more.out
```

You can also run it directly:

```bash
python3 bnp_na_lib/helical_axis_info.py -i model.pdb \
  --chains "C D" \
  --vector "0 0 1" \
  --helix-bp 24 \
  --workdir output/tmp_file \
  --bild model_chains_CD_axis.bild \
  --reference-vector-length 20
```

## XYZ Axes BILD Tool

`bnp_na` includes `bnp_na_lib/xyz_bild.py`, a small utility for writing a coordinate-axis `.bild` file for Chimera or ChimeraX. It draws:

- A sphere at the origin.
- A red X-axis arrow.
- A yellow Y-axis arrow.
- A blue Z-axis arrow.

In the main `bnp_na` GUI, use the `Other tools` section near the bottom, immediately above `Log output`, and click:

```text
Write XYZ axes BILD
```

The dialog lets you choose the output file and set:

- `Origin x y z (Å)`, default `0 0 0`.
- `Arrow length (Å)`, default `20`.
- `Arrow width (Å)`, default `1`.
- `Origin sphere radius (Å)`, default `0.5`.

The arrow width is used as the BILD arrow shaft radius. The arrow head radius is `2.5 x arrow width`.

The default output is equivalent to this BILD geometry, with explicit arrow-width parameters added:

```text
.sphere 0 0 0 0.5
.color 1 0 0
.arrow 0 0 0 20 0 0
.color 1 1 0
.arrow 0 0 0 0 20 0
.color 0 0 1
.arrow 0 0 0 0 0 20
```

You can also run it directly:

```bash
python3 bnp_na_lib/xyz_bild.py -o xyz_axes.bild
python3 bnp_na_lib/xyz_bild.py -o xyz_axes.bild --length 30 --width 0.75
python3 bnp_na_lib/xyz_bild.py -o shifted_axes.bild --origin "10 0 0" --length 15 --width 0.5
```

## Delete Hydrogens

`Delete hydrogens from generated PDB` removes hydrogen atoms during the PDB naming-normalization step. This is useful when downstream refinement, minimization, or model placement should start from heavy atoms only.

This option affects B-DNA, A-DNA, A-RNA, and Z-DNA after the initial DSSR build/fiber step.

## Placement And Orientation

Before final placement, the app aligns the generated helix so that its axis starts at the origin and points along +Z. The placement controls then transform that aligned model.

Units:

- `x`, `y`, `z`, and `delta_z` are in Å.
- `roll`, `phi`, and `theta` are in °.

Transform order:

1. Shift along local +Z by `delta_z`.
2. Roll about the local +Z axis by `roll`.
3. Rotate by `phi` around the Y axis.
4. Rotate by `theta` around the Z axis.
5. Translate by `x`, `y`, and `z`.

`delta_z` should usually stay `0`. Use it only when you need to move the helix along its own aligned axis before the angular rotations and final translation.

`roll` changes rotation around the helix axis. The GUI note says:

```text
For GIDEON: roll = roll at GIDEON - 111.25
```

This means that if you are copying a roll angle from GIDEON, subtract `111.25` before entering it here.

`phi` tilts the helix away from +Z through a Y-axis rotation.

`theta` turns the tilted helix around the global Z axis.

`x`, `y`, and `z` are the final translation applied after all rotations.

## Log Output

The embedded log is not just a message box. It records the build recipe and is the best place to debug a failed run.

The log includes:

- DSSR startup check output.
- The selected nucleic-acid type.
- Output and intermediate folders.
- Whether Phenix minimization was enabled.
- Whether hydrogens were deleted.
- Exact DSSR and Phenix command lines.
- Paths to generated tables, rebuilt PDB files, normalized PDB files, minimized PDB files, aligned PDB files, and final placed PDB files.
- The L-form inv/rot mirror operation and mirrored intermediate PDB path when enabled.
- The final PDB `REMARK` annotations, including provenance and L-form residue records.
- Parameter overrides applied from the GUI.
- The final orientation/placement transform.
- B-Z builder input order, axis mode/source, Z-DNA auto-trim status, junction-fit messages, final B-Z PDB path, and raw aligned PDB path.
- Triplex converter chain/sequence preview, selected strand-I range check, mode, strand IDs, template-fit RMSD summary, and final triplex PDB path.

## Generated Files

For a helix named `B-DNA25`, typical intermediate files include:

```text
B-DNA25.txt
B-DNA25-rb10.5.pdb
B-DNA25-rb10.5_out.pdb
B-DNA25-rb10.5_out_minimized.pdb
B-DNA25-rb10.5_out_minimized_aligned2Z.pdb
dssr-frames.txt
dssr-helices.pdb
dssr-pairs.pdb
```

For B-DNA, the `rb` number in intermediate filenames is the helical repeat
computed from the effective twist parameter: `360 / h-Twist`, rounded to at
most two decimal places. With the default B-DNA `h-Twist = 34.2857`, this is
`rb10.5`; with `h-Twist = 30`, future generated files use `rb12`.

Not every file appears for every run. For example, minimized files appear only when Phenix minimization is enabled.

The final placed model is written outside `tmp_file/`:

```text
<output folder>/<helix-name>_oriented_placed.pdb
```

B-Z builder outputs are written to the path chosen in the B-Z dialog. They are not placed in `tmp_file/` unless you explicitly choose that folder. The paired raw aligned file is written beside the final B-Z PDB with `_raw` added to the filename.

## Common Problems

### `x3dna-dssr: NOT FOUND`

Install DSSR and make sure the executable is named `x3dna-dssr` on `PATH`, or place it at:

```text
/usr/local/bin/x3dna-dssr
```

### DSSR Rebuild Failed

Check that the sequence uses the correct alphabet for the selected nucleic-acid type. DNA uses `T`; RNA uses `U`. Also inspect the generated `.txt` helical table in `tmp_file/`.

### Phenix Minimization Failed

Check that Phenix is installed and that the params file exists. If Phenix is installed but not on `PATH`, set `PHENIX_ENV` to the environment script before launching the GUI.

### Invalid Z-DNA Length

Z-DNA length must be positive and even. Examples that work: `10`, `20`, `42`. Examples that fail: `0`, `-2`, `15`, `abc`.

### B-Z Builder Input Order Or Z-DNA Phase Failed

The B-Z builder assumes the input order is B, Z, B, Z, starting with B. It does not infer B or Z from filenames. Check the file list in the B-Z dialog before running.

If the error mentions Z-DNA terminal phase, turn on `Auto-trim terminal Z-DNA bp if needed`, or regenerate the Z-DNA input 2 bp longer than the target final Z segment. The builder needs selected Z chains to start with a pyrimidine and end with a purine.

### Triplex Converter Range Or Partner Strand Failed

For `antiparallel` mode, the selected strand-I range must be all `G`. For `parallel` mode, the selected strand-I range must be all `A`. Use `Refresh strand info` to preview the selected sequence before running conversion.

If strand-II auto-detection fails, specify the `Strand II chain` explicitly. The selected partner segment must be a contiguous `C` segment for antiparallel `G·G-C` or a contiguous `T` segment for parallel `T·A-T`.

### The Final PDB Is Not Where Expected

The final placed PDB goes directly in the selected output folder. Intermediate files go in `tmp_file/` below that folder.

## Make The Script Directly Executable

The script already has a Python shebang. On macOS or Linux, mark it executable once:

```bash
chmod +x bnp_na.py
./bnp_na.py
```

This still uses your local Python environment and installed dependencies.

## Build A Standalone Executable

Install PyInstaller in the environment where NumPy is available:

```bash
python3 -m pip install pyinstaller
```

Build a single-folder app on macOS or Linux:

```bash
python3 -m PyInstaller --name bnp_na --windowed --paths bnp_na_lib --add-data "bnp_na_lib:bnp_na_lib" --add-data "assets:assets" bnp_na.py
```

On Windows, use semicolons in `--add-data`:

```powershell
python -m PyInstaller --name bnp_na --windowed --paths bnp_na_lib --add-data "bnp_na_lib;bnp_na_lib" --add-data "assets;assets" --icon assets/bnp_na_icon.ico bnp_na.py
```

The GUI checks for `assets/bnp_na_icon.png` at startup and continues normally if the icon asset is missing. The tracked icon uses a right-handed helix and a +Z cue for the build-and-place workflow.

## Repository Layout

```text
bnp_na.py                  Main GUI/controller
CHANGELOG.md               Version-by-version change log
bnp_na_lib/                Build, alignment, placement, and PDB helpers
bnp_na_lib/angle_helical_axisV2_1.py Helical-axis radial-angle calculator and BILD writer
bnp_na_lib/build_bz.py      bnp_na wrapper for the B-Z junction builder
bnp_na_lib/build_triplex.py bnp_na wrapper for the triplex converter
bnp_na_lib/helical_axis_info.py DSSR selected-chain helical-axis reporter and BILD writer
bnp_na_lib/make_BZV2_3.py   Standalone B-Z structure builder CLI/GUI
bnp_na_lib/core_BZ.py       Bundled B-Z junction core structure data
bnp_na_lib/convert_to_triplex_pdbV2_1.py Standalone duplex-to-triplex converter CLI/GUI
bnp_na_lib/pdb_inv_rotV2.py Optional inversion/reflection helper for L-form mirror models
bnp_na_lib/pdb_name_standard.py PDB residue/atom-name standardization helper
bnp_na_lib/xyz_bild.py      Coordinate-axis BILD writer
bnp_na_lib/min_P_C5.params Default Phenix minimization params file
assets/                    Optional app/taskbar icon assets
requirements.txt           Python package dependency list
```

Generated run outputs are intentionally ignored by git:

```text
output/
tmp_file/
test/
```

## License

This project is released under the MIT License. See `LICENSE`.
