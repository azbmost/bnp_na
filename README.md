# bnp_na

`bnp_na` is a Tkinter GUI for building and placing nucleic acid helices. It can generate B-DNA, A-DNA, A-RNA, and Z-DNA models, normalize PDB atom/residue names, align the helix to the +Z axis, and write a final oriented/placed PDB file.

The current app version is `V13.1`.

## What It Does

- Builds B-DNA, A-DNA, and A-RNA from a 5' to 3' sequence using DSSR helical-parameter tables.
- Builds Z-DNA from a positive even base-pair length using DSSR fiber generation.
- Lets you customize the 12 DSSR base-pair/helical parameters for B-DNA, A-DNA, and A-RNA.
- Optionally runs `phenix.geometry_minimization` for B-DNA, A-DNA, and A-RNA.
- Aligns the generated helix to +Z, then applies roll, phi, theta, x, y, z, and delta_z placement values.
- Writes final PDB files to the selected output folder and intermediate files to `<output folder>/tmp_file/`.

## Requirements

- Python 3.9 or newer.
- Tkinter, normally included with Python from python.org and many system Python installs.
- NumPy, installed with `pip install -r requirements.txt`.
- `x3dna-dssr` available on `PATH`, or installed at `/usr/local/bin/x3dna-dssr`.
- Optional: `phenix.geometry_minimization` on `PATH`, or `PHENIX_ENV` pointing to a valid Phenix environment script.

## Clone The Repository

`git clone` downloads a complete local copy of the repository:

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

## Pull Updates

After cloning once, use `git pull` inside the repository folder to fetch and merge the latest changes from GitHub into your local copy:

```bash
cd bnp_na
git pull
```

If you have local edits, commit or stash them before pulling so Git can merge cleanly.

## Basic Use

1. Enter a 5' to 3' sequence for B-DNA, A-DNA, or A-RNA. Compact syntax such as `A10T5C2G` is supported for DNA and `A10U5C2G` for RNA.
2. For Z-DNA, choose `Z-DNA` and enter an even helix length. The sequence field is ignored for Z-DNA.
3. Choose an output folder. Final placed PDB files are written there, and intermediate files go to `tmp_file/` below it.
4. Adjust DSSR parameters if needed.
5. Turn `phenix.geometry_minimization` on or off. B-DNA defaults to on; A-DNA and A-RNA default to off.
6. Set placement/orientation values, then click `Generate`.

## Output Files

The final placed PDB is named:

```text
<helix-name>_oriented_placed.pdb
```

Intermediate DSSR, normalized, minimized, and aligned files are placed in:

```text
<output folder>/tmp_file/
```

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

Build a single-folder app:

```bash
python3 -m PyInstaller --name bnp_na --windowed --paths bnp_na_lib --add-data "bnp_na_lib:bnp_na_lib" --add-data "assets:assets" bnp_na.py
```

On Windows, use semicolons in `--add-data`:

```powershell
python -m PyInstaller --name bnp_na --windowed --paths bnp_na_lib --add-data "bnp_na_lib;bnp_na_lib" --add-data "assets;assets" --icon assets/bnp_na_icon.ico bnp_na.py
```

If you do not want to include an icon, omit the `--add-data "assets:assets"` or `--icon` option. The GUI checks for `assets/bnp_na_icon.png` at startup and continues normally if the icon asset is missing.

## Repository Layout

```text
bnp_na.py                 Main GUI/controller
bnp_na_lib/               Build, alignment, placement, and PDB helpers
bnp_na_lib/min_P_C5.params Default Phenix minimization params file
assets/                   Optional app/taskbar icon assets
```

## License

This project is released under the MIT License. See `LICENSE`.
