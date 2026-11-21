#!/usr/bin/env python3
"""bnp_na V13.23: Building and placing nucleic acid helices.

Top-level GUI/controller. All helper modules live in ./bnp_na_lib/.
"""
from __future__ import annotations

import ast
import math
import queue
import re
import shlex
import sys
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__version__ = "V13.23"
APP_NAME = "bnp_na"

HELP_TEXT = f"""usage: {APP_NAME}.py [-h] [-v]

{APP_NAME} {__version__} - build and place nucleic acid helices.

Runs as a Tkinter GUI. Start it with no options:

  python3 {APP_NAME}.py

options:
  -h, --help     show this message and exit
  -v, --version  print the version and exit

The GUI builds B-DNA, A-DNA, A-RNA, Z-DNA, and P-DNA helices, and hosts the
B-Z builder, the B-P builder, the triplex converter, and the tools under
`Other tools`.

Many of those tools also run standalone from the helper folder, each with
its own --help, for example:

  python3 bnp_na_lib/combine_pdb.py --help
  python3 bnp_na_lib/regularize_phosphates.py --help
  python3 bnp_na_lib/build_pdna.py --help
  python3 bnp_na_lib/make_BP.py --help

x3dna-dssr must be on PATH for model generation; phenix.geometry_minimization
is optional. See README.md for the full guide."""


def _answer_cli_request(args: List[str]) -> None:
    """Answer -h/--help and -v/--version before the GUI toolkit is imported.

    Staying ahead of the tkinter import keeps both working on systems with no
    Tk installed, and stops --help from silently opening a window instead of
    printing anything. The GUI itself takes no arguments, so anything else is
    a usage error rather than a reason to launch.
    """
    if not args:
        return
    if any(arg in ("-h", "--help") for arg in args):
        print(HELP_TEXT)
        sys.exit(0)
    if any(arg in ("-v", "--version") for arg in args):
        print(f"{APP_NAME} {__version__}")
        sys.exit(0)
    print(f"{APP_NAME}: unrecognized argument: {args[0]}", file=sys.stderr)
    print(f"Try 'python3 {APP_NAME}.py --help' for usage.", file=sys.stderr)
    sys.exit(2)


if __name__ == "__main__":
    _answer_cli_request(sys.argv[1:])

import tkinter as tk  # noqa: E402
from tkinter import filedialog, messagebox, ttk  # noqa: E402

APP_DIR = Path(__file__).resolve().parent
LIB_DIR = APP_DIR / "bnp_na_lib"

from bnp_na_lib.build_adna import build_adna  # noqa: E402
from bnp_na_lib.build_arna import build_arna  # noqa: E402
from bnp_na_lib.build_bp import build_bp_structure  # noqa: E402
from bnp_na_lib.build_bz import build_bz_structure  # noqa: E402
from bnp_na_lib.build_bdna import build_bdna  # noqa: E402
from bnp_na_lib.build_pdna import build_pdna, resolve_pdna_sequences  # noqa: E402
from bnp_na_lib import core_PDNA  # noqa: E402
from bnp_na_lib.build_common import (  # noqa: E402
    DEFAULT_PARAMS,
    PARAM_KEYS,
    PipelineError,
    check_dssr_installation,
    expand_sequence,
    sanitize_basename,
    sequence_alphabet,
)
from bnp_na_lib.build_triplex import (  # noqa: E402
    TMP_DIR_NAME as TRIPLEX_TMP_DIR_NAME,
    build_triplex_from_duplex,
    default_triplex_output_path,
    describe_triplex_input,
)
from bnp_na_lib.build_zdna import build_zdna  # noqa: E402
from bnp_na_lib.align2z import align_pdb_to_Z, format_alignment_report  # noqa: E402
from bnp_na_lib.combine_pdb import (  # noqa: E402
    combine_pdb_files,
    default_combine_pdb_output_path,
    format_chain_id as format_combine_chain_id,
    format_chain_token as format_combine_chain_token,
    list_pdb_chains,
    parse_chain_selection as parse_combine_chain_selection,
)
from bnp_na_lib.add_phosphates import (  # noqa: E402
    add_terminal_phosphates,
    analyze_phosphate_termini,
    default_add_phosphate_output_path,
    format_phosphate_report,
    parse_chain_selection,
)
from bnp_na_lib.regularize_phosphates import (  # noqa: E402
    default_regularized_output_path,
    regularize_phosphates,
)
from bnp_na_lib.opposing_phosphate_xdisp import (  # noqa: E402
    OpposingPhosphateCancelled,
    find_opposing_phosphate_xdisp,
    format_result_report,
)
from bnp_na_lib.angle_helical_axisV2_2 import launch_gui as launch_axis_angle_gui  # noqa: E402
from bnp_na_lib.helical_axis_info import format_axis_info_report, get_helical_axis_info, parse_chain_ids  # noqa: E402
from bnp_na_lib.pdb_inv_rotV2 import InvRotError, apply_inv_rot_to_pdb, parse_operation  # noqa: E402
from bnp_na_lib.xyz_bild import write_xyz_bild  # noqa: E402
from bnp_na_lib.na_placer import PlacerError, place_after_Z  # noqa: E402
from bnp_na_lib import azbmost_ui as ui  # noqa: E402
from bnp_na_lib import bnp_na_icons  # noqa: E402,F401  (registers the B-P builder glyph)


PARAM_BASE_LABELS = {
    "Shear": "Shear",
    "Stretch": "Stretch",
    "Stagger": "Stagger",
    "Buckle": "Buckle",
    "Propeller": "Propeller",
    "Opening": "Opening",
    "X-disp": "X-disp",
    "Y-disp": "Y-disp",
    "h-Rise": "h-Rise",
    "Incl.": "Incl.",
    "Tip": "Tip",
    "h-Twist": "h-Twist",
}
PARAM_UNITS = {
    "Shear": "Å",
    "Stretch": "Å",
    "Stagger": "Å",
    "Buckle": "°",
    "Propeller": "°",
    "Opening": "°",
    "X-disp": "Å",
    "Y-disp": "Å",
    "h-Rise": "Å",
    "Incl.": "°",
    "Tip": "°",
    "h-Twist": "°",
}
PARAM_LABELS = {
    key: f"{PARAM_BASE_LABELS.get(key, key)} ({PARAM_UNITS[key]})" for key in PARAM_BASE_LABELS
}

DEFAULT_PARAMS_FILE = LIB_DIR / "min_P_C5.params"
DEFAULT_OUTPUT_DIR = APP_DIR / "output"
DEFAULT_ICON_FILE = APP_DIR / "assets" / "bnp_na_icon.png"
NA_TYPES = ("B-DNA", "A-DNA", "A-RNA", "Z-DNA", "P-DNA")
NA_TYPES_WITH_TABLE = ("B-DNA", "A-DNA", "A-RNA")
DEFAULT_MINIMIZE_BY_TYPE = {"B-DNA": True, "A-DNA": False, "A-RNA": False, "P-DNA": False}
# Shown read-only for P-DNA: the 1PFI template fixes them (core_PDNA.py).
PDNA_PARAM_COLUMNS = (
    ("h-Rise (Å/nt)", core_PDNA.RISE),
    ("h-Twist (°/nt)", core_PDNA.TWIST_DEG),
    ("nt / turn", core_PDNA.NT_PER_TURN),
    ("Pitch (Å)", core_PDNA.PITCH),
    ("P radius (Å)", core_PDNA.P_RADIUS),
    ("C1' radius (Å)", core_PDNA.C1_RADIUS),
)
PDNA_PARAM_NOTE = (
    "P-DNA (Pf1 phage DNA, PDB 1PFI) is built from the template's backbone: phosphates at the center, "
    "bases pointing outward, no base pairs. Its helical parameters are fixed and cannot be customized."
)
INV_ROT_OPERATIONS = ("oyz", "oxz", "oxy", "i", "ix", "iy", "iz", "ixy", "ixz", "iyz", "ixyz")

_ARITHMETIC_BIN_OPS = {
    ast.Add: lambda a, b: a + b,
    ast.Sub: lambda a, b: a - b,
    ast.Mult: lambda a, b: a * b,
    ast.Div: lambda a, b: a / b,
    ast.Pow: lambda a, b: a**b,
    ast.Mod: lambda a, b: a % b,
}
_ARITHMETIC_UNARY_OPS = {
    ast.UAdd: lambda a: a,
    ast.USub: lambda a: -a,
}


def _path_text(path: Path) -> str:
    try:
        return str(path.resolve())
    except Exception:
        return str(path)


def _eval_arithmetic_node(node: ast.AST) -> float:
    if isinstance(node, ast.Expression):
        return _eval_arithmetic_node(node.body)
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return float(node.value)
    if isinstance(node, ast.BinOp):
        op_type = type(node.op)
        if op_type not in _ARITHMETIC_BIN_OPS:
            raise ValueError("unsupported operator")
        left = _eval_arithmetic_node(node.left)
        right = _eval_arithmetic_node(node.right)
        return float(_ARITHMETIC_BIN_OPS[op_type](left, right))
    if isinstance(node, ast.UnaryOp):
        op_type = type(node.op)
        if op_type not in _ARITHMETIC_UNARY_OPS:
            raise ValueError("unsupported unary operator")
        return float(_ARITHMETIC_UNARY_OPS[op_type](_eval_arithmetic_node(node.operand)))
    raise ValueError("only numbers and +, -, *, /, %, **, and parentheses are allowed")


def _parse_float_expression(text: str, field_name: str, *, default: Optional[float] = None) -> float:
    expr = text.strip()
    if not expr:
        if default is not None:
            return float(default)
        raise ValueError(f"{field_name} is blank.")
    try:
        value = _eval_arithmetic_node(ast.parse(expr, mode="eval"))
    except Exception as exc:
        raise ValueError(f"{field_name} must be a number or simple arithmetic expression, for example 360/10.5.") from exc
    if not math.isfinite(value):
        raise ValueError(f"{field_name} must evaluate to a finite number.")
    return float(value)


def _parse_int_expression(text: str, field_name: str) -> int:
    value = _parse_float_expression(text, field_name)
    rounded = round(value)
    if abs(value - rounded) > 1e-9:
        raise ValueError(f"{field_name} must evaluate to an integer.")
    return int(rounded)


def _format_number(value: float) -> str:
    return f"{float(value):.4f}"


def _ensure_output_dirs(output_dir_text: str) -> Tuple[Path, Path]:
    if not output_dir_text.strip():
        raise ValueError("Please choose an output folder.")
    output_dir = Path(output_dir_text).expanduser()
    if not output_dir.is_absolute():
        output_dir = output_dir.resolve()
    tmp_dir = output_dir / "tmp_file"
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return output_dir, tmp_dir


def _final_placed_path(output_dir: Path, base_name: str, suffix: str = "") -> Path:
    safe = sanitize_basename(base_name or "bnp_na_helix") or "bnp_na_helix"
    safe_suffix = sanitize_basename(suffix)
    suffix_part = f"_{safe_suffix}" if safe_suffix else ""
    return output_dir / f"{safe}{suffix_part}_oriented_placed.pdb"


def _pdb_coord_residues(pdb_path: Path) -> list[Tuple[str, str, str, str]]:
    residues = []
    seen = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip() or "."
            chain = line[21].strip() or "."
            resseq = line[22:26].strip() or "."
            icode = line[26].strip() or "."
            key = (chain, resseq, icode, resname)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def _prepend_final_pdb_remarks(
    pdb_path: Path,
    na_type: str,
    l_form_enabled: bool,
    invrot_result: Optional[Dict[str, object]],
) -> str:
    residues = _pdb_coord_residues(pdb_path)
    l_kind = "L-RNA" if na_type == "A-RNA" else "L-DNA"
    lines = [
        f"REMARK BNP_NA bnp_na {__version__} from DiLiuLab's AZBMOST was used to create this file.",
        "REMARK BNP_NA_REPOSITORY https://github.com/azbmost/bnp_na",
        f"REMARK BNP_NA_NA_TYPE {na_type}",
    ]

    if l_form_enabled:
        label = str(invrot_result.get("label", "")) if invrot_result else ""
        mode = str(invrot_result.get("mode", "")) if invrot_result else ""
        axes = str(invrot_result.get("rotation_axes", "")) if invrot_result else ""
        signs = invrot_result.get("coordinate_signs") if invrot_result else None
        signs_text = ",".join(str(value) for value in signs) if signs else ""
        lines += [
            "REMARK BNP_NA_L_FORM YES",
            f"REMARK BNP_NA_L_FORM_KIND {l_kind}",
            f"REMARK BNP_NA_INV_ROT_LABEL {label}",
            f"REMARK BNP_NA_INV_ROT_MODE {mode}",
            f"REMARK BNP_NA_INV_ROT_AXES {axes if axes else 'none'}",
            f"REMARK BNP_NA_INV_ROT_SIGNS {signs_text}",
            f"REMARK BNP_NA_L_RESIDUES_BEGIN COUNT {len(residues)}",
        ]
        for chain, resseq, icode, resname in residues:
            lines.append(
                f"REMARK BNP_NA_L_RESIDUE KIND {l_kind} CHAIN {chain} "
                f"RESSEQ {resseq} ICODE {icode} RESNAME {resname}"
            )
        lines.append("REMARK BNP_NA_L_RESIDUES_END")
    else:
        lines += [
            "REMARK BNP_NA_L_FORM NO",
            "REMARK BNP_NA_L_RESIDUES NONE",
        ]

    original = pdb_path.read_text(encoding="utf-8", errors="ignore")
    remark_block = "\n".join(lines) + "\n"
    pdb_path.write_text(remark_block + original, encoding="utf-8")
    return "\n".join(["=== Final PDB remarks ===", *lines])


class App(tk.Tk):
    def __init__(self):
        super().__init__(baseName=APP_NAME, className=APP_NAME)
        self._set_app_identity()
        self.title(f"{APP_NAME} {__version__} - AZBMOST Package Module #1 - Build and Place Nucleic Acid")
        self._set_optional_window_icon()
        ui.apply_theme(self, "bnp_na")
        self.geometry("1440x1000")
        self.minsize(1120, 700)
        self._style = ttk.Style(self)
        # Style names the dialogs use; they inherit the shared AZBMOST theme.
        self._style.configure("Compact.Treeview", rowheight=26, font="AzbSmall", borderwidth=1,
                              bordercolor=ui.C.BORDER, lightcolor=ui.C.BG, darkcolor=ui.C.BG)
        self._style.configure("Compact.Treeview.Heading", font="AzbSmallBold", padding=(4, 5))

        self.param_values_by_type: Dict[str, Dict[str, str]] = {
            na_type: {key: "" for key in PARAM_KEYS} for na_type in NA_TYPES_WITH_TABLE
        }
        self.minimize_by_type: Dict[str, bool] = dict(DEFAULT_MINIMIZE_BY_TYPE)
        self.regularize_phosphates_by_type: Dict[str, bool] = dict(DEFAULT_MINIMIZE_BY_TYPE)
        self._active_na_type: Optional[str] = None

        field = {"pady": 4}
        label_pad = (0, 14)

        # --- Header ------------------------------------------------------------
        header = ui.Header(
            self,
            "bnp_na",
            version=__version__,
            subtitle="Build and place nucleic-acid helices  ·  AZBMOST module 1 of 3",
        )
        header.pack(side="top", fill="x")
        self.dssr_status_var = tk.StringVar(value="x3dna-dssr: checking after GUI starts...")
        self.dssr_chip = ui.StatusChip(header.right, text="x3dna-dssr: checking after GUI starts...", kind="busy",
                                       style="Header.TFrame")
        self.dssr_chip.pack(anchor="e")
        self.dssr_path_label = ttk.Label(header.right, text="", style="Hint.TLabel")
        self.dssr_path_label.pack(anchor="e", pady=(2, 0))
        ui.hairline(self).pack(side="top", fill="x")

        # --- Footer: run buttons -------------------------------------------------
        footer = ui.Footer(self, text="Generate builds the helix, aligns its axis to +Z, then places it.")
        footer.pack(side="bottom", fill="x")
        ui.hairline(self).pack(side="bottom", fill="x")
        footer.add_button("Quit", self.destroy)
        self.generate_btn = footer.add_button("Generate", self.on_generate, primary=True)

        # --- Body: form on the left, tools and log on the right ------------------
        body = ttk.Frame(self)
        body.pack(side="top", fill="both", expand=True)
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1, minsize=640)
        form_area = ui.ScrollFrame(body)
        form_area.grid(row=0, column=0, sticky="nsew")
        ui.hairline(body, orient="vertical").grid(row=0, column=1, sticky="ns")
        side = ttk.Frame(body, padding=(14, 12, 20, 14), width=440)
        side.grid(row=0, column=2, sticky="nsew")
        side.grid_propagate(False)
        form = form_area.body
        form.configure(padding=(20, 16, 18, 20))
        form.columnconfigure(0, weight=1, uniform="form")
        form.columnconfigure(1, weight=1, uniform="form")

        # Helix ----------------------------------------------------------------
        helix = ui.Card(form, "Helix")
        helix.grid(row=0, column=0, columnspan=2, sticky="we", pady=(0, 12))
        helix.columnconfigure(1, weight=1)

        ttk.Label(helix, text="Nucleic acid type", style="Field.TLabel").grid(
            row=0, column=0, sticky="w", padx=label_pad, **field
        )
        self.na_type_var = tk.StringVar(value="B-DNA")
        self.type_selector = ui.Segmented(helix, self.na_type_var, list(NA_TYPES), command=self.on_type_changed)
        self.type_selector.grid(row=0, column=1, columnspan=2, sticky="w", **field)

        self.seq_label = ttk.Label(helix, text="Sequence (5'->3'):", style="Field.TLabel")
        self.seq_label.grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        self.seq_var = tk.StringVar(value="")
        self.seq_entry = ttk.Entry(helix, textvariable=self.seq_var)
        self.seq_entry.grid(row=1, column=1, sticky="we", **field)
        self.length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(helix, textvariable=self.length_var, style="Hint.TLabel").grid(
            row=1, column=2, sticky="w", padx=(12, 0), **field
        )

        # P-DNA only: strand II has its own sequence, as P-DNA forms no base pairs.
        self.seq2_label = ttk.Label(helix, text="Strand II (5'->3', optional):", style="Field.TLabel")
        self.seq2_var = tk.StringVar(value="")
        self.seq2_entry = ttk.Entry(helix, textvariable=self.seq2_var)
        self.seq2_hint = ttk.Label(
            helix,
            text=(
                "Blank = reverse complement of strand I, which puts complementary bases on the same level. "
                "Both strands need the same length."
            ),
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        )
        self.seq2_label.grid(row=2, column=0, sticky="w", padx=label_pad, **field)
        self.seq2_entry.grid(row=2, column=1, sticky="we", **field)
        self.seq2_hint.grid(row=3, column=1, columnspan=2, sticky="w", pady=(0, 4))
        self.seq2_label.grid_remove()
        self.seq2_entry.grid_remove()
        self.seq2_hint.grid_remove()

        self.z_len_label = ttk.Label(helix, text="Z-DNA helix length:", style="Field.TLabel")
        self.z_len_var = tk.StringVar(value="")
        self.z_len_entry = ttk.Entry(helix, textvariable=self.z_len_var, width=22, state="disabled")
        self.z_hint = ttk.Label(
            helix,
            text="Z-DNA fiber has a fixed sequence; the Sequence field and DSSR parameter table are ignored.",
            style="Hint.TLabel",
            wraplength=520,
        )
        self.z_len_label.grid(row=2, column=0, sticky="w", padx=label_pad, **field)
        self.z_len_entry.grid(row=2, column=1, sticky="w", **field)
        self.z_hint.grid(row=3, column=1, columnspan=2, sticky="w", pady=(0, 4))
        self.z_len_label.grid_remove()
        self.z_len_entry.grid_remove()
        self.z_hint.grid_remove()

        ttk.Label(helix, text="Helix name (optional):", style="Field.TLabel").grid(
            row=4, column=0, sticky="w", padx=label_pad, **field
        )
        self.name_var = tk.StringVar(value="")
        self.name_entry = ttk.Entry(helix, textvariable=self.name_var)
        self.name_entry.grid(row=4, column=1, columnspan=2, sticky="we", **field)

        ttk.Label(helix, text="Output folder:", style="Field.TLabel").grid(
            row=5, column=0, sticky="w", padx=label_pad, **field
        )
        out_row = ttk.Frame(helix)
        out_row.grid(row=5, column=1, columnspan=2, sticky="we", **field)
        out_row.columnconfigure(0, weight=1)
        self.output_dir_var = tk.StringVar(value=_path_text(DEFAULT_OUTPUT_DIR))
        self.output_dir_entry = ttk.Entry(out_row, textvariable=self.output_dir_var)
        self.output_dir_entry.grid(row=0, column=0, sticky="we")
        ttk.Button(out_row, text="Browse", command=self.browse_output_dir).grid(row=0, column=1, padx=(8, 0))

        self.deleteH_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            helix,
            text="Delete hydrogens from generated PDB",
            variable=self.deleteH_var,
        ).grid(row=6, column=1, columnspan=2, sticky="w", pady=(6, 0))

        # DSSR helical parameters ------------------------------------------------
        self.param_frame = ui.Card(form, "Current DSSR helical parameters")
        self.param_frame.grid(row=1, column=0, columnspan=2, sticky="we", pady=(0, 12))
        self.param_frame.columnconfigure(0, weight=1)

        self.param_controls = ttk.Frame(self.param_frame)
        self.param_controls.grid(row=0, column=0, sticky="we", pady=(0, 10))
        self.param_controls.columnconfigure(1, weight=1)
        self.custom_btn = ttk.Button(
            self.param_controls,
            text="Customize DSSR parameters",
            command=self.open_param_dialog,
        )
        self.custom_btn.grid(row=0, column=0, sticky="w")
        self.param_status_var = tk.StringVar(value="")
        ttk.Label(
            self.param_controls,
            textvariable=self.param_status_var,
            style="Hint.TLabel",
            justify="left",
            wraplength=520,
        ).grid(row=0, column=1, sticky="w", padx=(14, 0))

        self.param_table = ttk.Treeview(
            self.param_frame,
            columns=("row_label", *PARAM_KEYS),
            show="headings",
            height=3,
            selectmode="none",
            style="Compact.Treeview",
        )
        for column, heading, width, anchor in (("row_label", "", 60, "w"),):
            self.param_table.heading(column, text=heading)
            self.param_table.column(column, width=width, minwidth=56, anchor=anchor, stretch=False)
        for key in PARAM_KEYS:
            self.param_table.heading(key, text=PARAM_LABELS.get(key, key))
            self.param_table.column(key, width=82, minwidth=60, anchor="center", stretch=True)
        self.param_table.grid(row=1, column=0, sticky="we")

        # P-DNA: the template's fixed helical parameters, read-only.
        self.pdna_param_table = ttk.Treeview(
            self.param_frame,
            columns=("row_label", *[f"p{idx}" for idx in range(len(PDNA_PARAM_COLUMNS))]),
            show="headings",
            height=2,
            selectmode="none",
            style="Compact.Treeview",
        )
        self.pdna_param_table.heading("row_label", text="")
        self.pdna_param_table.column("row_label", width=60, minwidth=56, anchor="w", stretch=False)
        for idx, (heading, _value) in enumerate(PDNA_PARAM_COLUMNS):
            self.pdna_param_table.heading(f"p{idx}", text=heading)
            self.pdna_param_table.column(f"p{idx}", width=96, minwidth=72, anchor="center", stretch=True)
        self.pdna_param_table.insert(
            "", "end", values=("Fixed", *[_format_number(value) for _heading, value in PDNA_PARAM_COLUMNS])
        )
        self.pdna_param_table.insert(
            "", "end", values=("Source", *["1PFI" for _column in PDNA_PARAM_COLUMNS])
        )
        self.pdna_param_table.grid(row=1, column=0, sticky="we")
        self.pdna_param_table.grid_remove()

        left_col = ttk.Frame(form)
        left_col.grid(row=2, column=0, sticky="nsew", padx=(0, 6))
        left_col.columnconfigure(0, weight=1)
        right_col = ttk.Frame(form)
        right_col.grid(row=2, column=1, sticky="nsew", padx=(6, 0))
        right_col.columnconfigure(0, weight=1)

        def wrap_to(widget: tk.Misc, label: tk.Misc, margin: int = 30) -> None:
            widget.bind("<Configure>", lambda event: label.configure(wraplength=max(160, event.width - margin)),
                        add="+")

        # phenix.geometry_minimization ----------------------------------------------
        self.min_frame = ui.Card(left_col, "phenix.geometry_minimization")
        self.min_frame.grid(row=0, column=0, sticky="we", pady=(0, 12))
        self.min_frame.columnconfigure(1, weight=1)
        self.minimize_var = tk.BooleanVar(value=DEFAULT_MINIMIZE_BY_TYPE["B-DNA"])
        self.min_check = ttk.Checkbutton(
            self.min_frame,
            text="Run phenix.geometry_minimization",
            variable=self.minimize_var,
            command=self._on_minimize_toggled,
        )
        self.min_check.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 2))
        ttk.Label(self.min_frame, text="Params file (.eff / .params):", style="Field.TLabel").grid(
            row=1, column=0, columnspan=3, sticky="w", pady=(6, 2)
        )
        self.params_var = tk.StringVar(value=_path_text(DEFAULT_PARAMS_FILE))
        self.params_entry = ttk.Entry(self.min_frame, textvariable=self.params_var)
        self.params_entry.grid(row=2, column=0, columnspan=2, sticky="we", pady=(0, 4))
        self.params_browse_btn = ttk.Button(self.min_frame, text="Browse", command=self.browse_params)
        self.params_browse_btn.grid(row=2, column=2, sticky="w", padx=(8, 0), pady=(0, 4))
        self.regularize_phosphates_var = tk.BooleanVar(value=DEFAULT_MINIMIZE_BY_TYPE["B-DNA"])
        ttk.Checkbutton(
            self.min_frame,
            text="Regularize phosphates",
            variable=self.regularize_phosphates_var,
            command=self._on_regularize_phosphates_toggled,
        ).grid(row=3, column=0, columnspan=3, sticky="w", pady=(8, 0))
        regularize_hint = ttk.Label(
            self.min_frame,
            text="Use C1' atoms to restore helical symmetry to all atoms selected in min_P_C5.params.",
            style="Hint.TLabel",
            wraplength=420,
            justify="left",
        )
        regularize_hint.grid(row=4, column=0, columnspan=3, sticky="w", padx=(22, 0))
        self.min_hint_var = tk.StringVar(value="")
        min_hint = ttk.Label(self.min_frame, textvariable=self.min_hint_var, style="Hint.TLabel", wraplength=420,
                             justify="left")
        min_hint.grid(row=5, column=0, columnspan=3, sticky="w", pady=(8, 0))
        wrap_to(self.min_frame, regularize_hint, 60)
        wrap_to(self.min_frame, min_hint, 36)

        # Mirror-image L-form chirality ------------------------------------------------
        self.chirality_frame = ui.Card(left_col, "Mirror-image L-form chirality (L-DNA)")
        self.chirality_frame.grid(row=1, column=0, sticky="we")
        self.chirality_frame.columnconfigure(3, weight=1)
        self.invrot_enabled_var = tk.BooleanVar(value=False)
        self.invrot_operation_var = tk.StringVar(value="oyz")
        self.invrot_check = ttk.Checkbutton(
            self.chirality_frame,
            text="Apply inv/rot after align-to-Z",
            variable=self.invrot_enabled_var,
            command=self._on_invrot_toggled,
        )
        self.invrot_check.grid(row=0, column=0, columnspan=4, sticky="w", **field)
        ttk.Label(self.chirality_frame, text="Operation:", style="Field.TLabel").grid(
            row=1, column=0, sticky="w", padx=(0, 8), **field
        )
        self.invrot_operation_combo = ttk.Combobox(
            self.chirality_frame,
            textvariable=self.invrot_operation_var,
            values=INV_ROT_OPERATIONS,
            width=8,
            state="disabled",
        )
        self.invrot_operation_combo.grid(row=1, column=1, sticky="w", **field)
        self.invrot_operation_combo.bind("<<ComboboxSelected>>", lambda *_args: self._refresh_info_text())
        ttk.Button(self.chirality_frame, text="Help", command=self.open_invrot_help, style="Ghost.TButton").grid(
            row=1, column=2, sticky="w", padx=(10, 0), **field
        )
        self.invrot_hint_var = tk.StringVar(value="")
        invrot_hint = ttk.Label(self.chirality_frame, textvariable=self.invrot_hint_var, style="Hint.TLabel",
                                wraplength=420, justify="left")
        invrot_hint.grid(row=2, column=0, columnspan=4, sticky="w", pady=(4, 0))
        wrap_to(self.chirality_frame, invrot_hint, 36)
        self._refresh_invrot_state()

        # Placement / orientation ---------------------------------------------------------
        place = ui.Card(right_col, "Placement / Orientation")
        place.grid(row=0, column=0, sticky="we", pady=(0, 12))
        for col in (1, 3, 5):
            place.grid_columnconfigure(col, weight=1, uniform="place")

        self.x_var = tk.StringVar(value="0")
        self.y_var = tk.StringVar(value="0")
        self.z_var = tk.StringVar(value="0")
        self.roll_var = tk.StringVar(value="0")
        self.phi_var = tk.StringVar(value="0")
        self.theta_var = tk.StringVar(value="0")
        self.delta_z_var = tk.StringVar(value="0")

        def place_field(row: int, col: int, text: str, var: tk.StringVar) -> None:
            ttk.Label(place, text=text, style="Field.TLabel").grid(
                row=row, column=col, sticky="w", padx=(0 if col == 0 else 14, 6), **field
            )
            ttk.Entry(place, textvariable=var, width=7).grid(row=row, column=col + 1, sticky="we", **field)

        place_field(0, 0, "delta_z (Å)", self.delta_z_var)
        ttk.Label(place, text="delta_z should be 0 for most cases.", style="Hint.TLabel").grid(
            row=0, column=2, columnspan=4, sticky="w", padx=(14, 0), **field
        )
        place_field(1, 0, "x (Å)", self.x_var)
        place_field(1, 2, "y (Å)", self.y_var)
        place_field(1, 4, "z (Å)", self.z_var)
        place_field(2, 0, "roll (°)", self.roll_var)
        place_field(2, 2, "phi (°)", self.phi_var)
        place_field(2, 4, "theta (°)", self.theta_var)
        ttk.Label(
            place,
            text="For GIDEON: roll = roll at GIDEON - 111.25",
            style="Hint.TLabel",
            justify="left",
            wraplength=520,
        ).grid(row=3, column=0, columnspan=6, sticky="w", pady=(6, 0))

        # Pipeline summary ------------------------------------------------------------------
        summary = tk.Frame(right_col, background=ui.C.SUBTLE, highlightthickness=1,
                           highlightbackground=ui.C.BORDER, highlightcolor=ui.C.BORDER)
        summary.grid(row=1, column=0, sticky="we")
        self.info_var = tk.StringVar(value="")
        info_label = ttk.Label(summary, textvariable=self.info_var, style="Subtle.Hint.TLabel", justify="left",
                               wraplength=700)
        info_label.pack(fill="x", padx=14, pady=10)
        summary.bind("<Configure>", lambda event: info_label.configure(wraplength=max(200, event.width - 30)))

        # --- Tools ------------------------------------------------------------------------
        side.columnconfigure(0, weight=1)
        side.rowconfigure(2, weight=1)

        builders = ui.ToolList(side, title="Builders")
        builders.grid(row=0, column=0, sticky="we")
        builders.add("B-Z builder", command=self.open_bz_builder_dialog,
                     description="Join alternating B/Z PDBs.", icon="tool:bz_builder")
        builders.add("B-P builder", command=self.open_bp_builder_dialog,
                     description="Join alternating B/P-DNA PDBs.", icon="tool:bp_builder")
        builders.add("Triplex converter", command=self.open_triplex_converter_dialog,
                     description="Add strand III; tune X-disp (e.g. -1.9 Å) and h-Twist (e.g. 30°) first.",
                     icon="tool:triplex")

        def add_tool_button(parent, label, command, icon, message, help_title=None):
            """One entry of the Other tools list, with its ? help."""
            return parent.add(
                label,
                command=command,
                description=ui.tool_description(icon),
                icon=f"tool:{icon}",
                help_command=lambda: self._open_tool_help(help_title or label, message),
            )

        # --- Other tools frame ---
        tools = ui.ToolList(side, title="Other tools")
        tools.grid(row=1, column=0, sticky="we", pady=(12, 0))
        add_tool_button(
            tools,
            "Measure angle around axis",
            self.open_axis_angle_tool,
            "measure_angle",
            "Measure around-axis angles between two atom or XYZ points and write a Chimera/ChimeraX .bild drawing.",
        )
        add_tool_button(
            tools,
            "Align helix to z",
            self.open_align_to_z_dialog,
            "align_z",
            "Use DSSR --more axis endpoints to rotate and translate any helix PDB so its axis is aligned with the coordinate-system +Z axis.",
        )
        add_tool_button(
            tools,
            "Get helical-axis info",
            self.open_helical_axis_info_dialog,
            "axis_info",
            "Use DSSR on two selected PDB chains to report axis start/end points, direction, angle to a reference vector, and optional axis .bild output.",
        )
        add_tool_button(
            tools,
            "Add phosphates",
            self.open_add_phosphates_dialog,
            "add_phosphates",
            "Report and add terminal 5' and 3' phosphates by chain, with an option to add the O3' atom in residue n-1 before a 5' phosphate on residue n.",
        )
        add_tool_button(
            tools,
            "Write XYZ axes BILD",
            self.open_xyz_bild_dialog,
            "xyz_bild",
            "Create a coordinate-axis .bild helper with configurable origin, arrow length, arrow width, and origin marker.",
            help_title="XYZ axes BILD",
        )
        add_tool_button(
            tools,
            "Combine_PDB",
            self.open_combine_pdb_dialog,
            "combine_pdb",
            "Combine two or more PDB files, optionally using only selected chains from each file. Kept source chains are assigned consecutive chain IDs A, B, C, ... in input-file and first-appearance order; LINK and chain-bearing REMARK records are updated automatically, and metadata for unselected chains is dropped.",
        )
        add_tool_button(
            tools,
            "Regularize phosphates",
            self.open_regularize_phosphates_dialog,
            "regularize_phosphates",
            "Restore the regular repeating geometry of the phosphate backbone after Phenix minimization. "
            "The tool measures each chain's helical repeat from its C1' atoms, uses nonterminal residues "
            "to determine the typical backbone geometry, and adjusts only P, OP1, OP2, O5', C5', and O3'. "
            "Terminal atoms are corrected afterward using the internal geometry. All other atoms remain unchanged.",
        )

        log_frame = ttk.Frame(side)
        log_frame.grid(row=2, column=0, sticky="nsew", pady=(14, 0))
        log_frame.grid_columnconfigure(0, weight=1)
        log_frame.grid_rowconfigure(1, weight=1)
        ttk.Label(log_frame, text="Log output", style="Section.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 6))
        log_view = ui.LogView(log_frame, wrap="none", height=10, horizontal=True)
        log_view.grid(row=1, column=0, sticky="nsew")
        self.log_text = log_view.text
        self.log_text.configure(state="disabled")
        self.dssr_status_var.trace_add("write", lambda *_args: self._refresh_dssr_chip())

        self.seq_var.trace_add("write", lambda *_args: self._update_sequence_length())
        self.seq2_var.trace_add("write", lambda *_args: self._update_sequence_length())
        self.z_len_var.trace_add("write", lambda *_args: self._update_sequence_length())
        self.output_dir_var.trace_add("write", lambda *_args: self._refresh_info_text())
        self.invrot_operation_var.trace_add("write", lambda *_args: self._refresh_invrot_state())

        self.on_type_changed()
        self._refresh_info_text()
        self._set_log(f"=== {APP_NAME} {__version__} startup ===\nApp folder: {APP_DIR}\nHelper folder: {LIB_DIR}\n")
        self.after(100, self.check_dssr_on_startup)

    def _get_optional_window_icon(self) -> Optional[tk.PhotoImage]:
        icon = getattr(self, "_window_icon", None)
        if icon is not None:
            return icon
        if getattr(self, "_window_icon_unavailable", False):
            return None
        if not DEFAULT_ICON_FILE.exists():
            self._window_icon_unavailable = True
            return None
        try:
            self._window_icon = tk.PhotoImage(file=str(DEFAULT_ICON_FILE))
            return self._window_icon
        except tk.TclError:
            self._window_icon_unavailable = True
            return None

    def _set_optional_window_icon(self, window: Optional[tk.Misc] = None) -> None:
        target = window or self
        icon = self._get_optional_window_icon()
        if icon is None:
            return
        try:
            target.iconphoto(target is self, icon)
        except tk.TclError:
            pass

    def _set_app_identity(self) -> None:
        try:
            self.tk.call("tk", "appname", APP_NAME)
        except tk.TclError:
            pass

    # ------------------------------------------------------------------
    # Log helpers
    # ------------------------------------------------------------------
    def _set_log(self, text: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.insert("1.0", text)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")
        self.update_idletasks()

    def _append_log(self, text: str) -> None:
        self.log_text.configure(state="normal")
        if self.log_text.index("end-1c") != "1.0":
            self.log_text.insert("end", "\n")
        self.log_text.insert("end", text)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")
        self.update_idletasks()

    def _open_tool_help(self, title: str, message: str) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title(title)
        win.geometry("+260+180")
        win.minsize(420, 140)
        footer = ui.Footer(win)
        footer.pack(side="bottom", fill="x")
        ui.hairline(win).pack(side="bottom", fill="x")
        body = ttk.Frame(win, padding=(20, 16, 20, 18))
        body.pack(side="top", fill="both", expand=True)
        body.grid_columnconfigure(0, weight=1)
        ttk.Label(body, text=title, style="Title.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 8))
        ttk.Label(body, text=message, wraplength=480, justify="left").grid(row=1, column=0, sticky="nw")
        footer.add_button("Close", win.destroy)
        win.lift()
        win.focus_set()

    def _make_help_button(self, parent: tk.Misc, title: str, message: str) -> tk.Misc:
        return ui.HelpButton(parent, command=lambda: self._open_tool_help(title, message))

    def _refresh_dssr_chip(self) -> None:
        status = self.dssr_status_var.get()
        head, _sep, tail = status.partition(" — ")
        if "NOT FOUND" in status:
            kind = "error"
        elif "FOUND" in status:
            kind = "ok"
        else:
            kind = "busy"
        self.dssr_chip.set(kind, head)
        self.dssr_path_label.configure(text=tail)

    # ------------------------------------------------------------------
    # Browse / state helpers
    # ------------------------------------------------------------------
    def browse_output_dir(self) -> None:
        path = filedialog.askdirectory(title="Choose output folder", initialdir=self.output_dir_var.get() or str(APP_DIR))
        if path:
            self.output_dir_var.set(path)

    def browse_params(self) -> None:
        path = filedialog.askopenfilename(
            title="Choose phenix.geometry_minimization params file",
            initialdir=str(DEFAULT_PARAMS_FILE.parent),
            filetypes=[("Params/eff", "*.params *.eff"), ("All files", "*.*")],
        )
        if path:
            self.params_var.set(path)

    def _on_minimize_toggled(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type in self.minimize_by_type:
            self.minimize_by_type[na_type] = bool(self.minimize_var.get())
        self._refresh_minimization_hint()

    def _on_regularize_phosphates_toggled(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type in self.regularize_phosphates_by_type:
            self.regularize_phosphates_by_type[na_type] = bool(self.regularize_phosphates_var.get())
        self._refresh_info_text()

    def _refresh_minimization_hint(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type == "B-DNA":
            self.min_hint_var.set("Default for B-DNA is to run phenix.geometry_minimization, matching the previous B-DNA workflow.")
        elif na_type in ("A-DNA", "A-RNA"):
            self.min_hint_var.set("Default for A-DNA/A-RNA is to skip minimization; check the box to run phenix.geometry_minimization.")
        elif na_type == "P-DNA":
            self.min_hint_var.set(
                "Default for P-DNA is to keep the 1PFI template geometry as built; check the box to run "
                "phenix.geometry_minimization."
            )
        else:
            self.min_hint_var.set("")

    def _on_invrot_toggled(self) -> None:
        self._refresh_invrot_state()
        if hasattr(self, "info_var"):
            self._refresh_info_text()

    def _refresh_invrot_state(self) -> None:
        if not hasattr(self, "invrot_operation_combo"):
            return
        enabled = bool(self.invrot_enabled_var.get())
        self.invrot_operation_combo.configure(state="readonly" if enabled else "disabled")

        try:
            mode, label, rotation_axes = parse_operation(self.invrot_operation_var.get())
            if mode == "o":
                plane = label.replace("o_", "")
                hint = f"Reflection across {plane} plane; implemented as inversion + 180-degree rotation around {rotation_axes}."
            else:
                axes = rotation_axes if rotation_axes else "none"
                hint = f"Point inversion mode; 180-degree rotation axes: {axes}."
        except Exception as exc:
            hint = f"Invalid inv/rot operation: {exc}"
        self.invrot_hint_var.set(hint if enabled else "Disabled. Enable to generate a mirror-image L-form before placement.")

    def open_invrot_help(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("L-form inv/rot help")
        win.geometry("780x640+220+140")
        win.minsize(680, 460)
        win.transient(self)

        footer = ui.Footer(win)
        footer.pack(side="bottom", fill="x")
        ui.hairline(win).pack(side="bottom", fill="x")
        body = ttk.Frame(win, padding=(20, 16, 20, 16))
        body.pack(side="top", fill="both", expand=True)
        body.grid_columnconfigure(0, weight=1)
        body.grid_rowconfigure(0, weight=1)
        text_view = ui.LogView(body, wrap="word", height=22)
        text_view.grid(row=0, column=0, sticky="nsew")
        text = text_view.text
        text.configure(font="AzbBody", background=ui.C.BG, padx=16, pady=12, spacing1=2, spacing3=2)
        text.tag_configure("heading", font="AzbSection", spacing1=4, spacing3=2)
        text.tag_configure("code", font="AzbMono", foreground=ui.C.TEXT_2)

        help_text = """L-form mirror step

This option runs after DSSR align-to-Z and before final placement/orientation. At that point the helix axis is standardized, so the mirror operation is predictable and the normal roll/phi/theta/x/y/z placement controls still work afterward.

Why inv+rot changes chirality

A proper rotation has determinant +1, so it can turn an object but cannot turn a right-handed model into its mirror image. Point inversion maps:

    (x, y, z) -> (-x, -y, -z)

This has determinant -1, so it changes handedness. If you then add any 180-degree rotation, that rotation has determinant +1, so the combined operation still has determinant -1 and still changes chirality. Geometrically, inversion plus a 180-degree rotation is equivalent to reflection across a coordinate plane.

i mode

    i     inversion only: (-x, -y, -z)
    ix    inversion + 180-degree rotation around x
    iy    inversion + 180-degree rotation around y
    iz    inversion + 180-degree rotation around z
    ixy, ixz, iyz, ixyz are also accepted.

o mode

    oxy   reflection across the xy plane, equivalent to i + Rz(180)
    oyz   reflection across the yz plane, equivalent to i + Rx(180)
    oxz   reflection across the xz plane, equivalent to i + Ry(180)

After align-to-Z, oyz and oxz keep the z coordinate sign, so the aligned +Z direction is preserved. oxy and plain i change the z sign, so they reverse the aligned helix direction before placement.

The GUI default is oyz because it changes chirality while keeping the +Z axis direction."""

        text.insert("1.0", help_text)
        # Same text; only the look: short unindented lines are headings, indented lines are code.
        for line_number, line in enumerate(help_text.splitlines(), start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if line.startswith("    "):
                text.tag_add("code", f"{line_number}.0", f"{line_number}.end")
            elif len(stripped) < 40 and not stripped.endswith("."):
                text.tag_add("heading", f"{line_number}.0", f"{line_number}.end")
        text.configure(state="disabled")

        footer.add_button("Close", win.destroy)

    def open_bp_builder_help(self) -> None:
        """Help window of the B-P builder: short sections with bullets instead of a wall of text."""
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("B-P builder help")
        win.geometry("760x660+240+120")
        win.minsize(600, 420)
        win.transient(self)

        footer = ui.Footer(win)
        footer.pack(side="bottom", fill="x")
        ui.hairline(win).pack(side="bottom", fill="x")
        body = ttk.Frame(win, padding=(20, 16, 20, 16))
        body.pack(side="top", fill="both", expand=True)
        body.grid_columnconfigure(0, weight=1)
        body.grid_rowconfigure(0, weight=1)
        text_view = ui.LogView(body, wrap="word", height=24)
        text_view.grid(row=0, column=0, sticky="nsew")
        text = text_view.text
        text.configure(font="AzbBody", background=ui.C.BG, padx=16, pady=12, spacing1=2, spacing3=2)
        text.tag_configure("title", font="AzbTitle", spacing3=6)
        text.tag_configure("heading", font="AzbSection", spacing1=12, spacing3=4)
        text.tag_configure("bullet", lmargin1=12, lmargin2=30, spacing1=1, spacing3=3)
        text.tag_configure("para", spacing3=4)

        sections = [
            ("What a junction does", [
                "Between every neighbouring pair of input helices the builder inserts a modeled B-P junction "
                "core and fits it to both helices on their sugar atoms. Nothing is relaxed or refined: a build "
                "is rigid fitting, trimming, base grafting and ligation, and takes a fraction of a second.",
                "The B-DNA end hands its two terminal base pairs to the core, which carries the same bases, so "
                "every B-DNA base pair is kept.",
                "The P-DNA end's T terminal levels (T = Transition) are rebuilt as the transition, in which the "
                "backbone moves from the B-DNA radius (P atoms 9.5 Å from the axis) into the P-DNA core "
                "(2.5 Å). They keep the P-DNA's bases.",
                "The output strands therefore carry exactly the input sequences, and chain A of the first input "
                "stays chain A.",
            ]),
            ("Transition: 2 or 3 nucleotides per strand", [
                "2 (default): the compact core. Its strand I swings back about the axis as it enters the P-DNA "
                "(one backbone step of −134°), so a junction adds only about 25° of twist.",
                "3: every backbone step on both strands turns forward, between the B-DNA twist (34°) and the "
                "P-DNA twist (132°); all sugars are C2′-endo and all bases anti. A junction adds about "
                "267° of twist.",
                "Both are models: there is no experimental structure of a B-P junction.",
            ]),
            ("Choosing lengths", [
                "The transition nucleotides count as P-DNA: a P-DNA of N levels keeps N − T P-form levels "
                "with one junction and N − 2T between two junctions. For k P-form levels enter k + T levels "
                "(or k + 2T).",
                "Minimum P-DNA: T + 2 levels for one junction (4 with Transition 2, 5 with Transition 3) and "
                "7 levels between two junctions with either transition.",
                "Minimum B-DNA: 3 bp for one junction, 5 bp with junctions at both ends. No B-DNA base pair is "
                "lost.",
                "Twist of a B-P-B construct from its first to its last base pair: about 34.3° × "
                "(x + z − 2) + 131.8° × (y − 2T − 1) + 2 × 25° (Transition 2) "
                "or 2 × 267° (Transition 3), for x bp, y levels and z bp.",
            ]),
            ("Inputs and output", [
                "Use bnp_na-built helices. A default B-DNA (minimized and regularized) fits the core exactly, "
                "and every bnp_na P-DNA matches it exactly; other B-DNA models are fitted on their sugar atoms "
                "and the RMSD is reported in the log.",
                "Helix types: leave blank to detect B or P from the phosphate radius about the screw axis, or "
                "type one letter per file in construct order, e.g. BPB.",
                "Axis correction and the B-DNA axis source work as in the B-Z builder; P-DNA axes always come "
                "from the helix's own screw symmetry.",
                "Output: the ligated PDB with two continuous strands, and a raw aligned PDB (every piece in "
                "place, before trimming) with _raw before the extension.",
            ]),
        ]
        text.insert("end", "B-P structure builder\n", "title")
        text.insert("end", "Joins alternating B-DNA and P-DNA PDB files into one construct through modeled "
                           "junction cores.\n", "para")
        for heading, items in sections:
            text.insert("end", heading + "\n", "heading")
            for item in items:
                text.insert("end", "•  " + item + "\n", "bullet")
        text.configure(state="disabled")

        footer.add_button("Close", win.destroy)

    # ------------------------------------------------------------------
    # Shared look of the tool dialogs
    # ------------------------------------------------------------------
    def _tool_dialog_frame(
        self,
        win: tk.Toplevel,
        title: str,
        tool: Optional[str] = None,
        padding: Tuple[int, int, int, int] = (18, 14, 18, 14),
    ) -> Tuple[ttk.Frame, "ui.Footer"]:
        """Header and hairline on top, a footer for the buttons at the bottom,
        and the padded body in between. Returns (body, footer).

        App.__init__ has already themed the interpreter, and the theme's <Map>
        hook gives every new Toplevel its white background and light title bar
        once it is mapped. (ui.style_window() called before mapping cannot set
        the macOS appearance, and would stop that hook from doing it.)
        """
        header = ui.ToolHeader(win, title, tool=tool, module=None if tool else "bnp_na")
        header.pack(side="top", fill="x")
        ui.hairline(win).pack(side="top", fill="x")
        footer = ui.Footer(win)
        footer.pack(side="bottom", fill="x")
        ui.hairline(win).pack(side="bottom", fill="x")
        body = ttk.Frame(win, padding=padding)
        body.pack(side="top", fill="both", expand=True)
        footer.lift()  # stacking order is Tab order: the form's fields first, then the buttons
        return body, footer

    @staticmethod
    def _wrap_to_width(container: tk.Misc, label: tk.Misc, margin: int = 40, minimum: int = 240) -> None:
        """Re-wrap a hint label to the width of its container whenever that changes."""
        container.bind(
            "<Configure>",
            lambda event: label.configure(wraplength=max(minimum, event.width - margin)),
            add="+",
        )

    def open_axis_angle_tool(self) -> None:
        try:
            launch_axis_angle_gui(parent=self, log_callback=self._append_log)
        except Exception as exc:
            messagebox.showerror("Helical-axis angle tool", str(exc), parent=self)

    def open_align_to_z_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Align helix to z")
        win.geometry("840x560+240+160")
        win.minsize(700, 480)
        body, footer = self._tool_dialog_frame(win, "Align helix to z", tool="align_z")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(4, weight=1)

        def current_output_dir() -> Path:
            text = self.output_dir_var.get().strip()
            if not text:
                return DEFAULT_OUTPUT_DIR
            path = Path(text).expanduser()
            if not path.is_absolute():
                path = path.resolve()
            return path

        input_var = tk.StringVar(value="")
        output_var = tk.StringVar(value="")
        dssr_out_var = tk.StringVar(value="")
        last_default_output = {"path": ""}
        last_default_dssr = {"path": ""}

        def default_paths() -> Tuple[Path, Path]:
            pdb_text = input_var.get().strip()
            stem = Path(pdb_text).expanduser().stem if pdb_text else "helix"
            if not stem:
                stem = "helix"
            output_dir = current_output_dir()
            tmp_dir = output_dir / "tmp_file"
            return output_dir / f"{stem}_aligned2Z.pdb", tmp_dir / f"{stem}_align2z_dssr_more.out"

        def refresh_default_paths(force: bool = False) -> None:
            default_output, default_dssr = default_paths()
            default_output_text = str(default_output)
            default_dssr_text = str(default_dssr)
            if force or not output_var.get().strip() or output_var.get() == last_default_output["path"]:
                output_var.set(default_output_text)
            if force or not dssr_out_var.get().strip() or dssr_out_var.get() == last_default_dssr["path"]:
                dssr_out_var.set(default_dssr_text)
            last_default_output["path"] = default_output_text
            last_default_dssr["path"] = default_dssr_text

        def browse_input() -> None:
            path = filedialog.askopenfilename(
                title="Choose helix PDB",
                initialdir=str(current_output_dir()),
                filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                input_var.set(path)
                refresh_default_paths(force=True)

        def browse_output() -> None:
            current = Path(output_var.get().strip() or str(default_paths()[0])).expanduser()
            path = filedialog.asksaveasfilename(
                title="Save aligned PDB",
                initialdir=str(current.parent if current.parent else current_output_dir()),
                initialfile=current.name or "helix_aligned2Z.pdb",
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                output_var.set(path)

        def browse_dssr_output() -> None:
            current = Path(dssr_out_var.get().strip() or str(default_paths()[1])).expanduser()
            path = filedialog.asksaveasfilename(
                title="Save DSSR --more output",
                initialdir=str(current.parent if current.parent else current_output_dir()),
                initialfile=current.name or "align2z_dssr_more.out",
                defaultextension=".out",
                filetypes=[("DSSR output", "*.out"), ("Text files", "*.txt"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                dssr_out_var.set(path)

        def set_result_text(text: str) -> None:
            result_text.configure(state="normal")
            result_text.delete("1.0", "end")
            result_text.insert("1.0", text)
            result_text.see("1.0")
            result_text.configure(state="disabled")

        def run_alignment() -> None:
            try:
                input_text = input_var.get().strip()
                if not input_text:
                    raise ValueError("Please choose an input PDB file.")
                input_path = Path(input_text).expanduser()
                if not input_path.is_absolute():
                    input_path = input_path.resolve()
                if not input_path.exists():
                    raise ValueError(f"Input PDB not found: {input_path}")

                output_text = output_var.get().strip()
                if not output_text:
                    raise ValueError("Please choose an aligned output PDB file.")
                output_path = Path(output_text).expanduser()
                if not output_path.is_absolute():
                    output_path = output_path.resolve()
                output_path.parent.mkdir(parents=True, exist_ok=True)

                dssr_text = dssr_out_var.get().strip()
                if not dssr_text:
                    raise ValueError("Please choose a DSSR .out file.")
                dssr_out_path = Path(dssr_text).expanduser()
                if not dssr_out_path.is_absolute():
                    dssr_out_path = dssr_out_path.resolve()
                if dssr_out_path.suffix.lower() != ".out":
                    raise ValueError("The DSSR output filename should end with .out.")
                dssr_out_path.parent.mkdir(parents=True, exist_ok=True)

                _output_dir, tmp_dir = _ensure_output_dirs(self.output_dir_var.get())
                run_btn.configure(state="disabled")
                win.update_idletasks()
                result = align_pdb_to_Z(
                    str(input_path),
                    out_pdb=str(output_path),
                    dssr_out=str(dssr_out_path),
                    cwd=tmp_dir,
                )
            except Exception as exc:
                messagebox.showerror("Align helix to z", str(exc), parent=win)
                return
            finally:
                try:
                    run_btn.configure(state="normal")
                except Exception:
                    pass

            cli_args = [
                "python3",
                "bnp_na_lib/align2z.py",
                str(input_path),
                "-o",
                str(output_path),
                "--dssr-output",
                str(dssr_out_path),
                "--cwd",
                str(tmp_dir),
            ]
            report = format_alignment_report(result)
            log = "\n".join(
                [
                    "=== Align helix to z tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    report,
                ]
            )
            set_result_text(report)
            self._append_log(log)
            messagebox.showinfo("Align helix to z", f"Aligned PDB:\n{output_path}", parent=win)

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Input PDB:", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=input_var).grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_input).grid(row=0, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Aligned output PDB:", style="Field.TLabel").grid(
            row=1, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=output_var).grid(row=1, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_output).grid(row=1, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="DSSR .out:", style="Field.TLabel").grid(row=2, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=dssr_out_var).grid(row=2, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_dssr_output).grid(
            row=2, column=2, sticky="we", padx=(8, 0), **field
        )

        hint = ttk.Label(
            body,
            text="Uses DSSR --more point-one/point-two endpoints, then translates point-one to the origin and rotates the axis to +Z.",
            style="Hint.TLabel",
            wraplength=620,
            justify="left",
        )
        hint.grid(row=3, column=0, columnspan=3, sticky="w", pady=(6, 10))
        self._wrap_to_width(body, hint)

        result_view = ui.LogView(body, wrap="none", height=9, horizontal=True)
        result_view.grid(row=4, column=0, columnspan=3, sticky="nsew")
        result_text = result_view.text
        result_text.configure(state="disabled")

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Align to z", run_alignment, primary=True)

        refresh_default_paths(force=True)
        input_var.trace_add("write", lambda *_args: refresh_default_paths())

    def open_helical_axis_info_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Get helical-axis info")
        win.geometry("860x740+220+120")
        win.minsize(740, 620)
        body, footer = self._tool_dialog_frame(win, "Get helical-axis info", tool="axis_info")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(7, weight=1)

        def current_output_dir() -> Path:
            text = self.output_dir_var.get().strip()
            if not text:
                return DEFAULT_OUTPUT_DIR
            path = Path(text).expanduser()
            if not path.is_absolute():
                path = path.resolve()
            return path

        pdb_var = tk.StringVar(value="")
        chains_var = tk.StringVar(value="C D")
        vector_var = tk.StringVar(value="0 0 1")
        helix_bp_var = tk.StringVar(value="")
        write_bild_var = tk.BooleanVar(value=True)
        draw_ref_var = tk.BooleanVar(value=True)
        ref_length_var = tk.StringVar(value="")
        bild_var = tk.StringVar(value="")
        last_default_bild = {"path": ""}

        def default_bild_path() -> Path:
            pdb_text = pdb_var.get().strip()
            stem = Path(pdb_text).expanduser().stem if pdb_text else "helical_axis"
            if not stem:
                stem = "helical_axis"
            try:
                chains = parse_chain_ids(chains_var.get())
                chain_text = "".join(ch if ch.isalnum() else "_" for ch in chains)
            except Exception:
                chain_text = "chains"
            return current_output_dir() / f"{stem}_chains_{chain_text}_axis.bild"

        def refresh_default_bild(force: bool = False) -> None:
            new_default = str(default_bild_path())
            if force or not bild_var.get().strip() or bild_var.get() == last_default_bild["path"]:
                bild_var.set(new_default)
            last_default_bild["path"] = new_default

        def browse_input() -> None:
            path = filedialog.askopenfilename(
                title="Choose PDB file",
                initialdir=str(current_output_dir()),
                filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                pdb_var.set(path)
                refresh_default_bild(force=True)

        def browse_bild() -> None:
            path = filedialog.asksaveasfilename(
                title="Save helical-axis BILD",
                initialfile=Path(default_bild_path()).name,
                defaultextension=".bild",
                filetypes=[("BILD files", "*.bild"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                bild_var.set(path)

        def parse_reference_vector() -> Optional[Tuple[float, float, float]]:
            text = vector_var.get().strip()
            if not text:
                return None
            parts = text.replace(",", " ").split()
            if len(parts) != 3:
                raise ValueError("Reference vector must have three numbers, for example: 0 0 1")
            return (
                _parse_float_expression(parts[0], "Reference vector x"),
                _parse_float_expression(parts[1], "Reference vector y"),
                _parse_float_expression(parts[2], "Reference vector z"),
            )

        def parse_reference_vector_length() -> Optional[float]:
            text = ref_length_var.get().strip()
            if not text:
                return None
            length = _parse_float_expression(text, "Reference-vector drawing length")
            if length <= 0:
                raise ValueError("Reference-vector drawing length must be positive.")
            return length

        def parse_helix_bp_count() -> Optional[int]:
            text = helix_bp_var.get().strip()
            if not text:
                return None
            bp_count = _parse_int_expression(text, "Helix length in bp")
            if bp_count <= 1:
                raise ValueError("Helix length in bp must be greater than 1.")
            return bp_count

        def set_result_text(text: str) -> None:
            result_text.configure(state="normal")
            result_text.delete("1.0", "end")
            result_text.insert("1.0", text)
            result_text.see("1.0")
            result_text.configure(state="disabled")

        def refresh_bild_state() -> None:
            state = "normal" if write_bild_var.get() else "disabled"
            bild_entry.configure(state=state)
            bild_browse_btn.configure(state=state)
            draw_ref_check.configure(state=state)
            ref_length_entry.configure(state="normal" if write_bild_var.get() and draw_ref_var.get() else "disabled")
            if write_bild_var.get():
                refresh_default_bild()

        def run_info() -> None:
            try:
                pdb_text = pdb_var.get().strip()
                if not pdb_text:
                    raise ValueError("Please choose an input PDB file.")
                pdb_path = Path(pdb_text).expanduser()
                if not pdb_path.is_absolute():
                    pdb_path = pdb_path.resolve()
                chains = parse_chain_ids(chains_var.get())
                reference_vector = parse_reference_vector()
                reference_vector_length = parse_reference_vector_length()
                helix_bp_count = parse_helix_bp_count()
                output_dir, tmp_dir = _ensure_output_dirs(self.output_dir_var.get())
                if write_bild_var.get():
                    bild_text = bild_var.get().strip()
                    if not bild_text:
                        raise ValueError("Please choose an output .bild file, or turn off BILD output.")
                    bild_output: Optional[Path] = Path(bild_text).expanduser()
                    if not bild_output.is_absolute():
                        bild_output = bild_output.resolve()
                else:
                    bild_output = None

                run_btn.configure(state="disabled")
                win.update_idletasks()
                info = get_helical_axis_info(
                    pdb_path,
                    chains,
                    reference_vector=reference_vector,
                    workdir=tmp_dir,
                    bild_output=bild_output,
                    draw_reference_vector=bool(draw_ref_var.get()),
                    reference_vector_length=reference_vector_length,
                    helix_bp_count=helix_bp_count,
                )
            except Exception as exc:
                messagebox.showerror("Get helical-axis info", str(exc), parent=win)
                return
            finally:
                try:
                    run_btn.configure(state="normal")
                except Exception:
                    pass

            cli_args = [
                "python3",
                "bnp_na_lib/helical_axis_info.py",
                "-i",
                str(pdb_path),
                "--chains",
                f"{chains[0]} {chains[1]}",
                "--workdir",
                str(output_dir / "tmp_file"),
            ]
            if reference_vector is not None:
                cli_args.extend(["--vector", " ".join(_format_number(value) for value in reference_vector)])
            if helix_bp_count is not None:
                cli_args.extend(["--helix-bp", str(helix_bp_count)])
            if bild_output is not None:
                cli_args.extend(["--bild", str(bild_output)])
                if not draw_ref_var.get():
                    cli_args.append("--no-reference-vector-bild")
                if reference_vector_length is not None:
                    cli_args.extend(["--reference-vector-length", _format_number(reference_vector_length)])

            report = format_axis_info_report(info)
            log = "\n".join(
                [
                    "=== Get helical-axis info tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    report,
                ]
            )
            set_result_text(report)
            self._append_log(log)
            if info.bild_output is not None:
                messagebox.showinfo("Get helical-axis info", f"Wrote BILD:\n{info.bild_output}", parent=win)
            else:
                messagebox.showinfo("Get helical-axis info", "Helical-axis info has been written to the log.", parent=win)

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Input PDB:", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=pdb_var).grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_input).grid(row=0, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Chain IDs:", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=chains_var, width=14).grid(row=1, column=1, sticky="w", **field)

        ttk.Label(body, text="Reference vector:", style="Field.TLabel").grid(
            row=2, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=vector_var, width=18).grid(row=2, column=1, sticky="w", **field)

        ttk.Label(body, text="Helix length (bp):", style="Field.TLabel").grid(
            row=3, column=0, sticky="w", padx=label_pad, **field
        )
        helix_bp_frame = ttk.Frame(body)
        helix_bp_frame.grid(row=3, column=1, columnspan=2, sticky="w", **field)
        ttk.Entry(helix_bp_frame, textvariable=helix_bp_var, width=12).pack(side="left")
        ttk.Label(helix_bp_frame, text="Optional, for full-length estimate", style="Hint.TLabel").pack(
            side="left", padx=(12, 0)
        )

        ttk.Checkbutton(body, text="Write axis .bild", variable=write_bild_var, command=refresh_bild_state).grid(
            row=4, column=0, sticky="w", padx=label_pad, pady=(12, 4)
        )
        bild_entry = ttk.Entry(body, textvariable=bild_var)
        bild_entry.grid(row=4, column=1, sticky="we", pady=(12, 4))
        bild_browse_btn = ttk.Button(body, text="Browse", command=browse_bild)
        bild_browse_btn.grid(row=4, column=2, sticky="we", padx=(8, 0), pady=(12, 4))

        draw_ref_check = ttk.Checkbutton(
            body,
            text="Draw reference vector",
            variable=draw_ref_var,
            command=refresh_bild_state,
        )
        draw_ref_check.grid(row=5, column=0, sticky="w", padx=label_pad, **field)
        ref_length_frame = ttk.Frame(body)
        ref_length_frame.grid(row=5, column=1, columnspan=2, sticky="w", **field)
        ref_length_entry = ttk.Entry(ref_length_frame, textvariable=ref_length_var, width=14)
        ref_length_entry.pack(side="left")
        ttk.Label(ref_length_frame, text="Length (Å), blank = axis length", style="Hint.TLabel").pack(
            side="left", padx=(12, 0)
        )

        hint = ttk.Label(
            body,
            text="Selected-chain PDB and DSSR .out files are written to the selected output folder's tmp_file folder.",
            style="Hint.TLabel",
            wraplength=620,
            justify="left",
        )
        hint.grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 10))
        self._wrap_to_width(body, hint)

        result_view = ui.LogView(body, wrap="none", height=14, horizontal=True)
        result_view.grid(row=7, column=0, columnspan=3, sticky="nsew")
        result_text = result_view.text
        result_text.configure(state="disabled")

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Get info", run_info, primary=True)

        refresh_default_bild(force=True)
        pdb_var.trace_add("write", lambda *_args: refresh_default_bild(force=True))
        chains_var.trace_add("write", lambda *_args: refresh_default_bild(force=True))
        refresh_bild_state()

    def open_combine_pdb_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Combine_PDB")
        win.geometry("1000x820+220+90")
        win.minsize(860, 680)
        body, footer = self._tool_dialog_frame(win, "Combine_PDB", tool="combine_pdb")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(2, weight=1)
        body.grid_rowconfigure(4, weight=1)

        def current_output_dir() -> Path:
            text = self.output_dir_var.get().strip()
            if not text:
                return DEFAULT_OUTPUT_DIR
            path = Path(text).expanduser()
            if not path.is_absolute():
                path = path.resolve()
            return path

        file_count_var = tk.StringVar(value="2")
        input_vars = [tk.StringVar(value="") for _index in range(26)]
        chain_vars = [tk.StringVar(value="") for _index in range(26)]
        chain_hint_vars = [tk.StringVar(value="blank = all chains") for _index in range(26)]
        chain_scan_cache: Dict[Tuple[str, int], str] = {}
        output_var = tk.StringVar(value=str(default_combine_pdb_output_path(current_output_dir())))

        def set_dialog_log(text: str) -> None:
            result_text.configure(state="normal")
            result_text.delete("1.0", "end")
            result_text.insert("1.0", text)
            result_text.see("1.0")
            result_text.configure(state="disabled")
            win.update_idletasks()

        def refresh_chain_hint(index: int) -> None:
            """Show the chain IDs found in the chosen file next to its chain field."""
            text = input_vars[index].get().strip()
            path = Path(text).expanduser() if text else None
            try:
                key = (str(path), path.stat().st_mtime_ns) if path else None
            except OSError:
                key = None
            if key is None:
                chain_hint_vars[index].set("blank = all chains")
                return
            found = chain_scan_cache.get(key)
            if found is None:
                try:
                    found = ", ".join(format_combine_chain_id(chain) for chain in list_pdb_chains(path))
                except Exception:
                    found = ""
                chain_scan_cache[key] = found
            chain_hint_vars[index].set(f"in file: {found}" if found else "blank = all chains")

        def browse_input(index: int) -> None:
            current_text = input_vars[index].get().strip()
            initial_dir = current_output_dir()
            if current_text:
                current_path = Path(current_text).expanduser()
                if str(current_path.parent) not in ("", "."):
                    initial_dir = current_path.parent
            selected = filedialog.askopenfilename(
                title=f"Choose input PDB {index + 1}",
                initialdir=str(initial_dir),
                filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                input_vars[index].set(selected)

        def render_input_fields(*_args: object) -> None:
            for child in input_frame.winfo_children():
                child.destroy()
            count = int(file_count_var.get())
            input_frame.grid_columnconfigure(1, weight=1)
            for index in range(count):
                ttk.Label(input_frame, text=f"Input PDB {index + 1}:", style="Field.TLabel").grid(
                    row=index, column=0, sticky="w", padx=(0, 12), pady=4
                )
                ttk.Entry(input_frame, textvariable=input_vars[index]).grid(
                    row=index, column=1, sticky="we", pady=4
                )
                ttk.Button(
                    input_frame,
                    text="Browse",
                    command=lambda input_index=index: browse_input(input_index),
                ).grid(row=index, column=2, sticky="w", padx=(8, 18), pady=4)
                ttk.Label(input_frame, text="Chains:", style="Field.TLabel").grid(
                    row=index, column=3, sticky="w", padx=(0, 8), pady=4
                )
                ttk.Entry(input_frame, textvariable=chain_vars[index], width=12).grid(
                    row=index, column=4, sticky="w", pady=4
                )
                ttk.Label(
                    input_frame,
                    textvariable=chain_hint_vars[index],
                    style="Hint.TLabel",
                    width=22,
                ).grid(row=index, column=5, sticky="w", padx=(12, 6), pady=4)
                refresh_chain_hint(index)
            input_frame.update_idletasks()
            input_canvas.configure(scrollregion=input_canvas.bbox("all"))
            input_canvas.yview_moveto(0.0)

        def browse_output() -> None:
            current = Path(
                output_var.get().strip() or str(default_combine_pdb_output_path(current_output_dir()))
            ).expanduser()
            selected = filedialog.asksaveasfilename(
                title="Save combined PDB",
                initialdir=str(current.parent if str(current.parent) not in ("", ".") else current_output_dir()),
                initialfile=current.name or "combine_PDB_out.pdb",
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                output_var.set(selected)

        def run_combine() -> None:
            try:
                count = int(file_count_var.get())
                input_paths: List[Path] = []
                chain_selections: List[Optional[List[str]]] = []
                for index in range(count):
                    value = input_vars[index].get().strip()
                    if not value:
                        raise ValueError(f"Please choose input PDB {index + 1}.")
                    input_paths.append(Path(value).expanduser())
                    chain_selections.append(parse_combine_chain_selection(chain_vars[index].get()))

                output_text = output_var.get().strip()
                if not output_text:
                    raise ValueError("Please choose an output PDB file.")
                output_path = Path(output_text).expanduser()

                run_btn.configure(state="disabled")
                win.update_idletasks()
                result = combine_pdb_files(input_paths, output_path, chain_selections)
            except Exception as exc:
                messagebox.showerror("Combine_PDB", str(exc), parent=win)
                return
            finally:
                try:
                    run_btn.configure(state="normal")
                except Exception:
                    pass

            cli_args = [
                "python3",
                "bnp_na_lib/combine_pdb.py",
                *[str(path) for path in result.input_pdbs],
            ]
            if any(selection is not None for selection in chain_selections):
                for selection in chain_selections:
                    cli_args.extend(
                        [
                            "--chains",
                            " ".join(format_combine_chain_token(chain) for chain in selection)
                            if selection
                            else "all",
                        ]
                    )
            cli_args.extend(["-o", str(result.output_pdb)])
            log = "\n".join(
                [
                    "=== Combine_PDB tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    result.log_text,
                ]
            )
            set_dialog_log(result.log_text)
            self._append_log(log)
            messagebox.showinfo(
                "Combine_PDB",
                f"Combined {len(result.input_pdbs)} PDB files into {result.chain_count} chains.\n\n"
                f"Wrote:\n{result.output_pdb}",
                parent=win,
            )

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Number of input PDB files:", style="Field.TLabel").grid(
            row=0, column=0, sticky="w", padx=label_pad, **field
        )
        count_combo = ttk.Combobox(
            body,
            textvariable=file_count_var,
            values=tuple(str(value) for value in range(2, 27)),
            state="readonly",
            width=5,
        )
        count_combo.grid(row=0, column=1, sticky="w", **field)
        count_combo.bind("<<ComboboxSelected>>", render_input_fields)

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=output_var).grid(row=1, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_output).grid(row=1, column=2, sticky="we", padx=(8, 0), **field)

        input_box = ui.Card(body, "Input PDB files (chain assignment follows this order)")
        input_box.grid(row=2, column=0, columnspan=3, sticky="nsew", pady=(12, 8))
        input_box.grid_columnconfigure(0, weight=1)
        input_box.grid_rowconfigure(0, weight=1)
        # A scrolling area (the shared ScrollFrame: its own scroll region, a
        # scrollbar only when the rows overflow, and wheel scrolling).
        input_area = ui.ScrollFrame(input_box)
        input_area.grid(row=0, column=0, sticky="nsew")
        input_canvas = input_area.canvas
        input_canvas.configure(height=190)
        input_frame = input_area.body

        combine_hint = ttk.Label(
            body,
            text="Leave a Chains field blank (or type all) to use every chain of that file; otherwise list the "
            "chain IDs to keep, such as A B or A,B, using _ for a blank chain ID. Selected chains are detected "
            "in first-appearance order within each file and renamed A-Z across all inputs, so the order you type "
            "them in does not reorder them. Coordinates, LINK endpoints, and chain-bearing REMARK/HET metadata "
            "belonging to unselected chains are dropped. The same PDB may be selected more than once, with a "
            "different chain selection each time; every occurrence receives fresh chain IDs. The PDB format "
            "limits the combined output to 26 chains.",
            style="Hint.TLabel",
            wraplength=920,
            justify="left",
        )
        combine_hint.grid(row=3, column=0, columnspan=3, sticky="w", pady=(0, 10))
        self._wrap_to_width(body, combine_hint)

        result_view = ui.LogView(body, wrap="none", height=8, horizontal=True)
        result_view.grid(row=4, column=0, columnspan=3, sticky="nsew")
        result_text = result_view.text
        result_text.configure(state="disabled")

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Combine PDB files", run_combine, primary=True)

        for index in range(len(input_vars)):
            input_vars[index].trace_add(
                "write", lambda *_args, row=index: refresh_chain_hint(row)
            )
        render_input_fields()
        set_dialog_log(
            "Choose two or more PDB files. Input order determines the new chain-ID order, and each file's "
            "Chains field selects which of its chains are combined (blank = all)."
        )

    def open_regularize_phosphates_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Regularize phosphates")
        win.geometry("880x620+220+140")
        win.minsize(740, 520)
        body, footer = self._tool_dialog_frame(win, "Regularize phosphates", tool="regularize_phosphates")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(3, weight=1)

        def current_output_dir() -> Path:
            text = self.output_dir_var.get().strip()
            if not text:
                return DEFAULT_OUTPUT_DIR
            path = Path(text).expanduser()
            return path if path.is_absolute() else path.resolve()

        input_var = tk.StringVar(value="")
        output_var = tk.StringVar(value="")
        last_default_output = {"path": ""}

        def default_output_path() -> Path:
            input_text = input_var.get().strip()
            if input_text:
                return default_regularized_output_path(Path(input_text).expanduser())
            return current_output_dir() / "regularized_phosphates_out.pdb"

        def refresh_default_output(force: bool = False) -> None:
            default_text = str(default_output_path())
            current = output_var.get().strip()
            if force or not current or current == last_default_output["path"]:
                output_var.set(default_text)
            last_default_output["path"] = default_text

        def browse_input() -> None:
            selected = filedialog.askopenfilename(
                title="Choose input PDB",
                initialdir=str(current_output_dir()),
                filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                input_var.set(selected)
                refresh_default_output(force=True)

        def browse_output() -> None:
            refresh_default_output()
            current = Path(output_var.get().strip() or str(default_output_path())).expanduser()
            selected = filedialog.asksaveasfilename(
                title="Save phosphate-regularized PDB",
                initialdir=str(current.parent),
                initialfile=current.name,
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                output_var.set(selected)

        def set_dialog_log(text: str) -> None:
            result_text.configure(state="normal")
            result_text.delete("1.0", "end")
            result_text.insert("1.0", text)
            result_text.see("1.0")
            result_text.configure(state="disabled")
            win.update_idletasks()

        def run_regularization() -> None:
            try:
                input_text = input_var.get().strip()
                if not input_text:
                    raise ValueError("Please choose an input PDB file.")
                input_path = Path(input_text).expanduser()
                if not input_path.is_absolute():
                    input_path = input_path.resolve()
                refresh_default_output()
                output_text = output_var.get().strip()
                if not output_text:
                    raise ValueError("Please choose an output PDB file.")
                output_path = Path(output_text).expanduser()
                if not output_path.is_absolute():
                    output_path = output_path.resolve()

                run_btn.configure(state="disabled")
                win.update_idletasks()
                result = regularize_phosphates(input_path, output_path)
            except Exception as exc:
                messagebox.showerror("Regularize phosphates", str(exc), parent=win)
                return
            finally:
                try:
                    run_btn.configure(state="normal")
                except Exception:
                    pass

            cli_args = [
                "python3",
                "bnp_na_lib/regularize_phosphates.py",
                str(input_path),
                "-o",
                str(output_path),
            ]
            log = "\n".join(
                [
                    "=== Regularize phosphates tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    result.log_text,
                ]
            )
            set_dialog_log(result.log_text)
            self._append_log(log)
            messagebox.showinfo("Regularize phosphates", f"Wrote:\n{result.output_pdb}", parent=win)

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Input PDB:", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=input_var).grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_input).grid(row=0, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=output_var).grid(row=1, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_output).grid(row=1, column=2, sticky="we", padx=(8, 0), **field)

        hint = ttk.Label(
            body,
            text=(
                "Fits one helical step per chain from consecutive C1' atoms. Internal P/OP1/OP2/O5'/C5'/O3' "
                "positions define the consensus; affected terminal atoms are excluded initially and regularized afterward."
            ),
            style="Hint.TLabel",
            wraplength=790,
            justify="left",
        )
        hint.grid(row=2, column=0, columnspan=3, sticky="w", pady=(6, 10))
        self._wrap_to_width(body, hint)

        result_view = ui.LogView(body, wrap="none", height=14, horizontal=True)
        result_view.grid(row=3, column=0, columnspan=3, sticky="nsew")
        result_text = result_view.text
        result_text.configure(state="disabled")

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Regularize phosphates", run_regularization, primary=True)

        refresh_default_output(force=True)
        set_dialog_log("Choose an input PDB. The output defaults beside the source file.")
        input_var.trace_add("write", lambda *_args: refresh_default_output())

    def open_add_phosphates_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Add phosphates")
        win.geometry("920x740+220+110")
        win.minsize(780, 620)
        body, footer = self._tool_dialog_frame(win, "Add phosphates", tool="add_phosphates")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(4, weight=1)

        def current_output_dir() -> Path:
            text = self.output_dir_var.get().strip()
            if not text:
                return DEFAULT_OUTPUT_DIR
            path = Path(text).expanduser()
            if not path.is_absolute():
                path = path.resolve()
            return path

        input_var = tk.StringVar(value="")
        output_var = tk.StringVar(value="")
        chains_var = tk.StringVar(value="")
        add_5_var = tk.BooleanVar(value=True)
        add_3_var = tk.BooleanVar(value=True)
        add_5_o3_var = tk.BooleanVar(value=False)
        last_default_output = {"path": ""}

        def default_output_path() -> Path:
            input_text = input_var.get().strip()
            if input_text:
                try:
                    return default_add_phosphate_output_path(Path(input_text).expanduser())
                except Exception:
                    pass
            return current_output_dir() / "add_phosphates_out.pdb"

        def refresh_default_output(force: bool = False) -> None:
            default_path = default_output_path()
            default_text = str(default_path)
            current = output_var.get().strip()
            if force or not current or current == last_default_output["path"]:
                output_var.set(default_text)
            last_default_output["path"] = default_text

        def set_dialog_log(text: str) -> None:
            result_text.configure(state="normal")
            result_text.delete("1.0", "end")
            result_text.insert("1.0", text)
            result_text.see("1.0")
            result_text.configure(state="disabled")
            win.update_idletasks()

        def refresh_info() -> None:
            pdb_text = input_var.get().strip()
            if not pdb_text:
                set_dialog_log("Choose an input PDB to list terminal phosphate status by chain.")
                return
            refresh_default_output()
            try:
                pdb_path = Path(pdb_text).expanduser()
                if not pdb_path.is_absolute():
                    pdb_path = pdb_path.resolve()
                statuses = analyze_phosphate_termini(pdb_path)
                set_dialog_log(format_phosphate_report(statuses))
            except Exception as exc:
                set_dialog_log(f"Could not read terminal phosphate status:\n{exc}")

        def browse_input() -> None:
            selected = filedialog.askopenfilename(
                title="Choose input PDB",
                initialdir=str(current_output_dir()),
                filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                input_var.set(selected)
                refresh_default_output(force=True)
                refresh_info()

        def browse_output() -> None:
            refresh_default_output()
            current = Path(output_var.get().strip() or str(default_output_path())).expanduser()
            selected = filedialog.asksaveasfilename(
                title="Save phosphate-added PDB",
                initialdir=str(current.parent if str(current.parent) not in ("", ".") else current_output_dir()),
                initialfile=current.name or "add_phosphates_out.pdb",
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                output_var.set(selected)

        def ends_arg() -> str:
            add_5 = bool(add_5_var.get())
            add_3 = bool(add_3_var.get())
            if add_5 and add_3:
                return "both"
            if add_5:
                return "5"
            if add_3:
                return "3"
            return "none"

        def run_add_phosphates() -> None:
            try:
                input_text = input_var.get().strip()
                if not input_text:
                    raise ValueError("Please choose an input PDB file.")
                input_path = Path(input_text).expanduser()
                if not input_path.is_absolute():
                    input_path = input_path.resolve()
                if not input_path.exists():
                    raise ValueError(f"Input PDB not found: {input_path}")

                refresh_default_output()
                output_text = output_var.get().strip()
                if not output_text:
                    raise ValueError("Please choose an output PDB file.")
                output_path = Path(output_text).expanduser()
                if not output_path.is_absolute():
                    output_path = output_path.resolve()

                selected_chains = parse_chain_selection(chains_var.get())
                end_text = ends_arg()
                add_5 = end_text in ("both", "5")
                add_3 = end_text in ("both", "3")
                add_5_o3 = bool(add_5_o3_var.get())
                if not add_5 and not add_3 and not add_5_o3:
                    raise ValueError(
                        "Select at least one item to add: 5' phosphate, 3' phosphate, or preceding O3'."
                    )

                run_btn.configure(state="disabled")
                win.update_idletasks()
                result = add_terminal_phosphates(
                    input_path,
                    output_path,
                    chain_ids=selected_chains,
                    add_5prime=add_5,
                    add_3prime=add_3,
                    add_5prime_o3=add_5_o3,
                )
            except Exception as exc:
                messagebox.showerror("Add phosphates", str(exc), parent=win)
                return
            finally:
                try:
                    run_btn.configure(state="normal")
                except Exception:
                    pass

            cli_args = [
                "python3",
                "bnp_na_lib/add_phosphates.py",
                str(input_path),
                "-o",
                str(output_path),
                "--ends",
                end_text,
            ]
            if selected_chains:
                cli_args.extend(["--chains", chains_var.get().strip()])
            if add_5_o3:
                cli_args.append("--add-5prime-o3")
            log = "\n".join(
                [
                    "=== Add phosphates tool ===",
                    "Equivalent CLI command:",
                    "  " + " ".join(shlex.quote(part) for part in cli_args),
                    "",
                    result.log_text,
                ]
            )
            set_dialog_log(result.log_text)
            self._append_log(log)
            messagebox.showinfo("Add phosphates", f"Wrote:\n{result.output_pdb}", parent=win)

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Input PDB:", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=label_pad, **field)
        input_entry = ttk.Entry(body, textvariable=input_var)
        input_entry.grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse/load", command=browse_input).grid(
            row=0, column=2, sticky="we", padx=(8, 0), **field
        )
        input_entry.bind("<Return>", lambda *_args: refresh_info())
        input_entry.bind("<FocusOut>", lambda *_args: refresh_info())

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=output_var).grid(row=1, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_output).grid(row=1, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Chain IDs:", style="Field.TLabel").grid(row=2, column=0, sticky="w", padx=label_pad, **field)
        chain_frame = ttk.Frame(body)
        chain_frame.grid(row=2, column=1, columnspan=2, sticky="w", **field)
        ttk.Entry(chain_frame, textvariable=chains_var, width=18).pack(side="left")
        ttk.Label(chain_frame, text="Blank = all chains", style="Hint.TLabel").pack(side="left", padx=(12, 0))

        ends_frame = ttk.Frame(body)
        ends_frame.grid(row=3, column=1, columnspan=2, sticky="w", pady=(4, 12))
        ttk.Label(body, text="Add:", style="Field.TLabel").grid(row=3, column=0, sticky="w", padx=label_pad, pady=(4, 12))
        ttk.Checkbutton(ends_frame, text="5' phosphates", variable=add_5_var).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(ends_frame, text="3' phosphates", variable=add_3_var).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(ends_frame, text="O3' before 5' phosphate", variable=add_5_o3_var).pack(side="left")

        result_view = ui.LogView(body, wrap="none", height=16, horizontal=True)
        result_view.grid(row=4, column=0, columnspan=3, sticky="nsew")
        result_text = result_view.text
        result_text.configure(state="disabled")

        footer.add_button("Refresh info", refresh_info)
        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Add phosphates", run_add_phosphates, primary=True)

        refresh_default_output(force=True)
        refresh_info()

    def open_xyz_bild_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Write XYZ axes BILD")
        win.geometry("680x372+240+160")
        win.minsize(600, 360)
        win.transient(self)
        body, footer = self._tool_dialog_frame(win, "Write XYZ axes BILD", tool="xyz_bild")
        body.grid_columnconfigure(1, weight=1)

        try:
            out_dir = Path(self.output_dir_var.get()).expanduser()
            if not out_dir.is_absolute():
                out_dir = out_dir.resolve()
        except Exception:
            out_dir = DEFAULT_OUTPUT_DIR

        out_var = tk.StringVar(value=str(out_dir / "xyz_axes.bild"))
        origin_var = tk.StringVar(value="0 0 0")
        length_var = tk.StringVar(value="20")
        width_var = tk.StringVar(value="1")
        sphere_var = tk.StringVar(value="0.5")

        def browse_out() -> None:
            path = filedialog.asksaveasfilename(
                title="Save XYZ axes BILD",
                initialfile="xyz_axes.bild",
                defaultextension=".bild",
                filetypes=[("BILD files", "*.bild"), ("All files", "*.*")],
                parent=win,
            )
            if path:
                out_var.set(path)

        def parse_origin(text: str) -> Tuple[float, float, float]:
            if not text.strip():
                return 0.0, 0.0, 0.0
            parts = [p.strip() for p in text.split(",")] if "," in text else text.split()
            if len(parts) != 3:
                raise ValueError("Origin must have three numbers, for example: 0 0 0")
            return (
                _parse_float_expression(parts[0], "Origin x"),
                _parse_float_expression(parts[1], "Origin y"),
                _parse_float_expression(parts[2], "Origin z"),
            )

        def write_file() -> None:
            try:
                out_path = Path(out_var.get()).expanduser()
                if not str(out_path).strip():
                    raise ValueError("Please choose an output .bild file.")
                origin = parse_origin(origin_var.get())
                length = _parse_float_expression(length_var.get(), "Arrow length", default=20.0)
                width = _parse_float_expression(width_var.get(), "Arrow width", default=1.0)
                sphere_radius = _parse_float_expression(sphere_var.get(), "Origin sphere radius", default=0.5)
                written = write_xyz_bild(
                    out_path,
                    origin=origin,
                    length=length,
                    width=width,
                    sphere_radius=sphere_radius,
                )
            except Exception as exc:
                messagebox.showerror("XYZ axes BILD", str(exc), parent=win)
                return
            cli_args = [
                "python3",
                "bnp_na_lib/xyz_bild.py",
                "--out",
                str(written),
                "--origin",
                " ".join(_format_number(value) for value in origin),
                "--length",
                _format_number(length),
                "--width",
                _format_number(width),
                "--sphere-radius",
                _format_number(sphere_radius),
            ]
            self._append_log(
                "\n".join(
                    [
                        "=== XYZ axes BILD tool ===",
                        "Equivalent CLI command:",
                        "  " + " ".join(shlex.quote(part) for part in cli_args),
                        "",
                        f"Wrote XYZ axes BILD: {written}",
                        f"Origin (Å): {origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f}",
                        f"Arrow length (Å): {length:.4f}",
                        f"Arrow width (Å): {width:.4f}",
                        f"Origin sphere radius (Å): {sphere_radius:.4f}",
                    ]
                )
            )
            messagebox.showinfo("XYZ axes BILD", f"Wrote:\n{written}", parent=win)
            win.destroy()

        field = {"pady": 4}
        label_pad = (0, 14)
        ttk.Label(body, text="Output .bild:", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=out_var).grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_out).grid(row=0, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Origin x y z (Å):", style="Field.TLabel").grid(
            row=1, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=origin_var, width=18).grid(row=1, column=1, sticky="w", **field)

        ttk.Label(body, text="Arrow length (Å):", style="Field.TLabel").grid(
            row=2, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=length_var, width=12).grid(row=2, column=1, sticky="w", **field)

        ttk.Label(body, text="Arrow width (Å):", style="Field.TLabel").grid(
            row=3, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=width_var, width=12).grid(row=3, column=1, sticky="w", **field)

        ttk.Label(body, text="Origin sphere radius (Å):", style="Field.TLabel").grid(
            row=4, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=sphere_var, width=12).grid(row=4, column=1, sticky="w", **field)

        hint = ttk.Label(
            body,
            text="Colors: X red, Y yellow, Z blue. Arrow head radius is 2.5 x arrow width.",
            style="Hint.TLabel",
            wraplength=460,
            justify="left",
        )
        hint.grid(row=5, column=0, columnspan=3, sticky="w", pady=(10, 0))
        self._wrap_to_width(body, hint)

        footer.add_button("Cancel", win.destroy)
        footer.add_button("Write", write_file, primary=True)

    def open_bz_builder_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("B-Z structure builder")
        win.geometry("940x860+180+60")
        win.minsize(820, 720)
        win.transient(self)
        body, footer = self._tool_dialog_frame(win, "B-Z builder", tool="bz_builder")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(1, weight=1)
        body.grid_rowconfigure(8, weight=2)

        try:
            out_dir = Path(self.output_dir_var.get()).expanduser()
            if not out_dir.is_absolute():
                out_dir = out_dir.resolve()
        except Exception:
            out_dir = DEFAULT_OUTPUT_DIR

        files_view = ui.LogView(body, wrap="none", height=5, horizontal=True)
        files_text = files_view.text
        files_text.configure(background=ui.C.BG)
        out_var = tk.StringVar(value=str(out_dir / "make_BZ_out.pdb"))
        axis_mode_var = tk.StringVar(value="codirectional")
        axis_source_var = tk.StringVar(value="auto")
        auto_trim_var = tk.BooleanVar(value=True)
        out_default_mode = {"value": True}

        def _files() -> list[str]:
            return [line.strip() for line in files_text.get("1.0", "end").splitlines() if line.strip()]

        def _input_default_dir() -> str:
            files = _files()
            if files:
                try:
                    return str(Path(files[0]).expanduser().resolve().parent)
                except Exception:
                    return str(out_dir)
            return str(out_dir)

        def _out_is_default_like() -> bool:
            current = out_var.get().strip()
            if not current:
                return True
            return Path(current).name == "make_BZ_out.pdb" and out_default_mode["value"]

        def _update_default_output_from_inputs() -> None:
            files = _files()
            if files and _out_is_default_like():
                try:
                    out_var.set(str(Path(files[0]).expanduser().resolve().parent / "make_BZ_out.pdb"))
                    out_default_mode["value"] = True
                except Exception:
                    pass

        def add_files() -> None:
            selected = filedialog.askopenfilenames(
                title="Choose B/Z input PDB files in alternating order",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if not selected:
                return
            current = files_text.get("1.0", "end").strip()
            added = "\n".join(selected)
            files_text.delete("1.0", "end")
            files_text.insert("1.0", (current + "\n" + added).strip() if current else added)
            _update_default_output_from_inputs()

        def clear_files() -> None:
            files_text.delete("1.0", "end")

        def browse_out() -> None:
            current_out = Path(out_var.get().strip() or "make_BZ_out.pdb").expanduser()
            initialfile = current_out.name or "make_BZ_out.pdb"
            initialdir = str(current_out.parent) if str(current_out.parent) not in ("", ".") else _input_default_dir()
            selected = filedialog.asksaveasfilename(
                title="Choose B-Z output PDB",
                initialdir=initialdir,
                initialfile=initialfile,
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                out_default_mode["value"] = False
                out_var.set(selected)

        field = {"pady": 4}
        label_pad = (0, 14)
        files_hint = ttk.Label(
            body,
            text="Input PDB files must alternate by order: B1, Z1, B2, Z2, ...",
            style="Hint.TLabel",
            wraplength=820,
            justify="left",
        )
        files_hint.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        files_view.grid(row=1, column=0, columnspan=3, sticky="nsew")

        file_buttons = ttk.Frame(body)
        file_buttons.grid(row=2, column=0, columnspan=3, sticky="w", pady=(8, 12))
        ttk.Button(file_buttons, text="Add input PDBs", command=add_files).pack(side="left", padx=(0, 8))
        ttk.Button(file_buttons, text="Clear", command=clear_files).pack(side="left")

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=3, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=out_var).grid(row=3, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_out).grid(row=3, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Axis correction:", style="Field.TLabel").grid(
            row=4, column=0, sticky="w", padx=label_pad, **field
        )
        axis_mode_row = ttk.Frame(body)
        axis_mode_row.grid(row=4, column=1, columnspan=2, sticky="we", **field)
        ttk.Combobox(
            axis_mode_row,
            textvariable=axis_mode_var,
            values=("codirectional", "collinear", "none"),
            state="readonly",
            width=16,
        ).pack(side="left")
        ttk.Label(
            axis_mode_row,
            text="Default codirectional rotates helix axes parallel without lateral axis-line shift.",
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        ttk.Label(body, text="Axis source:", style="Field.TLabel").grid(row=5, column=0, sticky="w", padx=label_pad, **field)
        axis_source_row = ttk.Frame(body)
        axis_source_row.grid(row=5, column=1, columnspan=2, sticky="we", **field)
        ttk.Combobox(
            axis_source_row,
            textvariable=axis_source_var,
            values=("auto", "dssr", "pca"),
            state="readonly",
            width=16,
        ).pack(side="left")
        ttk.Label(
            axis_source_row,
            text="auto tries DSSR --more through align2z.py, then falls back to C1' PCA.",
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        ttk.Checkbutton(
            body,
            text=(
                "Auto-trim terminal Z-DNA bp if needed. For DSSR/bnp_na Z inputs, "
                "prepare Z-DNA 2 bp longer than the target final Z segment."
            ),
            variable=auto_trim_var,
        ).grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 2))

        builder_hint = ttk.Label(
            body,
            text=(
                "The builder inserts B-Z junction cores, trims overlapping base pairs, and writes both a final ligated "
                "PDB plus a raw aligned PDB with _raw before the file extension."
            ),
            style="Hint.TLabel",
            wraplength=820,
            justify="left",
        )
        builder_hint.grid(row=7, column=0, columnspan=3, sticky="w", pady=(2, 12))
        self._wrap_to_width(body, builder_hint)

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Build B-Z structure", primary=True)

        bz_log_view = ui.LogView(body, wrap="none", height=10, horizontal=True)
        bz_log_view.grid(row=8, column=0, columnspan=3, sticky="nsew")
        bz_log = bz_log_view.text

        def _set_dialog_log(text: str) -> None:
            bz_log.delete("1.0", "end")
            bz_log.insert("1.0", text)
            bz_log.see("end")
            win.update_idletasks()

        def run_bz() -> None:
            files = _files()
            if len(files) < 2:
                messagebox.showerror("B-Z structure builder", "Please provide at least two files: B1 and Z1.", parent=win)
                return
            if _out_is_default_like():
                _update_default_output_from_inputs()
            out_text = out_var.get().strip() or str(Path(_input_default_dir()) / "make_BZ_out.pdb")

            run_btn.configure(state="disabled")
            start_log = (
                f"=== {APP_NAME} {__version__} B-Z job started ===\n"
                f"Input files: {len(files)}\n"
                f"Output PDB: {out_text}\n"
                f"Axis correction: {axis_mode_var.get()}\n"
                f"Axis source: {axis_source_var.get()}\n"
                f"Auto-trim Z-DNA terminals: {'ON' if auto_trim_var.get() else 'OFF'}\n"
            )
            _set_dialog_log(start_log)
            self._set_log(start_log)

            try:
                result = build_bz_structure(
                    files,
                    out_text,
                    axis_mode=axis_mode_var.get(),
                    axis_source=axis_source_var.get(),
                    auto_trim_z=auto_trim_var.get(),
                )
                pdb_out = Path(str(result["pdb_out"]))
                remarks_log = _prepend_final_pdb_remarks(
                    pdb_out,
                    na_type="B-Z DNA",
                    l_form_enabled=False,
                    invrot_result=None,
                )
                full_log = "\n".join(
                    [
                        f"=== {APP_NAME} {__version__} B-Z job summary ===",
                        f"Final B-Z PDB: {pdb_out}",
                        f"Raw aligned PDB: {result['pdb_raw']}",
                        f"Axis correction: {result['axis_mode']}",
                        f"Axis source: {result['axis_source']}",
                        f"Auto-trim Z-DNA terminals: {'ON' if result['auto_trim_z'] else 'OFF'}",
                        "",
                        str(result.get("log_text", "")),
                        "",
                        remarks_log,
                    ]
                )
                _set_dialog_log(full_log)
                self._set_log(full_log)
                messagebox.showinfo(
                    "B-Z structure builder",
                    f"Final B-Z PDB:\n{pdb_out}\n\nRaw aligned PDB:\n{result['pdb_raw']}",
                    parent=win,
                )
            except PipelineError as exc:
                log_text = getattr(exc, "log_text", "") or str(exc)
                _set_dialog_log(log_text)
                self._set_log(log_text)
                messagebox.showerror("B-Z structure builder", str(exc), parent=win)
            finally:
                run_btn.configure(state="normal")

        run_btn.configure(command=run_bz)

    def open_bp_builder_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("B-P structure builder")
        win.geometry("940x860+200+70")
        win.minsize(820, 720)
        win.transient(self)
        body, footer = self._tool_dialog_frame(win, "B-P builder", tool="bp_builder")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(1, weight=1)
        body.grid_rowconfigure(9, weight=2)

        try:
            out_dir = Path(self.output_dir_var.get()).expanduser()
            if not out_dir.is_absolute():
                out_dir = out_dir.resolve()
        except Exception:
            out_dir = DEFAULT_OUTPUT_DIR

        default_name = "make_BP_out.pdb"
        files_view = ui.LogView(body, wrap="none", height=5, horizontal=True)
        files_text = files_view.text
        files_text.configure(background=ui.C.BG)
        out_var = tk.StringVar(value=str(out_dir / default_name))
        types_var = tk.StringVar(value="")
        axis_mode_var = tk.StringVar(value="codirectional")
        axis_source_var = tk.StringVar(value="auto")
        transition_var = tk.StringVar(value="2")
        out_default_mode = {"value": True}

        def _files() -> list[str]:
            return [line.strip() for line in files_text.get("1.0", "end").splitlines() if line.strip()]

        def _input_default_dir() -> str:
            files = _files()
            if files:
                try:
                    return str(Path(files[0]).expanduser().resolve().parent)
                except Exception:
                    return str(out_dir)
            return str(out_dir)

        def _out_is_default_like() -> bool:
            current = out_var.get().strip()
            if not current:
                return True
            return Path(current).name == default_name and out_default_mode["value"]

        def _update_default_output_from_inputs() -> None:
            files = _files()
            if files and _out_is_default_like():
                try:
                    out_var.set(str(Path(files[0]).expanduser().resolve().parent / default_name))
                    out_default_mode["value"] = True
                except Exception:
                    pass

        def add_files() -> None:
            selected = filedialog.askopenfilenames(
                title="Choose B-DNA and P-DNA input PDB files in construct order",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if not selected:
                return
            current = files_text.get("1.0", "end").strip()
            added = "\n".join(selected)
            files_text.delete("1.0", "end")
            files_text.insert("1.0", (current + "\n" + added).strip() if current else added)
            _update_default_output_from_inputs()

        def clear_files() -> None:
            files_text.delete("1.0", "end")

        def browse_out() -> None:
            current_out = Path(out_var.get().strip() or default_name).expanduser()
            initialfile = current_out.name or default_name
            initialdir = str(current_out.parent) if str(current_out.parent) not in ("", ".") else _input_default_dir()
            selected = filedialog.asksaveasfilename(
                title="Choose B-P output PDB",
                initialdir=initialdir,
                initialfile=initialfile,
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                out_default_mode["value"] = False
                out_var.set(selected)

        field = {"pady": 4}
        label_pad = (0, 14)
        files_hint = ttk.Label(
            body,
            text=(
                "Input PDB files in construct order, alternating B-DNA and P-DNA in either order: "
                "B1, P1, B2, ... or P1, B1, P2, ..."
            ),
            style="Hint.TLabel",
            wraplength=820,
            justify="left",
        )
        files_hint.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        self._wrap_to_width(body, files_hint)
        files_view.grid(row=1, column=0, columnspan=3, sticky="nsew")

        file_buttons = ttk.Frame(body)
        file_buttons.grid(row=2, column=0, columnspan=3, sticky="w", pady=(8, 12))
        ttk.Button(file_buttons, text="Add input PDBs", command=add_files).pack(side="left", padx=(0, 8))
        ttk.Button(file_buttons, text="Clear", command=clear_files).pack(side="left")

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=3, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=out_var).grid(row=3, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_out).grid(row=3, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Helix types:", style="Field.TLabel").grid(row=4, column=0, sticky="w", padx=label_pad, **field)
        types_row = ttk.Frame(body)
        types_row.grid(row=4, column=1, columnspan=2, sticky="we", **field)
        ttk.Entry(types_row, textvariable=types_var, width=18).pack(side="left")
        ttk.Label(
            types_row,
            text="Blank detects B or P from the phosphate radius; or type one letter per file, e.g. BPB.",
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        ttk.Label(body, text="Axis correction:", style="Field.TLabel").grid(
            row=5, column=0, sticky="w", padx=label_pad, **field
        )
        axis_mode_row = ttk.Frame(body)
        axis_mode_row.grid(row=5, column=1, columnspan=2, sticky="we", **field)
        ttk.Combobox(
            axis_mode_row,
            textvariable=axis_mode_var,
            values=("codirectional", "collinear", "none"),
            state="readonly",
            width=16,
        ).pack(side="left")
        ttk.Label(
            axis_mode_row,
            text="Default codirectional rotates helix axes parallel without lateral axis-line shift.",
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        ttk.Label(body, text="B-DNA axis source:", style="Field.TLabel").grid(
            row=6, column=0, sticky="w", padx=label_pad, **field
        )
        axis_source_row = ttk.Frame(body)
        axis_source_row.grid(row=6, column=1, columnspan=2, sticky="we", **field)
        ttk.Combobox(
            axis_source_row,
            textvariable=axis_source_var,
            values=("auto", "dssr", "screw", "pca"),
            state="readonly",
            width=16,
        ).pack(side="left")
        ttk.Label(
            axis_source_row,
            text=(
                "auto tries DSSR --more, then the helix's screw axis. P-DNA axes always come from its screw "
                "symmetry."
            ),
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        ttk.Label(body, text="Transition:", style="Field.TLabel").grid(
            row=7, column=0, sticky="w", padx=label_pad, **field
        )
        transition_row = ttk.Frame(body)
        transition_row.grid(row=7, column=1, columnspan=2, sticky="we", **field)
        ttk.Combobox(
            transition_row,
            textvariable=transition_var,
            values=("2", "3"),
            state="readonly",
            width=16,
        ).pack(side="left")
        ttk.Label(
            transition_row,
            text=(
                "Nucleotides per strand rebuilt at each junction. 2: compact core. 3: every backbone step on "
                "both strands turns forward, between the B-DNA and P-DNA twists."
            ),
            style="Hint.TLabel",
            wraplength=520,
            justify="left",
        ).pack(side="left", padx=(12, 0))

        hint_row = ttk.Frame(body)
        hint_row.grid(row=8, column=0, columnspan=3, sticky="we", pady=(6, 12))
        ui.help_button(hint_row, command=self.open_bp_builder_help).pack(side="left", padx=(0, 10))
        ttk.Label(
            hint_row,
            text=(
                "Each junction inserts a modeled B-P core and keeps every input base. The ? explains the two "
                "transitions, how levels are counted and what is written."
            ),
            style="Hint.TLabel",
            wraplength=780,
            justify="left",
        ).pack(side="left", fill="x", expand=True)

        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Build B-P structure", primary=True)

        bp_log_view = ui.LogView(body, wrap="none", height=10, horizontal=True)
        bp_log_view.grid(row=9, column=0, columnspan=3, sticky="nsew")
        bp_log = bp_log_view.text

        def _set_dialog_log(text: str) -> None:
            bp_log.delete("1.0", "end")
            bp_log.insert("1.0", text)
            bp_log.see("end")
            win.update_idletasks()

        def run_bp() -> None:
            files = _files()
            if len(files) < 2:
                messagebox.showerror(
                    "B-P structure builder", "Please provide at least two files: one B-DNA and one P-DNA.", parent=win
                )
                return
            if _out_is_default_like():
                _update_default_output_from_inputs()
            out_text = out_var.get().strip() or str(Path(_input_default_dir()) / default_name)

            run_btn.configure(state="disabled")
            start_log = (
                f"=== {APP_NAME} {__version__} B-P job started ===\n"
                f"Input files: {len(files)}\n"
                f"Output PDB: {out_text}\n"
                f"Helix types: {types_var.get().strip() or 'auto-detect'}\n"
                f"Axis correction: {axis_mode_var.get()}\n"
                f"B-DNA axis source: {axis_source_var.get()}\n"
                f"Transition: {transition_var.get()} nt per strand\n"
            )
            _set_dialog_log(start_log)
            self._set_log(start_log)

            try:
                result = build_bp_structure(
                    files,
                    out_text,
                    types=types_var.get().strip(),
                    axis_mode=axis_mode_var.get(),
                    axis_source=axis_source_var.get(),
                    transition=transition_var.get(),
                )
                pdb_out = Path(str(result["pdb_out"]))
                remarks_log = _prepend_final_pdb_remarks(
                    pdb_out,
                    na_type="B-P DNA",
                    l_form_enabled=False,
                    invrot_result=None,
                )
                full_log = "\n".join(
                    [
                        f"=== {APP_NAME} {__version__} B-P job summary ===",
                        f"Final B-P PDB: {pdb_out}",
                        f"Raw aligned PDB: {result['pdb_raw']}",
                        f"Helix types: {result['types']}",
                        f"Axis correction: {result['axis_mode']}",
                        f"B-DNA axis source: {result['axis_source']}",
                        f"Transition: {result['transition']} nt per strand",
                        "",
                        str(result.get("log_text", "")),
                        "",
                        remarks_log,
                    ]
                )
                _set_dialog_log(full_log)
                self._set_log(full_log)
                messagebox.showinfo(
                    "B-P structure builder",
                    f"Final B-P PDB:\n{pdb_out}\n\nRaw aligned PDB:\n{result['pdb_raw']}",
                    parent=win,
                )
            except PipelineError as exc:
                log_text = getattr(exc, "log_text", "") or str(exc)
                _set_dialog_log(log_text)
                self._set_log(log_text)
                messagebox.showerror("B-P structure builder", str(exc), parent=win)
            finally:
                run_btn.configure(state="normal")

        run_btn.configure(command=run_bp)

    def open_triplex_converter_dialog(self) -> None:
        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title("Triplex converter")
        win.geometry("960x860+190+50")
        win.minsize(860, 780)
        win.transient(self)
        body, footer = self._tool_dialog_frame(win, "Triplex converter", tool="triplex")
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(10, weight=1)

        try:
            out_dir = Path(self.output_dir_var.get()).expanduser()
            if not out_dir.is_absolute():
                out_dir = out_dir.resolve()
        except Exception:
            out_dir = DEFAULT_OUTPUT_DIR

        duplex_var = tk.StringVar(value="")
        out_var = tk.StringVar(value=str(out_dir / "triplex_2TH.pdb"))
        strand_i_var = tk.StringVar(value="A")
        range_start_var = tk.StringVar(value="1")
        range_end_var = tk.StringVar(value="1")
        mode_var = tk.StringVar(value="antiparallel")
        strand_ii_var = tk.StringVar(value="")
        strand_iii_var = tk.StringVar(value="")
        strand_iii_start_var = tk.StringVar(value="1")
        out_default_mode = {"value": True}

        triplex_log_view = ui.LogView(body, wrap="word", height=9)
        triplex_log = triplex_log_view.text

        def _set_dialog_log(text: str) -> None:
            triplex_log.delete("1.0", "end")
            triplex_log.insert("1.0", text)
            triplex_log.see("end")
            win.update_idletasks()

        def _current_range() -> Tuple[int, int]:
            return (
                _parse_int_expression(range_start_var.get(), "Triplex range start"),
                _parse_int_expression(range_end_var.get(), "Triplex range end"),
            )

        def _default_output_for_input() -> Optional[Path]:
            text = duplex_var.get().strip()
            if not text:
                return None
            try:
                return default_triplex_output_path(text)
            except Exception:
                return None

        def _update_default_output(force: bool = False) -> None:
            default_out = _default_output_for_input()
            if default_out is None:
                return
            current = out_var.get().strip()
            if force or not current or out_default_mode["value"]:
                out_var.set(str(default_out))
                out_default_mode["value"] = True

        def browse_input() -> None:
            selected = filedialog.askopenfilename(
                title="Choose input duplex PDB",
                initialdir=str(out_dir),
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                duplex_var.set(selected)
                _update_default_output(force=True)
                refresh_info()

        def browse_out() -> None:
            _update_default_output(force=False)
            current = Path(out_var.get().strip() or str(out_dir / "triplex_2TH.pdb")).expanduser()
            selected = filedialog.asksaveasfilename(
                title="Choose triplex output PDB",
                initialdir=str(current.parent if str(current.parent) not in ("", ".") else out_dir),
                initialfile=current.name or "triplex_2TH.pdb",
                defaultextension=".pdb",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
                parent=win,
            )
            if selected:
                out_default_mode["value"] = False
                out_var.set(selected)

        def refresh_info() -> None:
            duplex = duplex_var.get().strip()
            if not duplex:
                _set_dialog_log(
                    "Choose an input duplex PDB to list detected chains and sequences.\n"
                    "Default output naming inserts _2TH before the input file extension."
                )
                return
            _update_default_output(force=False)
            try:
                text = describe_triplex_input(
                    duplex,
                    strand_i_chain=strand_i_var.get().strip(),
                    residue_range=_current_range(),
                    mode=mode_var.get().strip(),
                )
            except Exception as exc:
                text = f"Could not read triplex input information:\n{exc}"
            _set_dialog_log(text)

        field = {"pady": 3}
        label_pad = (0, 14)
        ttk.Label(body, text="Input duplex PDB:", style="Field.TLabel").grid(
            row=0, column=0, sticky="w", padx=label_pad, **field
        )
        duplex_entry = ttk.Entry(body, textvariable=duplex_var)
        duplex_entry.grid(row=0, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse/load", command=browse_input).grid(
            row=0, column=2, sticky="we", padx=(8, 0), **field
        )
        duplex_entry.bind("<Return>", lambda *_args: refresh_info())
        duplex_entry.bind("<FocusOut>", lambda *_args: refresh_info())

        ttk.Label(body, text="Output PDB:", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=label_pad, **field)
        ttk.Entry(body, textvariable=out_var).grid(row=1, column=1, sticky="we", **field)
        ttk.Button(body, text="Browse", command=browse_out).grid(row=1, column=2, sticky="we", padx=(8, 0), **field)

        ttk.Label(body, text="Strand I purine chain ID:", style="Field.TLabel").grid(
            row=2, column=0, sticky="w", padx=label_pad, pady=(10, 3)
        )
        strand_i_entry = ttk.Entry(body, textvariable=strand_i_var, width=10)
        strand_i_entry.grid(row=2, column=1, sticky="w", pady=(10, 3))
        strand_i_entry.bind("<Return>", lambda *_args: refresh_info())
        strand_i_entry.bind("<FocusOut>", lambda *_args: refresh_info())

        range_frame = ttk.Frame(body)
        range_frame.grid(row=3, column=1, sticky="w", **field)
        ttk.Label(body, text="Strand I residue range:", style="Field.TLabel").grid(
            row=3, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(range_frame, textvariable=range_start_var, width=10).pack(side="left")
        ttk.Label(range_frame, text=" to ", style="Field.TLabel").pack(side="left", padx=4)
        ttk.Entry(range_frame, textvariable=range_end_var, width=10).pack(side="left")

        ttk.Label(body, text="Triplex mode:", style="Field.TLabel").grid(row=4, column=0, sticky="w", padx=label_pad, **field)
        mode_row = ttk.Frame(body)
        mode_row.grid(row=4, column=1, columnspan=2, sticky="w", **field)
        mode_combo = ttk.Combobox(
            mode_row,
            textvariable=mode_var,
            values=("antiparallel", "parallel"),
            state="readonly",
            width=16,
        )
        mode_combo.pack(side="left")
        mode_combo.bind("<<ComboboxSelected>>", lambda *_args: refresh_info())
        ttk.Label(
            mode_row,
            text="antiparallel = G·G-C; parallel = T·A-T.",
            style="Hint.TLabel",
            wraplength=300,
        ).pack(side="left", padx=(12, 0))

        ttk.Label(body, text="Strand II chain ID (optional):", style="Field.TLabel").grid(
            row=5, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=strand_ii_var, width=10).grid(row=5, column=1, sticky="w", **field)

        ttk.Label(body, text="Strand III chain ID (optional):", style="Field.TLabel").grid(
            row=6, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=strand_iii_var, width=10).grid(row=6, column=1, sticky="w", **field)

        ttk.Label(body, text="Strand III first resSeq:", style="Field.TLabel").grid(
            row=7, column=0, sticky="w", padx=label_pad, **field
        )
        ttk.Entry(body, textvariable=strand_iii_start_var, width=10).grid(row=7, column=1, sticky="w", **field)

        convention_hint = ttk.Label(
            body,
            text=(
                "Convention: triples are Z·X-Y. Strand I is the purine duplex strand X; "
                "strand II is the Watson-Crick partner Y; strand III is the added Hoogsteen strand Z. "
                "Only the selected strand-I range must be all G for antiparallel or all A for parallel."
            ),
            style="Hint.TLabel",
            wraplength=850,
            justify="left",
        )
        convention_hint.grid(row=8, column=0, columnspan=3, sticky="w", pady=(6, 12))
        self._wrap_to_width(body, convention_hint)

        min_frame = ui.Card(body, "Post-processing")
        min_frame.grid(row=9, column=0, columnspan=3, sticky="we", pady=(0, 12))
        min_frame.grid_columnconfigure(2, weight=1)

        triplex_minimize_var = tk.BooleanVar(value=True)
        triplex_regularize_var = tk.BooleanVar(value=True)
        triplex_params_var = tk.StringVar(value=self.params_var.get().strip() or _path_text(DEFAULT_PARAMS_FILE))

        min_check = ttk.Checkbutton(
            min_frame,
            text="Run phenix.geometry_minimization",
            variable=triplex_minimize_var,
        )
        min_check.grid(row=0, column=0, sticky="w", pady=4)
        ttk.Label(min_frame, text="Params file (.eff / .params):", style="Field.TLabel").grid(
            row=0, column=1, sticky="w", padx=(18, 10), pady=4
        )
        params_entry = ttk.Entry(min_frame, textvariable=triplex_params_var)
        params_entry.grid(row=0, column=2, sticky="we", pady=4)

        def browse_triplex_params() -> None:
            path = filedialog.askopenfilename(
                parent=win,
                title="Select params file",
                filetypes=[("Params files", "*.eff *.params"), ("All files", "*.*")],
            )
            if path:
                triplex_params_var.set(path)

        ttk.Button(min_frame, text="Browse", command=browse_triplex_params).grid(
            row=0, column=3, sticky="w", padx=(8, 0), pady=4
        )

        ttk.Checkbutton(
            min_frame,
            text="Regularize phosphates",
            variable=triplex_regularize_var,
        ).grid(row=1, column=0, sticky="w", pady=(2, 4))

        def on_triplex_minimize_toggled() -> None:
            state = "normal" if triplex_minimize_var.get() else "disabled"
            params_entry.configure(state=state)

        min_check.configure(command=on_triplex_minimize_toggled)

        post_hint = ttk.Label(
            min_frame,
            text=(
                "Applied to the converted triplex. The output PDB above receives the final model; "
                f"the raw conversion and every intermediate go into a `{TRIPLEX_TMP_DIR_NAME}` "
                "folder beside it."
            ),
            style="Hint.TLabel",
            wraplength=850,
            justify="left",
        )
        post_hint.grid(row=2, column=0, columnspan=4, sticky="w", pady=(4, 0))
        self._wrap_to_width(min_frame, post_hint, margin=36)

        footer.add_button("Refresh strand info", refresh_info)
        footer.add_button("Close", win.destroy)
        run_btn = footer.add_button("Convert to triplex", primary=True)

        triplex_log_view.grid(row=10, column=0, columnspan=3, sticky="nsew")

        def run_triplex() -> None:
            duplex = duplex_var.get().strip()
            _update_default_output(force=False)
            out_text = out_var.get().strip()
            run_phenix = bool(triplex_minimize_var.get())
            run_regularize = bool(triplex_regularize_var.get())
            try:
                residue_range = _current_range()
                strand_iii_start = _parse_int_expression(strand_iii_start_var.get(), "Strand III first resSeq")
                params_path = None
                if run_phenix:
                    params_text = triplex_params_var.get().strip()
                    if not params_text:
                        raise ValueError(
                            "Please specify a params file for phenix.geometry_minimization."
                        )
                    candidate = Path(params_text).expanduser()
                    if not candidate.is_absolute():
                        candidate = candidate.resolve()
                    if not candidate.exists():
                        raise ValueError(f"Params file not found: {candidate}")
                    params_path = str(candidate)
            except Exception as exc:
                messagebox.showerror("Triplex converter", str(exc), parent=win)
                return

            run_btn.configure(state="disabled")
            start_log = (
                f"=== {APP_NAME} {__version__} triplex job started ===\n"
                f"Input duplex PDB: {duplex}\n"
                f"Output PDB: {out_text or '<default _2TH path>'}\n"
                f"Mode: {mode_var.get()}\n"
                f"Strand I: {strand_i_var.get().strip()}\n"
                f"Residue range: {residue_range[0]}:{residue_range[1]}\n"
                f"phenix.geometry_minimization: {run_phenix}\n"
                f"Regularize phosphates: {run_regularize}\n"
            )
            _set_dialog_log(start_log)
            self._set_log(start_log)

            try:
                result = build_triplex_from_duplex(
                    duplex,
                    out_text,
                    strand_i_chain=strand_i_var.get().strip(),
                    residue_range=residue_range,
                    mode=mode_var.get(),
                    strand_ii_chain=(strand_ii_var.get().strip() or None),
                    strand_iii_chain=(strand_iii_var.get().strip() or None),
                    strand_iii_start_resseq=strand_iii_start,
                    run_phenix=run_phenix,
                    params_file=params_path,
                    run_regularize_phosphates=run_regularize,
                )
                pdb_out = Path(str(result["pdb_out"]))
                remarks_log = _prepend_final_pdb_remarks(
                    pdb_out,
                    na_type="Triplex DNA",
                    l_form_enabled=False,
                    invrot_result=None,
                )
                full_log = "\n".join(
                    [
                        f"=== {APP_NAME} {__version__} triplex job summary ===",
                        f"Triplex PDB: {pdb_out}",
                        f"Mode: {result['mode']} ({result['motif_label']})",
                        f"Strand I: {result['strand_i']}",
                        f"Strand II: {result['strand_ii']}",
                        f"Strand III: {result['strand_iii']}",
                        "",
                        str(result.get("log_text", "")),
                        "",
                        remarks_log,
                    ]
                )
                _set_dialog_log(full_log)
                self._set_log(full_log)
                if out_var.get().strip() == "":
                    out_var.set(str(pdb_out))
                messagebox.showinfo("Triplex converter", f"Triplex PDB:\n{pdb_out}", parent=win)
            except PipelineError as exc:
                log_text = getattr(exc, "log_text", "") or str(exc)
                _set_dialog_log(log_text)
                self._set_log(log_text)
                messagebox.showerror("Triplex converter", str(exc), parent=win)
            finally:
                run_btn.configure(state="normal")

        run_btn.configure(command=run_triplex)
        refresh_info()

    def _refresh_info_text(self) -> None:
        out = self.output_dir_var.get().strip() or "<not selected>"
        tmp = str(Path(out).expanduser() / "tmp_file") if out != "<not selected>" else "<not selected>"
        na_type = self.na_type_var.get().strip() if hasattr(self, "na_type_var") else ""
        if na_type == "P-DNA":
            pipeline = (
                "Pipeline: P-DNA from the 1PFI template (bases grafted per strand) -> normalize names -> "
                "optional phenix.geometry_minimization"
            )
        else:
            pipeline = (
                "Pipeline: build by selected type -> normalize names -> optional phenix.geometry_minimization "
                "for B-DNA/A-DNA/A-RNA/P-DNA"
            )
        if hasattr(self, "regularize_phosphates_var") and self.regularize_phosphates_var.get() and na_type != "Z-DNA":
            pipeline += " -> regularize phosphates from C1' helical symmetry"
        if na_type == "P-DNA":
            pipeline += " -> align to +Z by superposing onto the ideal template model (DSSR finds no P-DNA helix)"
        else:
            pipeline += " -> DSSR --more axis extraction -> align to +Z"
        if hasattr(self, "invrot_enabled_var") and self.invrot_enabled_var.get():
            pipeline += f" -> inv/rot mirror ({self.invrot_operation_var.get().strip() or 'i'})"
        pipeline += " -> orient/place."
        self.info_var.set(
            f"{pipeline}\n"
            f"Final placed PDB folder: {out}\n"
            f"Intermediate files folder: {tmp}"
        )

    def _current_defaults(self):
        return DEFAULT_PARAMS.get(self.na_type_var.get().strip())

    def _current_param_store(self) -> Optional[Dict[str, str]]:
        return self.param_values_by_type.get(self.na_type_var.get().strip())

    def _refresh_param_values_display(self) -> None:
        na_type = self.na_type_var.get().strip()
        defaults = self._current_defaults()
        store = self._current_param_store()
        if na_type == "P-DNA":
            self.param_frame.set_title("P-DNA helical parameters (fixed)")
            self.param_status_var.set(PDNA_PARAM_NOTE)
            self.param_table.delete(*self.param_table.get_children())
            self.param_table.grid_remove()
            self.pdna_param_table.grid()
            return
        self.param_frame.set_title("Current DSSR helical parameters")
        self.pdna_param_table.grid_remove()
        if defaults is None or store is None:
            self.param_status_var.set("Z-DNA uses DSSR fiber; the 12-parameter table is not used.")
            self.param_table.delete(*self.param_table.get_children())
            self.param_table.grid_remove()
            return

        self.param_table.grid()
        self.param_table.delete(*self.param_table.get_children())
        changed = []
        current_values = []
        default_values = []
        sources = []

        def shown_value_and_source(key: str, default: float) -> Tuple[str, str]:
            raw = store.get(key, "").strip()
            numeric = _parse_float_expression(raw, key, default=default)
            shown = _format_number(numeric)
            if raw and abs(numeric - float(default)) > 1e-10:
                changed.append(key)
                source = "custom"
            else:
                source = "default"
            return shown, source

        for key, default in zip(PARAM_KEYS, defaults):
            value, source = shown_value_and_source(key, default)
            current_values.append(value)
            default_values.append(_format_number(default))
            sources.append(source)

        self.param_table.insert("", "end", values=("Current", *current_values))
        self.param_table.insert("", "end", values=("Default", *default_values))
        self.param_table.insert("", "end", values=("Source", *sources))

        if changed:
            status = f"{na_type} table values; customized fields: {', '.join(changed)}."
        else:
            status = f"{na_type} table values; all fields use defaults."
        self.param_status_var.set(status)

    def _update_sequence_length(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type == "Z-DNA":
            text = self.z_len_var.get().strip()
            if not text:
                self.length_var.set("Length: 0 bp")
                return
            try:
                n = _parse_int_expression(text, "Z-DNA helix length")
                if n <= 0 or n % 2 != 0:
                    self.length_var.set("Length: invalid Z-DNA length")
                else:
                    self.length_var.set(f"Length: {n} bp; Z-DNA repeat = {n // 2}")
            except Exception:
                self.length_var.set("Length: invalid Z-DNA length")
            return

        seq_text = self.seq_var.get().strip()
        if na_type == "P-DNA":
            if not seq_text:
                self.length_var.set("Length: 0 nt per strand")
                return
            try:
                seq1 = expand_sequence(seq_text, alphabet="ATCG")
                seq2_text = self.seq2_var.get().strip()
                if seq2_text:
                    seq2 = expand_sequence(seq2_text, alphabet="ATCG")
                    if len(seq2) != len(seq1):
                        self.length_var.set(f"Length: strand I {len(seq1)} nt, strand II {len(seq2)} nt (must match)")
                        return
                resolve_pdna_sequences(seq_text, seq2_text)
                note = "; strand II = reverse complement" if not seq2_text else ""
                self.length_var.set(f"Length: {len(seq1)} nt per strand{note}")
            except Exception as exc:
                self.length_var.set(f"Length: invalid sequence ({exc})")
            return
        if not seq_text:
            self.length_var.set("Length: 0 bp")
            return
        try:
            seq = expand_sequence(seq_text, alphabet=sequence_alphabet(na_type))
            self.length_var.set(f"Length: {len(seq)} bp")
        except Exception as exc:
            self.length_var.set(f"Length: invalid sequence ({exc})")

    def open_param_dialog(self) -> None:
        na_type = self.na_type_var.get().strip()
        defaults = self._current_defaults()
        store = self._current_param_store()
        if na_type == "P-DNA":
            messagebox.showinfo("P-DNA", PDNA_PARAM_NOTE, parent=self)
            return
        if defaults is None or store is None:
            messagebox.showinfo("Z-DNA", "Z-DNA uses DSSR fiber; the 12-parameter table is not used.", parent=self)
            return

        win = tk.Toplevel(self)
        self._set_optional_window_icon(win)
        win.title(f"Customize DSSR parameters — {na_type}")
        win.geometry("980x616+180+100")
        win.minsize(900, 596)
        win.transient(self)
        win.grab_set()
        body, footer = self._tool_dialog_frame(win, f"Customize DSSR parameters — {na_type}")
        body.grid_columnconfigure(0, weight=1, uniform="param_groups")
        body.grid_columnconfigure(1, weight=1, uniform="param_groups")

        param_hint = ttk.Label(
            body,
            text=(
                "Leave a field blank to use the default value shown beside it. "
                "Previously saved custom values are kept. Values are written to four digits after the decimal."
            ),
            wraplength=900,
            style="Hint.TLabel",
            justify="left",
        )
        param_hint.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 10))
        self._wrap_to_width(body, param_hint)

        # The twelve fields in two groups of six, side by side.
        param_groups = []
        for group_index in range(2):
            group = ttk.Frame(body)
            group.grid(row=1, column=group_index, sticky="nw", padx=(0, 12) if group_index == 0 else (12, 0))
            param_groups.append(group)

        local_vars: Dict[str, tk.StringVar] = {}
        for idx, (key, default) in enumerate(zip(PARAM_KEYS, defaults)):
            group = param_groups[0 if idx < 6 else 1]
            row = idx if idx < 6 else idx - 6
            ttk.Label(group, text=PARAM_LABELS.get(key, key) + ":", style="Field.TLabel").grid(
                row=row, column=0, sticky="w", padx=(0, 14), pady=4
            )
            initial = store.get(key, "").strip()
            local_vars[key] = tk.StringVar(value=initial)
            ttk.Entry(group, textvariable=local_vars[key], width=18).grid(
                row=row, column=1, sticky="w", pady=4
            )
            ttk.Label(group, text=f"default {_format_number(default)} {PARAM_UNITS.get(key, '')}", style="Hint.TLabel").grid(
                row=row, column=2, sticky="w", padx=(12, 0), pady=4
            )

        self._build_xdisp_finder(win, na_type, defaults, local_vars, parent=body)

        def clear_all() -> None:
            for key in PARAM_KEYS:
                local_vars[key].set("")

        def restore_defaults() -> None:
            for key, default in zip(PARAM_KEYS, defaults):
                local_vars[key].set(_format_number(default))

        def save_close() -> None:
            for key in PARAM_KEYS:
                txt = local_vars[key].get().strip()
                if txt:
                    try:
                        _parse_float_expression(txt, key)
                    except Exception:
                        messagebox.showerror(
                            "Invalid value",
                            f"{PARAM_LABELS.get(key, key)} must be a number or simple arithmetic expression.",
                            parent=win,
                        )
                        return
                store[key] = txt
            self._refresh_param_values_display()
            win.destroy()

        footer.add_button("Clear fields", clear_all)
        footer.add_button("Restore defaults", restore_defaults)
        footer.add_button("Cancel", win.destroy)
        footer.add_button("Save", save_close, primary=True)

    def _build_xdisp_finder(
        self,
        win: tk.Toplevel,
        na_type: str,
        defaults: List[float],
        local_vars: Dict[str, tk.StringVar],
        parent: Optional[tk.Misc] = None,
    ) -> None:
        """Add the opposing-phosphate X-disp search panel to the parameter dialog.

        The search builds many helices, so it runs on a worker thread and reports
        back through a queue drained on the Tk thread. The panel is gridded into
        `parent` (the dialog's body frame), or into `win` when none is given.
        """
        host = parent if parent is not None else win
        # The card's title is its enable switch; give it the card-title font.
        ttk.Style(win).configure("CardTitle.TCheckbutton", font="AzbSection")
        frame = ttk.Labelframe(host, style="Card.TLabelframe", padding=(14, 6, 14, 12))
        frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(16, 0))
        frame.grid_columnconfigure(0, weight=1)

        # The enable switch is the group's own title, so ticking the box that
        # names the panel is what turns the panel on.
        enabled_var = tk.BooleanVar(value=False)
        enable_check = ttk.Checkbutton(
            frame,
            text="Opposing phosphate X-disp",
            variable=enabled_var,
            style="CardTitle.TCheckbutton",
        )
        frame.configure(labelwidget=enable_check)

        warn_label = ttk.Label(
            frame,
            text=(
                "Rarely needed. This is only useful for building the \"ideal helix\" presumed by many "
                "DNA nanostructure designs; ordinary helix building does not require it."
            ),
            wraplength=880,
            style="Warn.TLabel",
            justify="left",
        )
        warn_label.grid(row=0, column=0, sticky="w", pady=(4, 4))
        self._wrap_to_width(frame, warn_label, margin=36)

        search_hint = ttk.Label(
            frame,
            text=(
                "Searches for the X-disp that places opposing phosphate P atoms across the helix axis, "
                "holding every other parameter at the value shown above. The result is written into the "
                "X-disp field."
            ),
            wraplength=880,
            style="Hint.TLabel",
            justify="left",
        )
        search_hint.grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._wrap_to_width(frame, search_hint, margin=36)

        refine_var = tk.BooleanVar(
            value=bool(self.minimize_var.get() and self.regularize_phosphates_var.get())
        )
        refine_check = ttk.Checkbutton(
            frame,
            text="Refine with phenix.geometry_minimization + phosphate regularization",
            variable=refine_var,
            state="disabled",
        )
        refine_check.grid(row=2, column=0, sticky="w", pady=(0, 2))

        refine_hint = ttk.Label(
            frame,
            text=(
                "Off: raw DSSR rebuild only, a few seconds. "
                "On: matches the minimized build pipeline and takes a few minutes."
            ),
            wraplength=880,
            style="Hint.TLabel",
            justify="left",
        )
        refine_hint.grid(row=3, column=0, sticky="w", padx=(22, 0), pady=(0, 4))
        self._wrap_to_width(frame, refine_hint, margin=58)

        status_var = tk.StringVar(value="")
        # Tk drops a variable as soon as Python does, and on the disabled paths
        # below no closure keeps these alive. Anchor them to the dialog.
        frame.xdisp_vars = (enabled_var, refine_var, status_var)  # type: ignore[attr-defined]

        actions = ttk.Frame(frame)
        actions.grid(row=4, column=0, sticky="we", pady=(10, 0))
        run_btn = ttk.Button(actions, text="Find X-disp", state="disabled")
        stop_btn = ttk.Button(actions, text="Stop", state="disabled")
        run_btn.pack(side="left")
        stop_btn.pack(side="left", padx=(8, 0))
        ttk.Label(actions, textvariable=status_var, style="Hint.TLabel").pack(side="left", padx=(14, 0))

        if na_type != "B-DNA":
            enable_check.configure(state="disabled")
            status_var.set(f"Available for B-DNA only; {na_type} is selected.")
            return

        state: Dict[str, object] = {"worker": None, "cancel": False, "queue": queue.Queue()}

        def collect_fixed_overrides() -> Dict[str, float]:
            overrides: Dict[str, float] = {}
            for key, default in zip(PARAM_KEYS, defaults):
                if key == "X-disp":
                    continue
                text = local_vars[key].get().strip()
                if not text:
                    continue
                value = _parse_float_expression(text, PARAM_LABELS.get(key, key))
                if abs(value - float(default)) > 1e-10:
                    overrides[key] = value
            return overrides

        def set_running(running: bool) -> None:
            run_btn.configure(state="disabled" if running else "normal")
            stop_btn.configure(state="normal" if running else "disabled")
            refine_check.configure(state="disabled" if running else "normal")
            # Locked while a search is in flight, so enabled_var cannot change
            # underneath the run that is already going.
            enable_check.configure(state="disabled" if running else "normal")

        def apply_enabled_state() -> None:
            on = bool(enabled_var.get())
            refine_check.configure(state="normal" if on else "disabled")
            run_btn.configure(state="normal" if on else "disabled")
            if not on:
                status_var.set("")

        def finish(result: Dict[str, object]) -> None:
            self._append_log(format_result_report(result))
            x_value = float(result["recommended_x_disp"])
            local_vars["X-disp"].set(f"{x_value:.4f}")
            mode = "Phenix + regularized" if result["recommended_mode"] == "phenix" else "raw DSSR"
            status_var.set(f"X-disp set to {x_value:.4f} Å ({mode}).")

        def poll() -> None:
            finished = False
            try:
                while True:
                    kind, payload = state["queue"].get_nowait()  # type: ignore[attr-defined]
                    if kind == "log":
                        self._append_log(str(payload))
                        status_var.set(str(payload).strip()[:90])
                    elif kind == "done":
                        finish(payload)
                        finished = True
                    elif kind == "cancelled":
                        self._append_log("Opposing phosphate X-disp search cancelled.")
                        status_var.set("Cancelled.")
                        finished = True
                    elif kind == "error":
                        self._append_log(f"Opposing phosphate X-disp search failed: {payload}")
                        status_var.set("Search failed.")
                        finished = True
                        if win.winfo_exists():
                            messagebox.showerror("X-disp search", str(payload), parent=win)
            except queue.Empty:
                pass

            if not win.winfo_exists():
                return
            if finished:
                state["worker"] = None
                set_running(False)
                return
            win.after(150, poll)

        def start() -> None:
            if not enabled_var.get():
                return
            worker = state.get("worker")
            if worker is not None and worker.is_alive():  # type: ignore[union-attr]
                return
            try:
                seq = self.seq_var.get().strip()
                if not seq:
                    raise ValueError("Enter a sequence in the main window before searching.")
                fixed_overrides = collect_fixed_overrides()
                refine = bool(refine_var.get())
                params_file = self._validate_params_for_minimization(refine)
            except Exception as exc:
                messagebox.showerror("Input error", str(exc), parent=win)
                return

            state["cancel"] = False
            set_running(True)
            status_var.set("Searching...")
            self._append_log(
                "=== Opposing phosphate X-disp search started ===\n"
                f"Sequence: {seq}\n"
                f"Refine with Phenix + phosphate regularization: {refine}"
            )

            def work() -> None:
                try:
                    result = find_opposing_phosphate_xdisp(
                        seq,
                        params_file=params_file,
                        fixed_overrides=fixed_overrides,
                        refine_with_phenix=refine,
                        progress=lambda message: state["queue"].put(("log", message)),  # type: ignore[attr-defined]
                        should_cancel=lambda: bool(state["cancel"]),
                    )
                except OpposingPhosphateCancelled:
                    state["queue"].put(("cancelled", None))  # type: ignore[attr-defined]
                except Exception as exc:
                    state["queue"].put(("error", str(exc)))  # type: ignore[attr-defined]
                else:
                    state["queue"].put(("done", result))  # type: ignore[attr-defined]

            state["worker"] = threading.Thread(target=work, daemon=True)
            state["worker"].start()  # type: ignore[union-attr]
            win.after(150, poll)

        def stop() -> None:
            state["cancel"] = True
            status_var.set("Stopping after the current model...")

        enable_check.configure(command=apply_enabled_state)
        run_btn.configure(command=start)
        stop_btn.configure(command=stop)
        win.protocol("WM_DELETE_WINDOW", lambda: (state.__setitem__("cancel", True), win.destroy()))

    def _get_param_overrides(self) -> Dict[str, float]:
        defaults = self._current_defaults()
        store = self._current_param_store()
        if defaults is None or store is None:
            return {}
        overrides: Dict[str, float] = {}
        for key, default in zip(PARAM_KEYS, defaults):
            txt = store.get(key, "").strip()
            if txt:
                value = _parse_float_expression(txt, PARAM_LABELS.get(key, key))
                if abs(value - float(default)) > 1e-10:
                    overrides[key] = value
        return overrides

    def _get_placement_values(self) -> Dict[str, float]:
        return {
            "delta_z": _parse_float_expression(self.delta_z_var.get(), "delta_z (Å)", default=0.0),
            "x": _parse_float_expression(self.x_var.get(), "x (Å)", default=0.0),
            "y": _parse_float_expression(self.y_var.get(), "y (Å)", default=0.0),
            "z": _parse_float_expression(self.z_var.get(), "z (Å)", default=0.0),
            "roll": _parse_float_expression(self.roll_var.get(), "roll (°)", default=0.0),
            "phi": _parse_float_expression(self.phi_var.get(), "phi (°)", default=0.0),
            "theta": _parse_float_expression(self.theta_var.get(), "theta (°)", default=0.0),
        }

    def on_type_changed(self) -> None:
        old_type = self._active_na_type
        if old_type in self.minimize_by_type:
            self.minimize_by_type[old_type] = bool(self.minimize_var.get())
            self.regularize_phosphates_by_type[old_type] = bool(self.regularize_phosphates_var.get())

        na_type = self.na_type_var.get().strip()
        is_pdna = na_type == "P-DNA"
        self.seq_label.configure(text="Strand I (5'->3'):" if is_pdna else "Sequence (5'->3'):")
        for widget in (self.seq2_label, self.seq2_entry, self.seq2_hint):
            if is_pdna:
                widget.grid()
            else:
                widget.grid_remove()
        if na_type == "Z-DNA":
            self.seq_entry.configure(state="disabled")
            self.custom_btn.configure(state="disabled")
            self.z_len_label.grid()
            self.z_len_entry.grid()
            self.z_len_entry.configure(state="normal")
            self.z_hint.grid()
            self.min_frame.grid_remove()
        else:
            self.seq_entry.configure(state="normal")
            # P-DNA parameters come from the 1PFI template and cannot be customized.
            self.custom_btn.configure(state="disabled" if is_pdna else "normal")
            self.z_len_entry.configure(state="disabled")
            self.z_len_label.grid_remove()
            self.z_len_entry.grid_remove()
            self.z_hint.grid_remove()
            self.min_frame.grid()
            self.minimize_var.set(self.minimize_by_type.get(na_type, False))
            self.regularize_phosphates_var.set(self.regularize_phosphates_by_type.get(na_type, False))
            self._refresh_minimization_hint()

        self._active_na_type = na_type
        self._refresh_param_values_display()
        self._update_sequence_length()
        self._refresh_info_text()

    def check_dssr_on_startup(self) -> None:
        info = check_dssr_installation()
        if info.get("installed"):
            status = f"x3dna-dssr: FOUND — {info.get('executable')}"
        else:
            status = "x3dna-dssr: NOT FOUND"
        self.dssr_status_var.set(status)

        report = [
            "=== x3dna-dssr startup check ===",
            f"Installed : {info.get('installed')}",
            f"Executable: {info.get('executable')}",
            f"Command   : {info.get('command')}",
            f"Returncode: {info.get('returncode')}",
            "Output:",
            str(info.get("output") or ""),
        ]
        self._append_log("\n".join(report))

    # ------------------------------------------------------------------
    # Build/generate
    # ------------------------------------------------------------------
    def _validate_params_for_minimization(self, run_phenix: bool) -> Optional[str]:
        if not run_phenix:
            return None
        params = self.params_var.get().strip()
        if not params:
            raise ValueError("Please specify a params file for phenix.geometry_minimization.")
        path = Path(params).expanduser()
        if not path.is_absolute():
            path = path.resolve()
        if not path.exists():
            raise ValueError(f"Params file not found: {path}")
        return str(path)

    def on_generate(self) -> None:
        na_type = self.na_type_var.get().strip()
        try:
            param_overrides = self._get_param_overrides()
            placement_values = self._get_placement_values()
            output_dir, tmp_dir = _ensure_output_dirs(self.output_dir_var.get())
            invrot_enabled = bool(self.invrot_enabled_var.get())
            invrot_operation = self.invrot_operation_var.get().strip() if invrot_enabled else ""
            if invrot_enabled:
                parse_operation(invrot_operation)
        except Exception as exc:
            messagebox.showerror("Input error", str(exc), parent=self)
            return

        user_name = sanitize_basename(self.name_var.get().strip())
        deleteH = bool(self.deleteH_var.get())
        run_phenix = bool(self.minimize_var.get()) if na_type != "Z-DNA" else False
        run_regularize_phosphates = (
            bool(self.regularize_phosphates_var.get()) if na_type != "Z-DNA" else False
        )

        try:
            params_for_min = self._validate_params_for_minimization(run_phenix)
        except Exception as exc:
            messagebox.showerror("phenix.geometry_minimization", str(exc), parent=self)
            return

        self.generate_btn.configure(state="disabled")
        self._set_log(
            f"=== {APP_NAME} {__version__} job started ===\n"
            f"Type: {na_type}\n"
            f"Output folder: {output_dir}\n"
            f"Intermediate folder: {tmp_dir}\n"
            f"phenix.geometry_minimization: {'ON' if run_phenix else 'OFF'}\n"
            f"Regularize phosphates: {'ON' if run_regularize_phosphates else 'OFF'}\n"
            f"Delete hydrogens: {deleteH}\n"
            f"L-form inv/rot mirror: {'ON (' + invrot_operation + ')' if invrot_enabled else 'OFF'}\n"
        )

        try:
            if na_type == "B-DNA":
                seq_in = self.seq_var.get().strip()
                if not seq_in:
                    raise ValueError("Please enter a DNA sequence using A/T/C/G.")
                result = build_bdna(
                    seq_in,
                    user_name,
                    params_for_min,
                    tmp_dir,
                    param_overrides=param_overrides,
                    deleteH=deleteH,
                    run_phenix=run_phenix,
                    run_regularize_phosphates=run_regularize_phosphates,
                )

            elif na_type == "A-DNA":
                seq_in = self.seq_var.get().strip()
                if not seq_in:
                    raise ValueError("Please enter a DNA sequence using A/T/C/G.")
                result = build_adna(
                    seq_in,
                    user_name,
                    tmp_dir,
                    param_overrides=param_overrides,
                    deleteH=deleteH,
                    run_phenix=run_phenix,
                    params_file=params_for_min,
                    run_regularize_phosphates=run_regularize_phosphates,
                )

            elif na_type == "A-RNA":
                seq_in = self.seq_var.get().strip()
                if not seq_in:
                    raise ValueError("Please enter an RNA sequence using A/U/C/G.")
                result = build_arna(
                    seq_in,
                    user_name,
                    tmp_dir,
                    param_overrides=param_overrides,
                    deleteH=deleteH,
                    run_phenix=run_phenix,
                    params_file=params_for_min,
                    run_regularize_phosphates=run_regularize_phosphates,
                )

            elif na_type == "Z-DNA":
                try:
                    length_value = _parse_int_expression(self.z_len_var.get().strip(), "Z-DNA helix length")
                except Exception as exc:
                    raise ValueError("Z-DNA helix length must be a positive even integer.") from exc
                result = build_zdna(length_value, user_name, tmp_dir, deleteH=deleteH)

            elif na_type == "P-DNA":
                seq_in = self.seq_var.get().strip()
                if not seq_in:
                    raise ValueError("Please enter the strand I DNA sequence using A/T/C/G.")
                result = build_pdna(
                    seq_in,
                    user_name,
                    tmp_dir,
                    seq_strand2=self.seq2_var.get().strip() or None,
                    deleteH=deleteH,
                    run_phenix=run_phenix,
                    params_file=params_for_min,
                    run_regularize_phosphates=run_regularize_phosphates,
                )

            else:
                raise ValueError(f"Unsupported type: {na_type}")

            aligned_pdb = Path(str(result["pdb_aligned"]))
            placement_input_pdb = aligned_pdb
            invrot_result = None
            mirror_suffix = ""
            if invrot_enabled:
                invrot_result = apply_inv_rot_to_pdb(aligned_pdb, instruction=invrot_operation)
                placement_input_pdb = Path(str(invrot_result["pdb_out"]))
                mirror_suffix = f"L_{invrot_result['label']}"

            placed_pdb = _final_placed_path(
                output_dir,
                str(result.get("base_name", "bnp_na_helix")),
                suffix=mirror_suffix,
            )
            place_result = place_after_Z(
                str(placement_input_pdb),
                str(placed_pdb),
                roll_deg=placement_values["roll"],
                phi_deg=placement_values["phi"],
                theta_deg=placement_values["theta"],
                tx=placement_values["x"],
                ty=placement_values["y"],
                tz=placement_values["z"],
                delta_z=placement_values["delta_z"],
            )
            remarks_log = _prepend_final_pdb_remarks(
                placed_pdb,
                na_type=na_type,
                l_form_enabled=invrot_enabled,
                invrot_result=invrot_result,
            )

            length_unit = "nt per strand" if na_type == "P-DNA" else "bp"
            log_parts = [
                f"=== {APP_NAME} {__version__} job summary ===",
                f"Final placed PDB: {placed_pdb}",
                f"Intermediate files: {tmp_dir}",
                f"Length: {result.get('length')} {length_unit}",
                "",
                str(result.get("log_text", "")),
            ]
            if param_overrides and na_type in NA_TYPES_WITH_TABLE:
                log_parts += [
                    "",
                    "=== GUI parameter overrides applied ===",
                    ", ".join(f"{key}={value:.4f}" for key, value in param_overrides.items()),
                ]
            if invrot_result:
                log_parts += ["", str(invrot_result.get("log_text", ""))]
            log_parts += ["", str(place_result.get("log_text", "")), "", remarks_log]
            self._set_log("\n".join(log_parts))
            messagebox.showinfo("bnp_na complete", f"Final placed PDB:\n{placed_pdb}", parent=self)

        except (PipelineError, PlacerError, InvRotError) as exc:
            self._set_log(getattr(exc, "log_text", "") or str(exc))
            messagebox.showerror("Job failed", str(exc), parent=self)
        except Exception as exc:
            self._set_log(f"Unexpected error:\n{type(exc).__name__}: {exc}")
            messagebox.showerror("Unexpected error", str(exc), parent=self)
        finally:
            self.generate_btn.configure(state="normal")


def main(argv: Optional[list[str]] = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if "-v" in args or "--version" in args:
        print(f"{APP_NAME} {__version__}")
        return
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
