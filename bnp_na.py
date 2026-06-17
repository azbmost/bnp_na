#!/usr/bin/env python3
"""bnp_na V13.6: Building and placing nucleic acid helices.

Top-level GUI/controller. All helper modules live in ./bnp_na_lib/.
"""
from __future__ import annotations

import ast
import math
import re
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Dict, Optional, Tuple

__version__ = "V13.6"
APP_NAME = "bnp_na"

APP_DIR = Path(__file__).resolve().parent
LIB_DIR = APP_DIR / "bnp_na_lib"
sys.path.insert(0, str(LIB_DIR))

from build_adna import build_adna  # noqa: E402
from build_arna import build_arna  # noqa: E402
from build_bz import build_bz_structure  # noqa: E402
from build_bdna import build_bdna  # noqa: E402
from build_common import (  # noqa: E402
    DEFAULT_PARAMS,
    PARAM_KEYS,
    PipelineError,
    check_dssr_installation,
    expand_sequence,
    sanitize_basename,
    sequence_alphabet,
)
from build_zdna import build_zdna  # noqa: E402
from angle_helical_axisV2_1 import launch_gui as launch_axis_angle_gui  # noqa: E402
from pdb_inv_rotV2 import InvRotError, apply_inv_rot_to_pdb, parse_operation  # noqa: E402
from xyz_bild import write_xyz_bild  # noqa: E402
from na_placer import PlacerError, place_after_Z  # noqa: E402


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
NA_TYPES_WITH_TABLE = ("B-DNA", "A-DNA", "A-RNA")
DEFAULT_MINIMIZE_BY_TYPE = {"B-DNA": True, "A-DNA": False, "A-RNA": False}
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
        self.geometry("1240x1080")
        self.minsize(1080, 920)
        self._style = ttk.Style(self)
        self._style.configure("Bold.TLabelframe.Label", font=("Helvetica", 11, "bold"))
        self._style.configure("Hint.TLabel", foreground="#666", font=("Helvetica", 9))
        self._style.configure("Status.TLabel", foreground="#555", font=("Helvetica", 9))
        self._style.configure("Compact.Treeview", rowheight=18, font=("Helvetica", 9))
        self._style.configure("Compact.Treeview.Heading", font=("Helvetica", 9, "bold"))

        self.param_values_by_type: Dict[str, Dict[str, str]] = {
            na_type: {key: "" for key in PARAM_KEYS} for na_type in NA_TYPES_WITH_TABLE
        }
        self.minimize_by_type: Dict[str, bool] = dict(DEFAULT_MINIMIZE_BY_TYPE)
        self._active_na_type: Optional[str] = None

        self.grid_columnconfigure(0, weight=0, minsize=260)
        self.grid_columnconfigure(1, weight=1, minsize=460)
        self.grid_columnconfigure(2, weight=0, minsize=120)
        self.grid_columnconfigure(3, weight=1, minsize=300)
        self.grid_rowconfigure(14, weight=1)

        pad = {"padx": 10, "pady": 3}
        frame_pad = {"padx": 12, "pady": 4}
        inner_y = 3

        title = ttk.Label(
            self,
            text=f"{APP_NAME} {__version__} — Module #1 of AZBMOST package: building and placing nucleic acid",
            font=("Helvetica", 14, "bold"),
        )
        title.grid(row=0, column=0, columnspan=4, sticky="w", padx=12, pady=(10, 4))

        self.dssr_status_var = tk.StringVar(value="x3dna-dssr: checking after GUI starts...")
        ttk.Label(self, textvariable=self.dssr_status_var, style="Status.TLabel").grid(
            row=1, column=0, columnspan=4, sticky="w", padx=12, pady=(0, 2)
        )

        ttk.Label(self, text="Sequence (5'->3'):").grid(
            row=2, column=0, sticky="e", **pad
        )
        self.seq_var = tk.StringVar(value="")
        self.seq_entry = ttk.Entry(self, textvariable=self.seq_var)
        self.seq_entry.grid(row=2, column=1, columnspan=2, sticky="we", **pad)
        self.length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(self, textvariable=self.length_var, style="Status.TLabel").grid(row=2, column=3, sticky="w", **pad)

        self.z_len_label = ttk.Label(self, text="Z-DNA helix length:")
        self.z_len_var = tk.StringVar(value="")
        self.z_len_entry = ttk.Entry(self, textvariable=self.z_len_var, width=22, state="disabled")
        self.z_hint = ttk.Label(
            self,
            text="Z-DNA fiber has a fixed sequence; the Sequence field and DSSR parameter table are ignored.",
            style="Hint.TLabel",
            wraplength=460,
        )
        self.z_len_label.grid(row=3, column=0, sticky="e", **pad)
        self.z_len_entry.grid(row=3, column=1, sticky="w", **pad)
        self.z_hint.grid(row=3, column=2, columnspan=2, sticky="w", **pad)
        self.z_len_label.grid_remove()
        self.z_len_entry.grid_remove()
        self.z_hint.grid_remove()

        ttk.Label(self, text="Helix name (optional):").grid(row=4, column=0, sticky="e", **pad)
        self.name_var = tk.StringVar(value="")
        self.name_entry = ttk.Entry(self, textvariable=self.name_var)
        self.name_entry.grid(row=4, column=1, columnspan=2, sticky="we", **pad)

        ttk.Label(self, text="Output folder:").grid(row=5, column=0, sticky="e", **pad)
        self.output_dir_var = tk.StringVar(value=_path_text(DEFAULT_OUTPUT_DIR))
        self.output_dir_entry = ttk.Entry(self, textvariable=self.output_dir_var)
        self.output_dir_entry.grid(row=5, column=1, sticky="we", **pad)
        ttk.Button(self, text="Browse", command=self.browse_output_dir).grid(row=5, column=2, sticky="w", **pad)

        type_frame = ttk.LabelFrame(self, text="Nucleic acid type", style="Bold.TLabelframe")
        type_frame.grid(row=6, column=0, columnspan=4, sticky="we", **frame_pad)
        self.na_type_var = tk.StringVar(value="B-DNA")
        for idx, label in enumerate(["B-DNA", "A-DNA", "A-RNA", "Z-DNA"]):
            ttk.Radiobutton(
                type_frame,
                text=label,
                value=label,
                variable=self.na_type_var,
                command=self.on_type_changed,
            ).grid(row=0, column=idx, padx=8, pady=inner_y, sticky="w")
        ttk.Button(type_frame, text="B-Z builder", command=self.open_bz_builder_dialog).grid(
            row=0, column=4, sticky="w", padx=(18, 8), pady=inner_y
        )

        self.param_frame = ttk.LabelFrame(
            self,
            text="Current DSSR helical parameters",
            style="Bold.TLabelframe",
        )
        self.param_frame.grid(row=7, column=0, columnspan=4, sticky="we", **frame_pad)
        self.param_frame.grid_columnconfigure(1, weight=1)

        self.custom_btn = ttk.Button(
            self.param_frame,
            text="Customize DSSR parameters",
            command=self.open_param_dialog,
        )
        self.custom_btn.grid(row=0, column=0, sticky="w", padx=8, pady=(4, 3))
        self.param_status_var = tk.StringVar(value="")
        ttk.Label(
            self.param_frame,
            textvariable=self.param_status_var,
            style="Status.TLabel",
            justify="left",
            wraplength=820,
        ).grid(row=0, column=1, sticky="w", padx=8, pady=(4, 3))

        self.param_table = ttk.Treeview(
            self.param_frame,
            columns=("row_label", *PARAM_KEYS),
            show="headings",
            height=3,
            selectmode="none",
            style="Compact.Treeview",
        )
        for column, heading, width, anchor in (("row_label", "", 72, "w"),):
            self.param_table.heading(column, text=heading)
            self.param_table.column(column, width=width, anchor=anchor, stretch=False)
        for key in PARAM_KEYS:
            self.param_table.heading(key, text=PARAM_LABELS.get(key, key))
            self.param_table.column(key, width=80, anchor="center", stretch=True)
        self.param_table.grid(row=1, column=0, columnspan=2, sticky="we", padx=8, pady=(0, 4))

        self.min_frame = ttk.LabelFrame(self, text="phenix.geometry_minimization", style="Bold.TLabelframe")
        self.min_frame.grid(row=8, column=0, columnspan=4, sticky="we", **frame_pad)
        self.min_frame.grid_columnconfigure(1, weight=1)
        self.minimize_var = tk.BooleanVar(value=DEFAULT_MINIMIZE_BY_TYPE["B-DNA"])
        self.min_check = ttk.Checkbutton(
            self.min_frame,
            text="Run phenix.geometry_minimization",
            variable=self.minimize_var,
            command=self._on_minimize_toggled,
        )
        self.min_check.grid(row=0, column=0, sticky="w", padx=8, pady=inner_y)
        ttk.Label(self.min_frame, text="Params file (.eff / .params):").grid(row=0, column=1, sticky="e", padx=8, pady=inner_y)
        self.params_var = tk.StringVar(value=_path_text(DEFAULT_PARAMS_FILE))
        self.params_entry = ttk.Entry(self.min_frame, textvariable=self.params_var)
        self.params_entry.grid(row=0, column=2, sticky="we", padx=8, pady=inner_y)
        self.min_frame.grid_columnconfigure(2, weight=1)
        self.params_browse_btn = ttk.Button(self.min_frame, text="Browse", command=self.browse_params)
        self.params_browse_btn.grid(row=0, column=3, sticky="w", padx=8, pady=inner_y)
        self.min_hint_var = tk.StringVar(value="")
        ttk.Label(self.min_frame, textvariable=self.min_hint_var, style="Hint.TLabel", wraplength=960).grid(
            row=1, column=0, columnspan=4, sticky="w", padx=8, pady=(0, 3)
        )

        options_row = ttk.Frame(self)
        options_row.grid(row=9, column=0, columnspan=4, sticky="w", padx=12, pady=(0, 3))
        self.deleteH_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            options_row,
            text="Delete hydrogens from generated PDB",
            variable=self.deleteH_var,
        ).pack(side="left", padx=8)

        self.chirality_frame = ttk.LabelFrame(
            self,
            text="Mirror-image L-form chirality (L-DNA)",
            style="Bold.TLabelframe",
        )
        self.chirality_frame.grid(row=10, column=0, columnspan=4, sticky="we", **frame_pad)
        self.chirality_frame.grid_columnconfigure(4, weight=1)
        self.invrot_enabled_var = tk.BooleanVar(value=False)
        self.invrot_operation_var = tk.StringVar(value="oyz")
        self.invrot_check = ttk.Checkbutton(
            self.chirality_frame,
            text="Apply inv/rot after align-to-Z",
            variable=self.invrot_enabled_var,
            command=self._on_invrot_toggled,
        )
        self.invrot_check.grid(row=0, column=0, sticky="w", padx=8, pady=inner_y)
        ttk.Label(self.chirality_frame, text="Operation:").grid(row=0, column=1, sticky="e", padx=(18, 6), pady=inner_y)
        self.invrot_operation_combo = ttk.Combobox(
            self.chirality_frame,
            textvariable=self.invrot_operation_var,
            values=INV_ROT_OPERATIONS,
            width=8,
            state="disabled",
        )
        self.invrot_operation_combo.grid(row=0, column=2, sticky="w", padx=4, pady=inner_y)
        self.invrot_operation_combo.bind("<<ComboboxSelected>>", lambda *_args: self._refresh_info_text())
        ttk.Button(self.chirality_frame, text="Help", command=self.open_invrot_help).grid(
            row=0, column=3, sticky="w", padx=8, pady=inner_y
        )
        self.invrot_hint_var = tk.StringVar(value="")
        ttk.Label(self.chirality_frame, textvariable=self.invrot_hint_var, style="Hint.TLabel", wraplength=560).grid(
            row=0, column=4, sticky="w", padx=8, pady=inner_y
        )
        self._refresh_invrot_state()

        place = ttk.LabelFrame(self, text="Placement / Orientation", style="Bold.TLabelframe")
        place.grid(row=11, column=0, columnspan=4, sticky="we", **frame_pad)
        for col, minsize in {0: 80, 1: 120, 2: 80, 3: 120, 4: 80, 5: 120}.items():
            place.grid_columnconfigure(col, minsize=minsize)
        place.grid_columnconfigure(7, weight=1)

        self.x_var = tk.StringVar(value="0")
        self.y_var = tk.StringVar(value="0")
        self.z_var = tk.StringVar(value="0")
        self.roll_var = tk.StringVar(value="0")
        self.phi_var = tk.StringVar(value="0")
        self.theta_var = tk.StringVar(value="0")
        self.delta_z_var = tk.StringVar(value="0")

        ttk.Label(place, text="delta_z (Å)").grid(row=0, column=0, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.delta_z_var, width=10).grid(row=0, column=1, sticky="w", padx=2, pady=inner_y)
        ttk.Label(place, text="delta_z should be 0 for most cases.", style="Status.TLabel").grid(
            row=0, column=2, columnspan=4, sticky="w", padx=(20, 6), pady=inner_y
        )

        ttk.Label(place, text="x (Å)").grid(row=1, column=0, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.x_var, width=10).grid(row=1, column=1, sticky="w", padx=2, pady=inner_y)
        ttk.Label(place, text="y (Å)").grid(row=1, column=2, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.y_var, width=10).grid(row=1, column=3, sticky="w", padx=2, pady=inner_y)
        ttk.Label(place, text="z (Å)").grid(row=1, column=4, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.z_var, width=10).grid(row=1, column=5, sticky="w", padx=2, pady=inner_y)

        ttk.Label(place, text="roll (°)").grid(row=2, column=0, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.roll_var, width=10).grid(row=2, column=1, sticky="w", padx=2, pady=inner_y)
        ttk.Label(place, text="phi (°)").grid(row=2, column=2, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.phi_var, width=10).grid(row=2, column=3, sticky="w", padx=2, pady=inner_y)
        ttk.Label(place, text="theta (°)").grid(row=2, column=4, sticky="e", padx=6, pady=inner_y)
        ttk.Entry(place, textvariable=self.theta_var, width=10).grid(row=2, column=5, sticky="w", padx=2, pady=inner_y)

        ttk.Label(
            place,
            text="For GIDEON: roll = roll at GIDEON - 111.25",
            style="Status.TLabel",
            justify="left",
            wraplength=420,
        ).grid(row=3, column=0, columnspan=4, sticky="w", padx=(20, 6), pady=(0, 3))

        self.info_var = tk.StringVar(value="")
        ttk.Label(self, textvariable=self.info_var, wraplength=1120, style="Status.TLabel", justify="left").grid(
            row=12, column=0, columnspan=4, sticky="we", padx=12, pady=(0, 3)
        )

        tools = ttk.LabelFrame(self, text="Analysis tools", style="Bold.TLabelframe")
        tools.grid(row=13, column=0, columnspan=4, sticky="we", **frame_pad)
        tools.grid_columnconfigure(2, weight=1)
        ttk.Button(tools, text="Open helical-axis angle tool", command=self.open_axis_angle_tool).grid(
            row=0, column=0, sticky="w", padx=8, pady=inner_y
        )
        ttk.Button(tools, text="Write XYZ axes BILD", command=self.open_xyz_bild_dialog).grid(
            row=0, column=1, sticky="w", padx=8, pady=inner_y
        )
        ttk.Label(
            tools,
            text="Measure around-axis angles or create coordinate-axis .bild helpers for Chimera/ChimeraX.",
            style="Hint.TLabel",
            wraplength=650,
        ).grid(row=0, column=2, sticky="w", padx=8, pady=inner_y)

        log_frame = ttk.LabelFrame(self, text="Log output", style="Bold.TLabelframe")
        log_frame.grid(row=14, column=0, columnspan=4, sticky="nsew", padx=12, pady=5)
        log_frame.grid_columnconfigure(0, weight=1)
        log_frame.grid_rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, wrap="none", height=12)
        try:
            self.log_text.configure(font=("Menlo", 10))
        except Exception:
            self.log_text.configure(font=("Courier", 10))
        yscroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        xscroll = ttk.Scrollbar(log_frame, orient="horizontal", command=self.log_text.xview)
        self.log_text.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        self.log_text.configure(state="disabled")

        btn_row = ttk.Frame(self)
        btn_row.grid(row=15, column=0, columnspan=4, sticky="e", padx=12, pady=(4, 8))
        ttk.Button(btn_row, text="Quit", command=self.destroy).pack(side="left", padx=8)
        self.generate_btn = ttk.Button(btn_row, text="Generate", command=self.on_generate)
        self.generate_btn.pack(side="left", padx=8)

        self.seq_var.trace_add("write", lambda *_args: self._update_sequence_length())
        self.z_len_var.trace_add("write", lambda *_args: self._update_sequence_length())
        self.output_dir_var.trace_add("write", lambda *_args: self._refresh_info_text())
        self.invrot_operation_var.trace_add("write", lambda *_args: self._refresh_invrot_state())

        self.on_type_changed()
        self._refresh_info_text()
        self._set_log(f"=== {APP_NAME} {__version__} startup ===\nApp folder: {APP_DIR}\nHelper folder: {LIB_DIR}\n")
        self.after(100, self.check_dssr_on_startup)

    def _set_optional_window_icon(self) -> None:
        if not DEFAULT_ICON_FILE.exists():
            return
        try:
            self._window_icon = tk.PhotoImage(file=str(DEFAULT_ICON_FILE))
            self.iconphoto(True, self._window_icon)
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

    def _refresh_minimization_hint(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type == "B-DNA":
            self.min_hint_var.set("Default for B-DNA is to run phenix.geometry_minimization, matching the previous B-DNA workflow.")
        elif na_type in ("A-DNA", "A-RNA"):
            self.min_hint_var.set("Default for A-DNA/A-RNA is to skip minimization; check the box to run phenix.geometry_minimization.")
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
        win.title("L-form inv/rot help")
        win.geometry("760x520+220+140")
        win.minsize(680, 460)
        win.transient(self)

        text = tk.Text(win, wrap="word", height=22)
        try:
            text.configure(font=("Helvetica", 10))
        except Exception:
            pass
        text.grid(row=0, column=0, sticky="nsew", padx=12, pady=12)
        win.grid_columnconfigure(0, weight=1)
        win.grid_rowconfigure(0, weight=1)

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
        text.configure(state="disabled")

        ttk.Button(win, text="Close", command=win.destroy).grid(row=1, column=0, sticky="e", padx=12, pady=(0, 12))

    def open_axis_angle_tool(self) -> None:
        try:
            launch_axis_angle_gui(parent=self)
        except Exception as exc:
            messagebox.showerror("Helical-axis angle tool", str(exc), parent=self)

    def open_xyz_bild_dialog(self) -> None:
        win = tk.Toplevel(self)
        win.title("Write XYZ axes BILD")
        win.geometry("620x260+240+180")
        win.minsize(560, 240)
        win.transient(self)
        win.grid_columnconfigure(1, weight=1)

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
            messagebox.showinfo("XYZ axes BILD", f"Wrote:\n{written}", parent=win)
            win.destroy()

        pad = {"padx": 10, "pady": 4}
        ttk.Label(win, text="Output .bild:").grid(row=0, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=out_var).grid(row=0, column=1, sticky="we", **pad)
        ttk.Button(win, text="Browse", command=browse_out).grid(row=0, column=2, sticky="w", **pad)

        ttk.Label(win, text="Origin x y z (Å):").grid(row=1, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=origin_var, width=18).grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(win, text="Arrow length (Å):").grid(row=2, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=length_var, width=12).grid(row=2, column=1, sticky="w", **pad)

        ttk.Label(win, text="Arrow width (Å):").grid(row=3, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=width_var, width=12).grid(row=3, column=1, sticky="w", **pad)

        ttk.Label(win, text="Origin sphere radius (Å):").grid(row=4, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=sphere_var, width=12).grid(row=4, column=1, sticky="w", **pad)

        ttk.Label(
            win,
            text="Colors: X red, Y yellow, Z blue. Arrow head radius is 2.5 x arrow width.",
            style="Hint.TLabel",
            wraplength=460,
        ).grid(row=5, column=0, columnspan=3, sticky="w", padx=10, pady=(2, 6))

        buttons = ttk.Frame(win)
        buttons.grid(row=6, column=0, columnspan=3, sticky="e", padx=10, pady=(3, 8))
        ttk.Button(buttons, text="Cancel", command=win.destroy).pack(side="left", padx=6)
        ttk.Button(buttons, text="Write", command=write_file).pack(side="left", padx=6)

    def open_bz_builder_dialog(self) -> None:
        win = tk.Toplevel(self)
        win.title("B-Z structure builder")
        win.geometry("900x720+180+100")
        win.minsize(760, 620)
        win.transient(self)
        win.grid_columnconfigure(1, weight=1)
        win.grid_columnconfigure(2, weight=1)
        win.grid_rowconfigure(1, weight=1)
        win.grid_rowconfigure(9, weight=1)

        try:
            out_dir = Path(self.output_dir_var.get()).expanduser()
            if not out_dir.is_absolute():
                out_dir = out_dir.resolve()
        except Exception:
            out_dir = DEFAULT_OUTPUT_DIR

        files_text = scrolledtext.ScrolledText(win, width=90, height=7, wrap="none")
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

        pad = {"padx": 10, "pady": 4}
        ttk.Label(
            win,
            text="Input PDB files must alternate by order: B1, Z1, B2, Z2, ...",
            style="Hint.TLabel",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=10, pady=(10, 2))
        files_text.grid(row=1, column=0, columnspan=3, sticky="nsew", padx=10, pady=(0, 4))

        file_buttons = ttk.Frame(win)
        file_buttons.grid(row=2, column=0, columnspan=3, sticky="w", padx=10, pady=(0, 6))
        ttk.Button(file_buttons, text="Add input PDBs", command=add_files).pack(side="left", padx=(0, 8))
        ttk.Button(file_buttons, text="Clear", command=clear_files).pack(side="left")

        ttk.Label(win, text="Output PDB:").grid(row=3, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=out_var).grid(row=3, column=1, sticky="we", **pad)
        ttk.Button(win, text="Browse", command=browse_out).grid(row=3, column=2, sticky="w", **pad)

        ttk.Label(win, text="Axis correction:").grid(row=4, column=0, sticky="e", **pad)
        ttk.Combobox(
            win,
            textvariable=axis_mode_var,
            values=("codirectional", "collinear", "none"),
            state="readonly",
            width=16,
        ).grid(row=4, column=1, sticky="w", **pad)
        ttk.Label(
            win,
            text="Default codirectional rotates helix axes parallel without lateral axis-line shift.",
            style="Hint.TLabel",
            wraplength=420,
        ).grid(row=4, column=2, sticky="w", padx=10, pady=4)

        ttk.Label(win, text="Axis source:").grid(row=5, column=0, sticky="e", **pad)
        ttk.Combobox(
            win,
            textvariable=axis_source_var,
            values=("auto", "dssr", "pca"),
            state="readonly",
            width=16,
        ).grid(row=5, column=1, sticky="w", **pad)
        ttk.Label(
            win,
            text="auto tries DSSR --more through align2z.py, then falls back to C1' PCA.",
            style="Hint.TLabel",
            wraplength=420,
        ).grid(row=5, column=2, sticky="w", padx=10, pady=4)

        ttk.Checkbutton(
            win,
            text=(
                "Auto-trim terminal Z-DNA bp if needed. For DSSR/bnp_na Z inputs, "
                "prepare Z-DNA 2 bp longer than the target final Z segment."
            ),
            variable=auto_trim_var,
        ).grid(row=6, column=0, columnspan=3, sticky="w", padx=10, pady=(2, 4))

        ttk.Label(
            win,
            text=(
                "The builder inserts B-Z junction cores, trims overlapping base pairs, and writes both a final ligated "
                "PDB plus a raw aligned PDB with _raw before the file extension."
            ),
            style="Hint.TLabel",
            wraplength=820,
        ).grid(row=7, column=0, columnspan=3, sticky="w", padx=10, pady=(0, 4))

        buttons = ttk.Frame(win)
        buttons.grid(row=8, column=0, columnspan=3, sticky="e", padx=10, pady=(2, 4))
        ttk.Button(buttons, text="Close", command=win.destroy).pack(side="left", padx=6)
        run_btn = ttk.Button(buttons, text="Build B-Z structure")
        run_btn.pack(side="left", padx=6)

        bz_log = scrolledtext.ScrolledText(win, width=90, height=12, wrap="none")
        try:
            bz_log.configure(font=("Menlo", 10))
        except Exception:
            bz_log.configure(font=("Courier", 10))
        bz_log.grid(row=9, column=0, columnspan=3, sticky="nsew", padx=10, pady=(0, 10))

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

    def _refresh_info_text(self) -> None:
        out = self.output_dir_var.get().strip() or "<not selected>"
        tmp = str(Path(out).expanduser() / "tmp_file") if out != "<not selected>" else "<not selected>"
        pipeline = (
            "Pipeline: build by selected type -> normalize names -> optional phenix.geometry_minimization "
            "for B-DNA/A-DNA/A-RNA -> DSSR --more axis extraction -> align to +Z"
        )
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
        if defaults is None or store is None:
            messagebox.showinfo("Z-DNA", "Z-DNA uses DSSR fiber; the 12-parameter table is not used.", parent=self)
            return

        win = tk.Toplevel(self)
        win.title(f"Customize DSSR parameters — {na_type}")
        win.geometry("960x500+180+120")
        win.minsize(860, 430)
        win.transient(self)
        win.grab_set()

        ttk.Label(
            win,
            text=(
                "Leave a field blank to use the default value shown beside it. "
                "Previously saved custom values are kept. Values are written to four digits after the decimal."
            ),
            wraplength=900,
            style="Hint.TLabel",
            justify="left",
        ).grid(row=0, column=0, columnspan=6, sticky="w", padx=14, pady=(10, 6))

        local_vars: Dict[str, tk.StringVar] = {}
        for idx, (key, default) in enumerate(zip(PARAM_KEYS, defaults)):
            col_group = 0 if idx < 6 else 3
            row = idx + 1 if idx < 6 else idx - 5
            ttk.Label(win, text=PARAM_LABELS.get(key, key) + ":").grid(
                row=row, column=col_group, sticky="e", padx=10, pady=4
            )
            initial = store.get(key, "").strip()
            local_vars[key] = tk.StringVar(value=initial)
            ttk.Entry(win, textvariable=local_vars[key], width=18).grid(
                row=row, column=col_group + 1, sticky="w", padx=10, pady=4
            )
            ttk.Label(win, text=f"default {_format_number(default)} {PARAM_UNITS.get(key, '')}", style="Hint.TLabel").grid(
                row=row, column=col_group + 2, sticky="w", padx=(0, 12), pady=4
            )

        btns = ttk.Frame(win)
        btns.grid(row=8, column=0, columnspan=6, sticky="e", padx=12, pady=8)

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

        ttk.Button(btns, text="Clear fields", command=clear_all).pack(side="left", padx=6)
        ttk.Button(btns, text="Restore defaults", command=restore_defaults).pack(side="left", padx=6)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="left", padx=6)
        ttk.Button(btns, text="Save", command=save_close).pack(side="left", padx=6)

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

        na_type = self.na_type_var.get().strip()
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
            self.custom_btn.configure(state="normal")
            self.z_len_entry.configure(state="disabled")
            self.z_len_label.grid_remove()
            self.z_len_entry.grid_remove()
            self.z_hint.grid_remove()
            self.min_frame.grid()
            self.minimize_var.set(self.minimize_by_type.get(na_type, False))
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
                )

            elif na_type == "Z-DNA":
                try:
                    length_value = _parse_int_expression(self.z_len_var.get().strip(), "Z-DNA helix length")
                except Exception as exc:
                    raise ValueError("Z-DNA helix length must be a positive even integer.") from exc
                result = build_zdna(length_value, user_name, tmp_dir, deleteH=deleteH)

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

            log_parts = [
                f"=== {APP_NAME} {__version__} job summary ===",
                f"Final placed PDB: {placed_pdb}",
                f"Intermediate files: {tmp_dir}",
                f"Length: {result.get('length')} bp",
                "",
                str(result.get("log_text", "")),
            ]
            if param_overrides and na_type != "Z-DNA":
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
