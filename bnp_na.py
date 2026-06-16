#!/usr/bin/env python3
"""bnp_na V13.4: Building and placing nucleic acid helices.

Top-level GUI/controller. All helper modules live in ./bnp_na_lib/.
"""
from __future__ import annotations

import re
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Dict, Optional, Tuple

__version__ = "V13.4"
APP_NAME = "bnp_na"

APP_DIR = Path(__file__).resolve().parent
LIB_DIR = APP_DIR / "bnp_na_lib"
sys.path.insert(0, str(LIB_DIR))

from build_adna import build_adna  # noqa: E402
from build_arna import build_arna  # noqa: E402
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
from angle_helical_axisV2 import launch_gui as launch_axis_angle_gui  # noqa: E402
from pdb_inv_rotV2 import InvRotError, apply_inv_rot_to_pdb, parse_operation  # noqa: E402
from xyz_bild import write_xyz_bild  # noqa: E402
from na_placer import PlacerError, place_after_Z  # noqa: E402


PARAM_LABELS = {
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

DEFAULT_PARAMS_FILE = LIB_DIR / "min_P_C5.params"
DEFAULT_OUTPUT_DIR = APP_DIR / "output"
DEFAULT_ICON_FILE = APP_DIR / "assets" / "bnp_na_icon.png"
NA_TYPES_WITH_TABLE = ("B-DNA", "A-DNA", "A-RNA")
DEFAULT_MINIMIZE_BY_TYPE = {"B-DNA": True, "A-DNA": False, "A-RNA": False}
INV_ROT_OPERATIONS = ("oyz", "oxz", "oxy", "i", "ix", "iy", "iz", "ixy", "ixz", "iyz", "ixyz")


def _path_text(path: Path) -> str:
    try:
        return str(path.resolve())
    except Exception:
        return str(path)


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
        super().__init__()
        self.title(f"{APP_NAME} {__version__} - Build and Place Nucleic Acid")
        self._set_optional_window_icon()
        self.geometry("1240x1080")
        self.minsize(1080, 920)

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

        pad = {"padx": 12, "pady": 6}

        title = ttk.Label(
            self,
            text=f"{APP_NAME} {__version__} — Building and placing nucleic acid",
            font=("Helvetica", 14, "bold"),
        )
        title.grid(row=0, column=0, columnspan=4, sticky="w", padx=12, pady=(12, 6))

        self.dssr_status_var = tk.StringVar(value="x3dna-dssr: checking after GUI starts...")
        ttk.Label(self, textvariable=self.dssr_status_var, foreground="#444").grid(
            row=1, column=0, columnspan=4, sticky="w", padx=12, pady=(0, 4)
        )

        ttk.Label(self, text="Sequence (5'->3'):", font=("Helvetica", 11, "bold")).grid(
            row=2, column=0, sticky="e", **pad
        )
        self.seq_var = tk.StringVar(value="")
        self.seq_entry = ttk.Entry(self, textvariable=self.seq_var)
        self.seq_entry.grid(row=2, column=1, columnspan=2, sticky="we", **pad)
        self.length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(self, textvariable=self.length_var, foreground="#555").grid(row=2, column=3, sticky="w", **pad)

        self.z_len_label = ttk.Label(self, text="Z-DNA helix length:")
        self.z_len_var = tk.StringVar(value="")
        self.z_len_entry = ttk.Entry(self, textvariable=self.z_len_var, width=22, state="disabled")
        self.z_hint = ttk.Label(
            self,
            text="Z-DNA fiber has a fixed sequence; the Sequence field and DSSR parameter table are ignored.",
            foreground="#666",
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

        type_frame = ttk.LabelFrame(self, text="Nucleic acid type")
        type_frame.grid(row=6, column=0, columnspan=4, sticky="we", padx=12, pady=8)
        self.na_type_var = tk.StringVar(value="B-DNA")
        for idx, label in enumerate(["B-DNA", "A-DNA", "A-RNA", "Z-DNA"]):
            ttk.Radiobutton(
                type_frame,
                text=label,
                value=label,
                variable=self.na_type_var,
                command=self.on_type_changed,
            ).grid(row=0, column=idx, padx=8, pady=4, sticky="w")

        self.custom_btn = ttk.Button(
            self,
            text="Customize DSSR parameters",
            command=self.open_param_dialog,
        )
        self.custom_btn.grid(row=7, column=0, sticky="e", **pad)
        self.param_values_var = tk.StringVar(value="")
        self.param_values_lbl = ttk.Label(
            self,
            textvariable=self.param_values_var,
            foreground="#555",
            justify="left",
            wraplength=860,
        )
        self.param_values_lbl.grid(row=7, column=1, columnspan=3, sticky="w", **pad)

        self.min_frame = ttk.LabelFrame(self, text="phenix.geometry_minimization")
        self.min_frame.grid(row=8, column=0, columnspan=4, sticky="we", padx=12, pady=6)
        self.min_frame.grid_columnconfigure(1, weight=1)
        self.minimize_var = tk.BooleanVar(value=DEFAULT_MINIMIZE_BY_TYPE["B-DNA"])
        self.min_check = ttk.Checkbutton(
            self.min_frame,
            text="Run phenix.geometry_minimization",
            variable=self.minimize_var,
            command=self._on_minimize_toggled,
        )
        self.min_check.grid(row=0, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(self.min_frame, text="Params file (.eff / .params):").grid(row=0, column=1, sticky="e", padx=8, pady=6)
        self.params_var = tk.StringVar(value=_path_text(DEFAULT_PARAMS_FILE))
        self.params_entry = ttk.Entry(self.min_frame, textvariable=self.params_var)
        self.params_entry.grid(row=0, column=2, sticky="we", padx=8, pady=6)
        self.min_frame.grid_columnconfigure(2, weight=1)
        self.params_browse_btn = ttk.Button(self.min_frame, text="Browse", command=self.browse_params)
        self.params_browse_btn.grid(row=0, column=3, sticky="w", padx=8, pady=6)
        self.min_hint_var = tk.StringVar(value="")
        ttk.Label(self.min_frame, textvariable=self.min_hint_var, foreground="#666", wraplength=960).grid(
            row=1, column=0, columnspan=4, sticky="w", padx=8, pady=(0, 6)
        )

        options_row = ttk.Frame(self)
        options_row.grid(row=9, column=0, columnspan=4, sticky="w", padx=12, pady=(0, 6))
        self.deleteH_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            options_row,
            text="Delete hydrogens from generated PDB",
            variable=self.deleteH_var,
        ).pack(side="left", padx=8)

        self.chirality_frame = ttk.LabelFrame(self, text="Mirror-image L-form chirality")
        self.chirality_frame.grid(row=10, column=0, columnspan=4, sticky="we", padx=12, pady=6)
        self.chirality_frame.grid_columnconfigure(4, weight=1)
        self.invrot_enabled_var = tk.BooleanVar(value=False)
        self.invrot_operation_var = tk.StringVar(value="oyz")
        self.invrot_check = ttk.Checkbutton(
            self.chirality_frame,
            text="Apply inv/rot after align-to-Z",
            variable=self.invrot_enabled_var,
            command=self._on_invrot_toggled,
        )
        self.invrot_check.grid(row=0, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(self.chirality_frame, text="Operation:").grid(row=0, column=1, sticky="e", padx=(18, 6), pady=6)
        self.invrot_operation_combo = ttk.Combobox(
            self.chirality_frame,
            textvariable=self.invrot_operation_var,
            values=INV_ROT_OPERATIONS,
            width=8,
            state="disabled",
        )
        self.invrot_operation_combo.grid(row=0, column=2, sticky="w", padx=4, pady=6)
        self.invrot_operation_combo.bind("<<ComboboxSelected>>", lambda *_args: self._refresh_info_text())
        ttk.Button(self.chirality_frame, text="Help", command=self.open_invrot_help).grid(
            row=0, column=3, sticky="w", padx=8, pady=6
        )
        self.invrot_hint_var = tk.StringVar(value="")
        ttk.Label(self.chirality_frame, textvariable=self.invrot_hint_var, foreground="#666", wraplength=560).grid(
            row=0, column=4, sticky="w", padx=8, pady=6
        )
        self._refresh_invrot_state()

        place = ttk.LabelFrame(self, text="Placement / Orientation")
        place.grid(row=11, column=0, columnspan=4, sticky="we", padx=12, pady=8)
        for col, minsize in {0: 80, 1: 120, 2: 80, 3: 120, 4: 80, 5: 120}.items():
            place.grid_columnconfigure(col, minsize=minsize)
        place.grid_columnconfigure(7, weight=1)

        self.x_var = tk.DoubleVar(value=0.0)
        self.y_var = tk.DoubleVar(value=0.0)
        self.z_var = tk.DoubleVar(value=0.0)
        self.roll_var = tk.DoubleVar(value=0.0)
        self.phi_var = tk.DoubleVar(value=0.0)
        self.theta_var = tk.DoubleVar(value=0.0)
        self.delta_z_var = tk.DoubleVar(value=0.0)

        ttk.Label(place, text="delta_z (A)").grid(row=0, column=0, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.delta_z_var, width=10).grid(row=0, column=1, sticky="w", padx=2, pady=5)
        ttk.Label(place, text="delta_z should be 0 for most cases.", foreground="#555").grid(
            row=0, column=2, columnspan=4, sticky="w", padx=(20, 6), pady=5
        )

        ttk.Label(place, text="x (A)").grid(row=1, column=0, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.x_var, width=10).grid(row=1, column=1, sticky="w", padx=2, pady=5)
        ttk.Label(place, text="y (A)").grid(row=1, column=2, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.y_var, width=10).grid(row=1, column=3, sticky="w", padx=2, pady=5)
        ttk.Label(place, text="z (A)").grid(row=1, column=4, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.z_var, width=10).grid(row=1, column=5, sticky="w", padx=2, pady=5)

        ttk.Label(place, text="roll (deg)").grid(row=2, column=0, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.roll_var, width=10).grid(row=2, column=1, sticky="w", padx=2, pady=5)
        ttk.Label(place, text="phi (deg)").grid(row=2, column=2, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.phi_var, width=10).grid(row=2, column=3, sticky="w", padx=2, pady=5)
        ttk.Label(place, text="theta (deg)").grid(row=2, column=4, sticky="e", padx=6, pady=5)
        ttk.Entry(place, textvariable=self.theta_var, width=10).grid(row=2, column=5, sticky="w", padx=2, pady=5)

        ttk.Label(
            place,
            text="For GIDEON: roll = roll at GIDEON - 111.25",
            foreground="#555",
            justify="left",
            wraplength=420,
        ).grid(row=3, column=0, columnspan=4, sticky="w", padx=(20, 6), pady=(0, 6))

        self.info_var = tk.StringVar(value="")
        ttk.Label(self, textvariable=self.info_var, wraplength=1120, foreground="#444", justify="left").grid(
            row=12, column=0, columnspan=4, sticky="we", padx=12, pady=(0, 6)
        )

        tools = ttk.LabelFrame(self, text="Analysis tools")
        tools.grid(row=13, column=0, columnspan=4, sticky="we", padx=12, pady=6)
        tools.grid_columnconfigure(1, weight=1)
        ttk.Button(tools, text="Open helical-axis angle tool", command=self.open_axis_angle_tool).grid(
            row=0, column=0, sticky="w", padx=8, pady=6
        )
        ttk.Label(
            tools,
            text="Measure the around-axis angle between two atom or XYZ points and write a Chimera/ChimeraX .bild file.",
            foreground="#666",
            wraplength=850,
        ).grid(row=0, column=1, sticky="w", padx=8, pady=6)
        ttk.Button(tools, text="Write XYZ axes BILD", command=self.open_xyz_bild_dialog).grid(
            row=1, column=0, sticky="w", padx=8, pady=6
        )
        ttk.Label(
            tools,
            text="Create a red/yellow/blue coordinate-axis .bild helper with configurable arrow length and width.",
            foreground="#666",
            wraplength=850,
        ).grid(row=1, column=1, sticky="w", padx=8, pady=6)

        log_frame = ttk.LabelFrame(self, text="Log output")
        log_frame.grid(row=14, column=0, columnspan=4, sticky="nsew", padx=12, pady=8)
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
        btn_row.grid(row=15, column=0, columnspan=4, sticky="e", padx=12, pady=(6, 12))
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
            text.configure(font=("Helvetica", 11))
        except Exception:
            pass
        text.grid(row=0, column=0, sticky="nsew", padx=14, pady=14)
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

        ttk.Button(win, text="Close", command=win.destroy).grid(row=1, column=0, sticky="e", padx=14, pady=(0, 14))

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
            parts = text.replace(",", " ").split()
            if len(parts) != 3:
                raise ValueError("Origin must have three numbers, for example: 0 0 0")
            return float(parts[0]), float(parts[1]), float(parts[2])

        def write_file() -> None:
            try:
                out_path = Path(out_var.get()).expanduser()
                if not str(out_path).strip():
                    raise ValueError("Please choose an output .bild file.")
                origin = parse_origin(origin_var.get())
                length = float(length_var.get())
                width = float(width_var.get())
                sphere_radius = float(sphere_var.get())
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

        pad = {"padx": 10, "pady": 6}
        ttk.Label(win, text="Output .bild:").grid(row=0, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=out_var).grid(row=0, column=1, sticky="we", **pad)
        ttk.Button(win, text="Browse", command=browse_out).grid(row=0, column=2, sticky="w", **pad)

        ttk.Label(win, text="Origin x y z:").grid(row=1, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=origin_var, width=18).grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(win, text="Arrow length:").grid(row=2, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=length_var, width=12).grid(row=2, column=1, sticky="w", **pad)

        ttk.Label(win, text="Arrow width:").grid(row=3, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=width_var, width=12).grid(row=3, column=1, sticky="w", **pad)

        ttk.Label(win, text="Origin sphere radius:").grid(row=4, column=0, sticky="e", **pad)
        ttk.Entry(win, textvariable=sphere_var, width=12).grid(row=4, column=1, sticky="w", **pad)

        ttk.Label(
            win,
            text="Colors: X red, Y yellow, Z blue. Arrow head radius is 2.5 x arrow width.",
            foreground="#666",
            wraplength=460,
        ).grid(row=5, column=0, columnspan=3, sticky="w", padx=10, pady=(4, 8))

        buttons = ttk.Frame(win)
        buttons.grid(row=6, column=0, columnspan=3, sticky="e", padx=10, pady=(4, 10))
        ttk.Button(buttons, text="Cancel", command=win.destroy).pack(side="left", padx=6)
        ttk.Button(buttons, text="Write", command=write_file).pack(side="left", padx=6)

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
            self.param_values_var.set("Z-DNA uses DSSR fiber; the 12-parameter table is not used.")
            return

        parts = []
        changed = []
        for key, default in zip(PARAM_KEYS, defaults):
            raw = store.get(key, "").strip()
            value = raw if raw else f"{default:.4f}"
            try:
                numeric = float(value)
                shown = f"{numeric:.4f}"
                if raw and abs(numeric - float(default)) > 1e-10:
                    changed.append(key)
            except Exception:
                shown = value
            parts.append(f"{key}={shown}")
        prefix = f"{na_type} DSSR table values: "
        if changed:
            prefix += "customized fields: " + ", ".join(changed) + ". "
        else:
            prefix += "default values. "
        self.param_values_var.set(prefix + "    ".join(parts))

    def _update_sequence_length(self) -> None:
        na_type = self.na_type_var.get().strip()
        if na_type == "Z-DNA":
            text = self.z_len_var.get().strip()
            if not text:
                self.length_var.set("Length: 0 bp")
                return
            try:
                n = int(text)
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
        win.geometry("800x500+180+120")
        win.minsize(700, 430)
        win.transient(self)
        win.grab_set()

        ttk.Label(
            win,
            text=(
                "Default values are filled automatically when no previous value exists. "
                "Previously saved values are kept. Values are written to four digits after the decimal."
            ),
            wraplength=740,
            foreground="#555",
            justify="left",
        ).grid(row=0, column=0, columnspan=4, sticky="w", padx=14, pady=(14, 8))

        local_vars: Dict[str, tk.StringVar] = {}
        for idx, (key, default) in enumerate(zip(PARAM_KEYS, defaults)):
            col_group = 0 if idx < 6 else 2
            row = idx + 1 if idx < 6 else idx - 5
            ttk.Label(win, text=PARAM_LABELS.get(key, key) + ":").grid(
                row=row, column=col_group, sticky="e", padx=10, pady=6
            )
            initial = store.get(key, "").strip() or f"{default:.4f}"
            local_vars[key] = tk.StringVar(value=initial)
            ttk.Entry(win, textvariable=local_vars[key], width=18).grid(
                row=row, column=col_group + 1, sticky="w", padx=10, pady=6
            )

        btns = ttk.Frame(win)
        btns.grid(row=8, column=0, columnspan=4, sticky="e", padx=12, pady=12)

        def clear_all() -> None:
            for key in PARAM_KEYS:
                local_vars[key].set("")

        def restore_defaults() -> None:
            for key, default in zip(PARAM_KEYS, defaults):
                local_vars[key].set(f"{default:.4f}")

        def save_close() -> None:
            for key in PARAM_KEYS:
                txt = local_vars[key].get().strip()
                if txt:
                    try:
                        float(txt)
                    except Exception:
                        messagebox.showerror("Invalid value", f"{key} must be a number.", parent=win)
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
                try:
                    value = float(txt)
                except Exception as exc:
                    raise ValueError(f"{key} must be a number.") from exc
                if abs(value - float(default)) > 1e-10:
                    overrides[key] = value
        return overrides

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
                    length_value = int(self.z_len_var.get().strip())
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
                roll_deg=self.roll_var.get(),
                phi_deg=self.phi_var.get(),
                theta_deg=self.theta_var.get(),
                tx=self.x_var.get(),
                ty=self.y_var.get(),
                tz=self.z_var.get(),
                delta_z=self.delta_z_var.get(),
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
