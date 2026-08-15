#!/usr/bin/env python3
"""A-RNA builder for bnp_na.

Pipeline: write DSSR helical-parameter table -> DSSR rebuild with
--backbone=RNA --par-type=heli -> normalize PDB naming -> optional
phenix.geometry_minimization -> optional phosphate regularization ->
DSSR --more axis extraction -> align to +Z.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Union

from align2z import align_pdb_to_Z, format_alignment_report
from build_common import (
    DEFAULT_PARAMS,
    PARAM_KEYS,
    PipelineError,
    command_to_text,
    expected_phenix_minimized_path,
    expand_sequence,
    sanitize_basename,
    sequence_alphabet,
    stage_params_to_output_dir,
    run_dssr_rebuild,
    run_phenix_minimization,
    write_helical_table,
)
from pdb_name_standard import normalize_pdb_naming as normalize_nucleotide_pdb_naming
from regularize_phosphates import default_regularized_output_path, regularize_phosphates


BACKBONE = "A-RNA"
DSSR_BACKBONE = "RNA"

TEMPLATES = {
    base: (tag, list(DEFAULT_PARAMS[BACKBONE]))
    for base, tag in {
        "A": "A-U",
        "U": "U-A",
        "C": "C-G",
        "G": "G-C",
    }.items()
}


def _log_default_params() -> str:
    return ", ".join(f"{key}={value:.4f}" for key, value in zip(PARAM_KEYS, DEFAULT_PARAMS[BACKBONE]))


def build_arna(
    seq_53: str,
    base_name: Optional[str],
    output_dir: Union[str, Path],
    param_overrides: Optional[Dict[str, float]] = None,
    deleteH: bool = False,
    run_phenix: bool = False,
    params_file: Optional[Union[str, Path]] = None,
    run_regularize_phosphates: Optional[bool] = None,
) -> Dict[str, object]:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logs = []
    if run_regularize_phosphates is None:
        run_regularize_phosphates = bool(run_phenix)

    try:
        seq_expanded = expand_sequence(seq_53, alphabet=sequence_alphabet(BACKBONE))
    except Exception as exc:
        raise PipelineError(f"Invalid sequence: {exc}", str(exc)) from exc
    if not seq_expanded:
        raise PipelineError("Empty sequence.")

    base = sanitize_basename(base_name or f"A-RNA{len(seq_expanded)}")
    if not base:
        raise PipelineError("Invalid helix name after sanitization.")

    par_file = out_dir / f"{base}.txt"
    pdb_rebuild = out_dir / f"{base}-rb.pdb"

    try:
        table_path, seq_expanded = write_helical_table(
            seq_expanded,
            par_file,
            na_type=BACKBONE,
            param_overrides=param_overrides,
        )
        logs += [
            f"Wrote A-RNA helical-parameter table for {len(seq_expanded)} bp:",
            str(table_path),
            "Default 12 parameters used when not overridden:",
            _log_default_params(),
        ]
    except Exception as exc:
        raise PipelineError(f"Failed to write helical table: {exc}", "\n".join(logs)) from exc

    ok_dssr, out_dssr, cmd_dssr = run_dssr_rebuild(par_file, pdb_rebuild, cwd=out_dir, backbone=DSSR_BACKBONE)
    logs += [
        "\n=== DSSR rebuild ===",
        f"Command: (cwd={out_dir}) {command_to_text(cmd_dssr)}",
        f"Status : {'OK' if ok_dssr else 'FAILED'}",
        f"Output :\n{out_dssr}",
    ]
    if not ok_dssr or not pdb_rebuild.exists():
        raise PipelineError("DSSR rebuild failed.", "\n".join(logs))

    pdb_norm_target = out_dir / f"{base}-rb_out.pdb"
    try:
        pdb_norm = Path(
            normalize_nucleotide_pdb_naming(
                str(pdb_rebuild),
                str(pdb_norm_target),
                deleteH=deleteH,
            )
        )
        logs += [
            "\n=== Normalize PDB names (pdb_name_standard.py) ===",
            f"Input : {pdb_rebuild}",
            f"Output: {pdb_norm}",
            f"deleteH: {deleteH}",
        ]
    except Exception as exc:
        raise PipelineError(f"PDB naming normalization failed: {exc}", "\n".join(logs)) from exc

    align_source = pdb_norm
    pdb_min = None
    staged_params = None
    if run_phenix:
        if not params_file:
            raise PipelineError("phenix.geometry_minimization was requested, but no params file was specified.", "\n".join(logs))
        try:
            staged_params = stage_params_to_output_dir(params_file, out_dir)
        except Exception as exc:
            raise PipelineError(str(exc), "\n".join(logs)) from exc

        ok_phx, out_phx, cmd_phx = run_phenix_minimization(pdb_norm, staged_params)
        logs += [
            "\n=== phenix.geometry_minimization ===",
            f"Params : {staged_params}",
            f"Command: (cwd={out_dir}) {command_to_text(cmd_phx)}",
            f"Status : {'OK' if ok_phx else 'FAILED'}",
            f"Output :\n{out_phx}",
        ]
        if not ok_phx:
            raise PipelineError("phenix.geometry_minimization failed.", "\n".join(logs))

        pdb_min = expected_phenix_minimized_path(pdb_norm)
        if not pdb_min.exists():
            raise PipelineError(f"Expected minimized PDB was not found: {pdb_min}", "\n".join(logs))
        logs.append(f"Minimized PDB: {pdb_min}")
        align_source = pdb_min
    else:
        logs += [
            "\n=== phenix.geometry_minimization ===",
            "Skipped by user option.",
        ]

    pdb_regularized = None
    if run_regularize_phosphates:
        try:
            regularize_result = regularize_phosphates(
                align_source,
                default_regularized_output_path(align_source),
            )
            pdb_regularized = regularize_result.output_pdb
            align_source = pdb_regularized
            logs += ["\n" + regularize_result.log_text]
        except Exception as exc:
            raise PipelineError(f"Phosphate regularization failed: {exc}", "\n".join(logs)) from exc
    else:
        logs += ["\n=== Regularize phosphates ===", "Skipped by user option."]

    pdb_aligned = align_source.with_name(align_source.stem + "_aligned2Z" + align_source.suffix)
    try:
        align_report = align_pdb_to_Z(str(align_source), out_pdb=str(pdb_aligned), cwd=out_dir)
        logs += ["\n" + format_alignment_report(align_report)]
    except Exception as exc:
        raise PipelineError(f"DSSR align-to-Z failed: {exc}", "\n".join(logs)) from exc

    return {
        "na_type": BACKBONE,
        "base_name": base,
        "length": len(seq_expanded),
        "par_file": par_file,
        "pdb_rebuild": pdb_rebuild,
        "pdb_normalized": pdb_norm,
        "pdb_minimized": pdb_min,
        "pdb_regularized": pdb_regularized,
        "pdb_aligned": pdb_aligned,
        "run_phenix": run_phenix,
        "run_regularize_phosphates": run_regularize_phosphates,
        "params_file": staged_params,
        "log_text": "\n".join(logs),
    }


# Backward-compatible alias for older imports while bnp_na.py uses build_arna.
build_arna_align2z = build_arna
