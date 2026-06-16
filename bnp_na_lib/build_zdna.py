#!/usr/bin/env python3
"""Z-DNA builder for bnp_na.

Z-DNA is still generated with DSSR fiber because the current GUI specifies only
helix length. Axis extraction/alignment now uses DSSR --more via align2z for accurate point-one/point-two endpoints.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Union

from align2z import align_pdb_to_Z, format_alignment_report
from build_common import PipelineError, command_to_text, run_dssr_fiber_z, sanitize_basename
from pdb_make_dna_v3_2 import normalize_pdb_naming as normalize_nucleotide_pdb_naming


BACKBONE = "Z-DNA"


def build_zdna(
    length: int,
    base_name: Optional[str],
    output_dir: Union[str, Path],
    deleteH: bool = False,
) -> Dict[str, object]:
    """Build an ideal Z-DNA fiber of even base-pair length and align to +Z."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logs = []

    try:
        n_bp = int(length)
    except Exception as exc:
        raise PipelineError("Z-DNA helix length must be a positive even integer.") from exc
    if n_bp <= 0 or n_bp % 2 != 0:
        raise PipelineError("Z-DNA helix length must be a positive even integer.")

    repeat = n_bp // 2
    base = sanitize_basename(base_name or f"Z-DNA{n_bp}")
    if not base:
        raise PipelineError("Invalid helix name after sanitization.")

    fiber_pdb = out_dir / f"{base}.pdb"
    ok_fb, out_fb, cmd_fb = run_dssr_fiber_z(repeat, fiber_pdb, cwd=out_dir)
    logs += [
        "=== DSSR fiber (Z-DNA) ===",
        f"Command: (cwd={out_dir}) {command_to_text(cmd_fb)}",
        f"Status : {'OK' if ok_fb else 'FAILED'}",
        f"Output :\n{out_fb}",
    ]
    if not ok_fb or not fiber_pdb.exists():
        raise PipelineError("DSSR Z-DNA fiber failed.", "\n".join(logs))

    fiber_norm = out_dir / f"{base}_norm.pdb"
    try:
        fiber_norm = Path(
            normalize_nucleotide_pdb_naming(
                str(fiber_pdb),
                str(fiber_norm),
                deleteH=deleteH,
            )
        )
        logs += [
            "\n=== Normalize PDB naming v3.2 ===",
            f"Input : {fiber_pdb}",
            f"Output: {fiber_norm}",
            f"deleteH: {deleteH}",
        ]
    except Exception as exc:
        raise PipelineError(f"PDB naming normalization failed: {exc}", "\n".join(logs)) from exc

    aligned_pdb = out_dir / f"{base}_aligned2Z.pdb"
    try:
        align_report = align_pdb_to_Z(str(fiber_norm), out_pdb=str(aligned_pdb), cwd=out_dir)
        logs += ["\n" + format_alignment_report(align_report)]
    except Exception as exc:
        raise PipelineError(f"DSSR align-to-Z failed: {exc}", "\n".join(logs)) from exc

    return {
        "na_type": BACKBONE,
        "base_name": base,
        "length": n_bp,
        "repeat": repeat,
        "pdb_fiber": fiber_pdb,
        "pdb_normalized": fiber_norm,
        "pdb_aligned": aligned_pdb,
        "log_text": "\n".join(logs),
    }


# Backward-compatible alias for older imports while bnp_na.py uses build_zdna.
build_zdna_align2z = build_zdna
