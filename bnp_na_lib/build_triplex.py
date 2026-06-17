"""Wrapper for converting duplex PDB models into triplex PDB models."""
from __future__ import annotations

import shlex
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from build_common import PipelineError
from convert_to_triplex_pdbV2_1 import (
    convert_duplex_to_triplex,
    default_output_path,
    format_chain_sequence_info,
    format_result_summary,
    format_selection_preview,
    normalize_mode,
)


def _expanded_path(path: Union[str, Path]) -> Path:
    value = Path(path).expanduser()
    if not value.is_absolute():
        value = value.resolve()
    return value


def default_triplex_output_path(duplex_path: Union[str, Path]) -> Path:
    return _expanded_path(default_output_path(str(duplex_path)))


def describe_triplex_input(
    duplex_path: Union[str, Path],
    *,
    strand_i_chain: str = "",
    residue_range: Optional[Tuple[int, int]] = None,
    mode: str = "antiparallel",
) -> str:
    """Return chain/sequence information and an optional selection preview."""

    path = _expanded_path(duplex_path)
    if not path.exists():
        raise PipelineError(f"Input duplex PDB not found: {path}")
    text = format_chain_sequence_info(str(path))
    if strand_i_chain and residue_range is not None:
        text += format_selection_preview(
            pdb_path=str(path),
            strand_i_chain=strand_i_chain,
            residue_range=residue_range,
            mode=mode,
        )
    return text


def build_triplex_from_duplex(
    duplex_path: Union[str, Path],
    output_path: Optional[Union[str, Path]],
    *,
    strand_i_chain: str,
    residue_range: Tuple[int, int],
    mode: str,
    strand_ii_chain: Optional[str] = None,
    strand_iii_chain: Optional[str] = None,
    strand_iii_start_resseq: int = 1,
) -> Dict[str, object]:
    """Convert a duplex PDB into a triplex and return paths plus log text."""

    duplex = _expanded_path(duplex_path)
    if not duplex.exists():
        raise PipelineError(f"Input duplex PDB not found: {duplex}")
    if not strand_i_chain.strip():
        raise PipelineError("Please specify strand I, the purine chain in the input duplex.")

    normalized_mode = normalize_mode(mode)
    if output_path is None or not str(output_path).strip():
        out_path = default_triplex_output_path(duplex)
    else:
        out_path = _expanded_path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    start_res, end_res = residue_range
    cli_args = [
        "python3",
        "bnp_na_lib/convert_to_triplex_pdbV2_1.py",
        str(duplex),
        "--strand-I",
        strand_i_chain.strip(),
        "--range",
        f"{start_res}:{end_res}",
        "--mode",
        normalized_mode,
        "--strand-III-start",
        str(int(strand_iii_start_resseq)),
        "--out",
        str(out_path),
    ]
    if strand_ii_chain:
        cli_args.extend(["--strand-II", strand_ii_chain])
    if strand_iii_chain:
        cli_args.extend(["--strand-III", strand_iii_chain])

    log_parts = [
        "Equivalent CLI command:",
        "  " + " ".join(shlex.quote(part) for part in cli_args),
        "",
    ]

    try:
        log_parts.append(
            describe_triplex_input(
                duplex,
                strand_i_chain=strand_i_chain.strip(),
                residue_range=residue_range,
                mode=normalized_mode,
            )
        )
        log_parts.append("")
        log_parts.append("Running triplex conversion...")
        result = convert_duplex_to_triplex(
            duplex_path=str(duplex),
            output_path=str(out_path),
            strand_i_chain=strand_i_chain.strip(),
            residue_range=residue_range,
            mode=normalized_mode,
            strand_ii_chain=(strand_ii_chain.strip() if strand_ii_chain else None),
            strand_iii_chain=(strand_iii_chain.strip() if strand_iii_chain else None),
            strand_iii_start_resseq=int(strand_iii_start_resseq),
        )
        summary = format_result_summary(result)
        log_parts.extend(["", summary])
    except PipelineError:
        raise
    except Exception as exc:
        log_text = "\n".join(log_parts + ["", f"ERROR: {type(exc).__name__}: {exc}"])
        raise PipelineError(f"Triplex conversion failed: {exc}", log_text) from exc

    return {
        "na_type": "Triplex DNA",
        "pdb_out": Path(result.output_path),
        "mode": result.mode,
        "motif_label": result.motif_label,
        "strand_i": result.strand_i,
        "strand_ii": result.strand_ii,
        "strand_iii": result.strand_iii,
        "length": len(result.strand_iii_resseqs),
        "summary": summary,
        "log_text": "\n".join(log_parts),
    }
