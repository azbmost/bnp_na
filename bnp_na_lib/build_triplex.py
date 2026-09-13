"""Wrapper for converting duplex PDB models into triplex PDB models."""
from __future__ import annotations

import shlex
import shutil
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from build_common import (
    PipelineError,
    command_to_text,
    expected_phenix_minimized_path,
    run_phenix_minimization,
    stage_params_to_output_dir,
)
from convert_to_triplex_pdbV2_1 import (
    convert_duplex_to_triplex,
    default_output_path,
    format_chain_sequence_info,
    format_result_summary,
    format_selection_preview,
    normalize_mode,
)
from regularize_phosphates import default_regularized_output_path, regularize_phosphates


#: Intermediates from the minimization/regularization steps are kept beside the
#: requested output file, so the named output stays the final model.
TMP_DIR_NAME = "triplex_gen_tmp"

#: Used when minimization is requested without naming a params file.
DEFAULT_PARAMS_FILE = Path(__file__).resolve().parent / "min_P_C5.params"


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


def resolve_triplex_paths(
    output_path: Union[str, Path],
    *,
    post_process: bool,
) -> Tuple[Path, Optional[Path], Path]:
    """Return ``(final_path, tmp_dir, convert_target)`` for one conversion.

    Without post-processing the conversion writes straight to the requested
    path, which is what the tool has always done. With it, the conversion and
    every later stage live in ``triplex_gen_tmp`` beside that path, and the
    requested path receives the finished model at the end.
    """
    out_path = _expanded_path(output_path)
    if not post_process:
        return out_path, None, out_path
    tmp_dir = out_path.parent / TMP_DIR_NAME
    convert_target = tmp_dir / f"{out_path.stem}_triplex_raw{out_path.suffix or '.pdb'}"
    return out_path, tmp_dir, convert_target


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
    run_phenix: bool = True,
    params_file: Optional[Union[str, Path]] = None,
    run_regularize_phosphates: Optional[bool] = None,
) -> Dict[str, object]:
    """Convert a duplex PDB into a triplex and return paths plus log text.

    Minimization and phosphate regularization are both on by default, matching
    the triplex converter dialog. With either enabled the converted triplex is
    post-processed and ``output_path`` receives the final model, while every
    intermediate, including the raw conversion and anything Phenix writes
    alongside it, goes into a ``triplex_gen_tmp`` folder next to that output.
    Pass ``run_phenix=False, run_regularize_phosphates=False`` to write the raw
    conversion straight to ``output_path``. ``params_file`` falls back to the
    bundled ``min_P_C5.params`` when minimization runs without one.
    """

    duplex = _expanded_path(duplex_path)
    if not duplex.exists():
        raise PipelineError(f"Input duplex PDB not found: {duplex}")
    if not strand_i_chain.strip():
        raise PipelineError("Please specify strand I, the purine chain in the input duplex.")

    if run_regularize_phosphates is None:
        run_regularize_phosphates = bool(run_phenix)

    normalized_mode = normalize_mode(mode)
    if output_path is None or not str(output_path).strip():
        out_path = default_triplex_output_path(duplex)
    else:
        out_path = _expanded_path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    post_process = bool(run_phenix or run_regularize_phosphates)
    params_source = None
    if run_phenix:
        params_source = DEFAULT_PARAMS_FILE if params_file is None else params_file
        if not str(params_source).strip():
            raise PipelineError(
                "phenix.geometry_minimization was requested, but no params file was specified."
            )
        params_source = _expanded_path(params_source)
        if not params_source.exists():
            raise PipelineError(f"Params file not found: {params_source}")
    out_path, tmp_dir, convert_target = resolve_triplex_paths(out_path, post_process=post_process)
    if tmp_dir is not None:
        tmp_dir.mkdir(parents=True, exist_ok=True)

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
        str(convert_target),
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
            output_path=str(convert_target),
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

    pdb_raw = Path(result.output_path)
    pdb_min = None
    pdb_regularized = None
    final_source = pdb_raw

    if run_phenix:
        assert tmp_dir is not None  # post_process is true whenever run_phenix is
        try:
            staged_params = stage_params_to_output_dir(params_source, tmp_dir)
        except Exception as exc:
            raise PipelineError(str(exc), "\n".join(log_parts)) from exc

        ok_phx, out_phx, cmd_phx = run_phenix_minimization(final_source, staged_params)
        log_parts += [
            "",
            "=== phenix.geometry_minimization ===",
            f"Params : {staged_params}",
            f"Command: (cwd={tmp_dir}) {command_to_text(cmd_phx)}",
            f"Status : {'OK' if ok_phx else 'FAILED'}",
            f"Output :\n{out_phx}",
        ]
        if not ok_phx:
            raise PipelineError("phenix.geometry_minimization failed.", "\n".join(log_parts))

        pdb_min = expected_phenix_minimized_path(final_source)
        if not pdb_min.exists():
            raise PipelineError(
                f"Expected minimized PDB was not found: {pdb_min}", "\n".join(log_parts)
            )
        log_parts.append(f"Minimized PDB: {pdb_min}")
        final_source = pdb_min
    else:
        log_parts += ["", "=== phenix.geometry_minimization ===", "Skipped by user option."]

    if run_regularize_phosphates:
        try:
            regularize_result = regularize_phosphates(
                final_source,
                default_regularized_output_path(final_source),
            )
        except Exception as exc:
            raise PipelineError(
                f"Phosphate regularization failed: {exc}", "\n".join(log_parts)
            ) from exc
        pdb_regularized = Path(regularize_result.output_pdb)
        final_source = pdb_regularized
        log_parts += ["", regularize_result.log_text]
    else:
        log_parts += ["", "=== Regularize phosphates ===", "Skipped by user option."]

    if post_process:
        try:
            shutil.copyfile(final_source, out_path)
        except Exception as exc:
            raise PipelineError(
                f"Could not write the final triplex PDB to {out_path}: {exc}",
                "\n".join(log_parts),
            ) from exc
        log_parts += [
            "",
            "=== Final triplex PDB ===",
            f"Intermediates: {tmp_dir}",
            f"Final model  : {out_path}",
        ]

    return {
        "na_type": "Triplex DNA",
        "pdb_out": out_path,
        "pdb_raw": pdb_raw,
        "pdb_minimized": pdb_min,
        "pdb_regularized": pdb_regularized,
        "tmp_dir": tmp_dir,
        "mode": result.mode,
        "motif_label": result.motif_label,
        "strand_i": result.strand_i,
        "strand_ii": result.strand_ii,
        "strand_iii": result.strand_iii,
        "length": len(result.strand_iii_resseqs),
        "summary": summary,
        "log_text": "\n".join(log_parts),
    }
