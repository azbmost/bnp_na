"""Wrapper for building B/Z DNA junction constructs from bnp_na."""
from __future__ import annotations

import argparse
import contextlib
import io
import shlex
from pathlib import Path
from typing import Dict, Sequence, Union

from build_common import PipelineError
from make_BZV2_3 import (
    DEFAULT_OUTPUT_NAME,
    effective_axis_source,
    normalize_axis_mode,
    raw_output_path,
    resolve_output_path,
    run_pipeline,
)


def _expanded_path(path: Union[str, Path]) -> Path:
    value = Path(path).expanduser()
    if not value.is_absolute():
        value = value.resolve()
    return value


def build_bz_structure(
    input_files: Sequence[Union[str, Path]],
    out_pdb: Union[str, Path],
    *,
    axis_mode: str = "codirectional",
    axis_source: str = "auto",
    auto_trim_z: bool = True,
) -> Dict[str, object]:
    """Build a multi-segment B/Z DNA model and capture the make_BZ log."""

    paths = [_expanded_path(path) for path in input_files if str(path).strip()]
    if len(paths) < 2:
        raise PipelineError("B-Z building requires at least two input files: B1 and Z1.")
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise PipelineError("B-Z input file(s) not found:\n" + "\n".join(missing))

    try:
        normalized_axis_mode = normalize_axis_mode(axis_mode)
        axis_source_arg = (axis_source or "auto").strip().lower()
        effective_axis_source(argparse.Namespace(axis_source=axis_source_arg))
    except Exception as exc:
        raise PipelineError(str(exc)) from exc

    out_path = Path(out_pdb).expanduser()
    if not str(out_path).strip():
        out_path = paths[0].parent / DEFAULT_OUTPUT_NAME
    if not out_path.is_absolute():
        out_path = out_path.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    args = argparse.Namespace(
        dna_files=[str(path) for path in paths],
        out=str(out_path),
        axis_mode=normalized_axis_mode,
        collinear=None,
        axis_source=axis_source_arg,
        dist_cutoff=2.2,
        no_z_auto_trim=not bool(auto_trim_z),
    )
    cli_args = [
        "--out",
        str(out_path),
        "--axis-mode",
        normalized_axis_mode,
        "--axis-source",
        axis_source_arg,
    ]
    if not auto_trim_z:
        cli_args.append("--no-z-auto-trim")
    cli_args.extend(str(path) for path in paths)

    buf = io.StringIO()
    buf.write(
        "Equivalent CLI command:\n  "
        + " ".join(shlex.quote(part) for part in ["python3", "bnp_na_lib/make_BZV2_3.py", *cli_args])
        + "\n\n"
    )
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            run_pipeline(args)
    except SystemExit as exc:
        message = str(exc) or "B-Z structure build stopped."
        log_text = buf.getvalue() + f"\nERROR: {message}\n"
        raise PipelineError(message, log_text) from exc
    except Exception as exc:
        message = f"B-Z structure build failed: {exc}"
        log_text = buf.getvalue() + f"\nERROR: {type(exc).__name__}: {exc}\n"
        raise PipelineError(message, log_text) from exc

    resolved_out = resolve_output_path(args.out, paths)
    return {
        "na_type": "B-Z DNA",
        "input_files": paths,
        "pdb_out": resolved_out,
        "pdb_raw": raw_output_path(resolved_out),
        "axis_mode": normalized_axis_mode,
        "axis_source": axis_source_arg,
        "auto_trim_z": bool(auto_trim_z),
        "log_text": buf.getvalue(),
    }
