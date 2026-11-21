"""Wrapper for building B-DNA/P-DNA junction constructs from bnp_na."""
from __future__ import annotations

import argparse
import contextlib
import io
import shlex
from pathlib import Path
from typing import Dict, Sequence, Union

from .build_common import PipelineError
from .make_BP import (
    AXIS_SOURCES,
    DEFAULT_OUTPUT_NAME,
    DEFAULT_TRANSITION,
    TRANSITIONS,
    raw_output_path,
    resolve_output_path,
    run_pipeline,
)
from .make_BZV2_3 import normalize_axis_mode


def _expanded_path(path: Union[str, Path]) -> Path:
    value = Path(path).expanduser()
    if not value.is_absolute():
        value = value.resolve()
    return value


def build_bp_structure(
    input_files: Sequence[Union[str, Path]],
    out_pdb: Union[str, Path],
    *,
    types: str = "",
    axis_mode: str = "codirectional",
    axis_source: str = "auto",
    transition: str = DEFAULT_TRANSITION,
) -> Dict[str, object]:
    """Build a multi-segment B/P DNA model and capture the make_BP log.

    ``transition`` is the number of transition nucleotides per strand at each
    junction, "2" (compact core) or "3" (longer, forward-winding transition).
    """

    paths = [_expanded_path(path) for path in input_files if str(path).strip()]
    if len(paths) < 2:
        raise PipelineError("B-P building requires at least two input files: one B-DNA and one P-DNA.")
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise PipelineError("B-P input file(s) not found:\n" + "\n".join(missing))

    try:
        normalized_axis_mode = normalize_axis_mode(axis_mode)
    except Exception as exc:
        raise PipelineError(str(exc)) from exc
    axis_source_arg = (axis_source or "auto").strip().lower()
    if axis_source_arg not in AXIS_SOURCES:
        raise PipelineError(f"Unknown axis source: {axis_source}")
    types_arg = "".join(ch for ch in (types or "").upper() if ch in "BP")
    transition_arg = str(transition or DEFAULT_TRANSITION).strip()
    if transition_arg not in TRANSITIONS:
        raise PipelineError(f"Transition must be one of {', '.join(TRANSITIONS)} nucleotides, got {transition!r}")

    out_path = Path(out_pdb).expanduser() if str(out_pdb).strip() else paths[0].parent / DEFAULT_OUTPUT_NAME
    if not out_path.is_absolute():
        out_path = out_path.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    args = argparse.Namespace(
        dna_files=[str(path) for path in paths],
        out=str(out_path),
        types=types_arg,
        axis_mode=normalized_axis_mode,
        axis_source=axis_source_arg,
        transition=transition_arg,
    )
    cli_args = ["--out", str(out_path), "--axis-mode", normalized_axis_mode, "--axis-source", axis_source_arg,
                "--transition", transition_arg]
    if types_arg:
        cli_args += ["--types", types_arg]
    cli_args.extend(str(path) for path in paths)

    buf = io.StringIO()
    buf.write(
        "Equivalent CLI command:\n  "
        + " ".join(shlex.quote(part) for part in ["python3", "bnp_na_lib/make_BP.py", *cli_args])
        + "\n\n"
    )
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            result = run_pipeline(args)
    except SystemExit as exc:
        message = str(exc) or "B-P structure build stopped."
        raise PipelineError(message, buf.getvalue() + f"\nERROR: {message}\n") from exc
    except Exception as exc:
        message = f"B-P structure build failed: {exc}"
        raise PipelineError(message, buf.getvalue() + f"\nERROR: {type(exc).__name__}: {exc}\n") from exc

    resolved_out = resolve_output_path(args.out, paths)
    return {
        "na_type": "B-P DNA",
        "input_files": paths,
        "types": "".join(result.get("types", [])),
        "pdb_out": resolved_out,
        "pdb_raw": raw_output_path(resolved_out),
        "axis_mode": normalized_axis_mode,
        "axis_source": axis_source_arg,
        "transition": transition_arg,
        "log_text": buf.getvalue(),
    }
