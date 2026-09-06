#!/usr/bin/env python3
"""Find the B-DNA X-disp that places opposing phosphate P atoms across the helix axis.

For an ``N``-bp helix the search compares ``A_i.P`` with ``B_(N+2-i).P`` for
``i = 2..N``. After alignment to +Z the DSSR helix axis is the z-axis through
``x = 0, y = 0``, so an opposing pair lies across the axis when the signed xy
cross product ``A_x * B_y - A_y * B_x`` is zero, the xy-projected P-P line
passes through the axis, and the around-axis angle is 180 degrees.

Two answers exist for the same helix because Phenix minimization relaxes the
phosphate/sugar geometry after the DSSR rebuild:

* raw DSSR rebuild only, and
* the bnp_na default workflow of Phenix minimization plus phosphate
  regularization.

This module grew out of ``test/oppo_phos_azbmost/find_opposing_phosphate_xdisp.py``.
For the ``A31`` reference helix at the bnp_na B-DNA defaults it reproduces
``X-disp = 3.0872`` (raw DSSR) and ``X-disp = 2.7992`` (Phenix + regularized).
"""
from __future__ import annotations

import argparse
import math
import tempfile
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

try:  # bnp_na_lib on sys.path, as the GUI and the bundled CLIs arrange it
    from build_common import PARAM_KEYS
except ImportError:  # imported as part of the bnp_na_lib package
    from bnp_na_lib.build_common import PARAM_KEYS


TOOL_VERSION = "V13.17"

#: X-disp is written to the DSSR helical table with four decimals, so the
#: search resolution can never be finer than this.
X_DISP_QUANTUM = 0.0001

DEFAULT_DSSR_X_MIN = 0.0
DEFAULT_DSSR_X_MAX = 8.0
DEFAULT_DSSR_X_STEP = 0.25

#: Half-width of the final four-decimal grid sweep, in angstrom. The raw DSSR
#: rebuild is fast enough to afford a wide sweep; every Phenix point costs a
#: full minimization, so its sweep is deliberately narrow.
DEFAULT_DSSR_GRID_HALFWIDTH = 0.0020
DEFAULT_PHENIX_GRID_HALFWIDTH = 0.0002

#: Safety cap on Phenix evaluations, which cost roughly 15 s each.
DEFAULT_MAX_PHENIX_BUILDS = 24

BISECTION_STEPS = 18
SECANT_STEPS = 8

MODE_DSSR = "dssr"
MODE_PHENIX = "phenix"

PAIR_FIELDNAMES = [
    "A_res",
    "B_res",
    "A_x",
    "A_y",
    "A_z",
    "B_x",
    "B_y",
    "B_z",
    "signed_cross_xy",
    "axis_line_distance_A",
    "radial_angle_deg",
    "axis_fraction_on_segment",
    "axis_between_atoms",
]

SUMMARY_FIELDNAMES = [
    "mode",
    "stage",
    "selected",
    "x_disp",
    "fixed_parameter_overrides",
    "mean_signed_cross_xy",
    "mean_abs_cross_xy",
    "max_abs_cross_xy",
    "rms_axis_line_distance_A",
    "mean_axis_line_distance_A",
    "max_axis_line_distance_A",
    "mean_radial_angle_deg",
    "min_radial_angle_deg",
    "max_radial_angle_deg",
    "min_axis_fraction_on_segment",
    "max_axis_fraction_on_segment",
]


class OpposingPhosphateError(Exception):
    """Raised when the opposing-phosphate X-disp search cannot be completed."""


class OpposingPhosphateCancelled(OpposingPhosphateError):
    """Raised when a caller-supplied cancel callback stops the search."""


def quantize_x(value: float) -> float:
    """Match bnp_na's four-decimal DSSR table formatting."""
    return float(Decimal(str(float(value))).quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP))


def x_label(value: float) -> str:
    """Filename-safe tag for an X-disp value, for example ``3p0872``."""
    return f"{quantize_x(value):.4f}".replace("-", "m").replace(".", "p")


def overrides_text(overrides: Optional[Dict[str, float]]) -> str:
    if not overrides:
        return "none"
    return "; ".join(f"{key}={value:.4f}" for key, value in sorted(overrides.items()))


def validate_fixed_overrides(overrides: Optional[Dict[str, object]]) -> Dict[str, float]:
    """Return a clean copy of non-X parameter overrides, or raise."""
    clean: Dict[str, float] = {}
    for key, value in (overrides or {}).items():
        if key not in PARAM_KEYS:
            raise OpposingPhosphateError(
                f"Unknown parameter {key!r}; valid names are: {', '.join(PARAM_KEYS)}"
            )
        if key == "X-disp":
            raise OpposingPhosphateError("X-disp is the searched parameter and cannot be held fixed.")
        try:
            clean[key] = float(value)
        except (TypeError, ValueError) as exc:
            raise OpposingPhosphateError(f"Parameter {key!r} value must be numeric: {value!r}") from exc
    return clean


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------
def parse_p_atoms(pdb_path: Union[str, Path]) -> Dict[Tuple[str, int], Tuple[float, float, float]]:
    """Return ``{(chain, resseq): (x, y, z)}`` for every phosphorus atom."""
    atoms: Dict[Tuple[str, int], Tuple[float, float, float]] = {}
    with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "P":
                continue
            chain = line[21].strip()
            resseq = int(line[22:26])
            atoms[(chain, resseq)] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    return atoms


def phosphate_pair_rows(pdb_path: Union[str, Path]) -> List[Dict[str, object]]:
    """Return per-pair metrics for ``A_i.P`` versus ``B_(N+2-i).P``, ``i = 2..N``."""
    p_atoms = parse_p_atoms(pdb_path)
    a_res = sorted(res for chain, res in p_atoms if chain == "A")
    b_res = sorted(res for chain, res in p_atoms if chain == "B")
    if not a_res or not b_res:
        raise OpposingPhosphateError(f"Could not find both chain A and chain B P atoms in {pdb_path}")
    n_bp = max(a_res)

    rows: List[Dict[str, object]] = []
    for i in range(2, n_bp + 1):
        j = n_bp + 2 - i
        try:
            ax, ay, az = p_atoms[("A", i)]
            bx, by, bz = p_atoms[("B", j)]
        except KeyError as exc:
            raise OpposingPhosphateError(f"Missing expected opposing P atom for A{i}/B{j}") from exc

        cross = ax * by - ay * bx
        dot = ax * bx + ay * by
        ra = math.hypot(ax, ay)
        rb = math.hypot(bx, by)
        if ra == 0.0 or rb == 0.0:
            raise OpposingPhosphateError(f"P atom for A{i}/B{j} lies on the helix axis.")
        angle = math.degrees(math.acos(max(-1.0, min(1.0, dot / (ra * rb)))))

        vx = bx - ax
        vy = by - ay
        segment_sq = (vx * vx) + (vy * vy)
        if segment_sq == 0.0:
            raise OpposingPhosphateError(f"Opposing P atoms A{i}/B{j} project onto the same xy point.")
        distance = abs(cross) / math.sqrt(segment_sq)
        t = -((ax * vx) + (ay * vy)) / segment_sq

        rows.append(
            {
                "A_res": i,
                "B_res": j,
                "A_x": ax,
                "A_y": ay,
                "A_z": az,
                "B_x": bx,
                "B_y": by,
                "B_z": bz,
                "signed_cross_xy": cross,
                "axis_line_distance_A": distance,
                "radial_angle_deg": angle,
                "axis_fraction_on_segment": t,
                "axis_between_atoms": 0.0 <= t <= 1.0,
            }
        )
    if not rows:
        raise OpposingPhosphateError(
            f"{pdb_path} has no opposing phosphate pairs; the helix needs at least 2 base pairs."
        )
    return rows


def summarize_rows(rows: Sequence[Dict[str, object]]) -> Dict[str, float]:
    """Collapse per-pair metrics into the summary used to score one X-disp."""
    n_rows = len(rows)
    if not n_rows:
        raise OpposingPhosphateError("No opposing phosphate pairs to summarize.")
    crosses = [float(row["signed_cross_xy"]) for row in rows]
    distances = [float(row["axis_line_distance_A"]) for row in rows]
    angles = [float(row["radial_angle_deg"]) for row in rows]
    fractions = [float(row["axis_fraction_on_segment"]) for row in rows]
    return {
        "mean_signed_cross_xy": sum(crosses) / n_rows,
        "mean_abs_cross_xy": sum(abs(value) for value in crosses) / n_rows,
        "max_abs_cross_xy": max(abs(value) for value in crosses),
        "rms_axis_line_distance_A": math.sqrt(sum(value * value for value in distances) / n_rows),
        "mean_axis_line_distance_A": sum(distances) / n_rows,
        "max_axis_line_distance_A": max(distances),
        "mean_radial_angle_deg": sum(angles) / n_rows,
        "min_radial_angle_deg": min(angles),
        "max_radial_angle_deg": max(angles),
        "min_axis_fraction_on_segment": min(fractions),
        "max_axis_fraction_on_segment": max(fractions),
    }


def _score(record: Dict[str, object]) -> Tuple[float, float]:
    """Ranking key: smallest axis-line distance, then smallest residual cross."""
    return (
        float(record["rms_axis_line_distance_A"]),
        abs(float(record["mean_signed_cross_xy"])),
    )


# ----------------------------------------------------------------------
# Model building
# ----------------------------------------------------------------------
def _load_build_bdna() -> Callable[..., Dict[str, object]]:
    """Import the B-DNA builder on first use.

    Keeping it out of module import lets the geometry helpers and the search
    driver be imported, and tested with an injected ``measure``, without pulling
    in the DSSR and Phenix build stack.
    """
    try:
        from build_bdna import build_bdna
    except ImportError:  # imported as part of the bnp_na_lib package
        from bnp_na_lib.build_bdna import build_bdna
    return build_bdna


def build_and_measure(
    *,
    seq: str,
    x_disp: float,
    output_dir: Union[str, Path],
    run_phenix: bool,
    params_file: Optional[Union[str, Path]],
    fixed_overrides: Optional[Dict[str, float]] = None,
) -> Tuple[Dict[str, float], List[Dict[str, object]], Dict[str, object]]:
    """Build one B-DNA model at ``x_disp`` and measure its opposing phosphates."""
    build_bdna = _load_build_bdna()
    x_value = quantize_x(x_disp)
    mode = MODE_PHENIX if run_phenix else MODE_DSSR
    param_overrides = dict(fixed_overrides or {})
    param_overrides["X-disp"] = x_value
    result = build_bdna(
        seq,
        f"{seq}_xdisp_{x_label(x_value)}_{mode}",
        params_file if run_phenix else None,
        Path(output_dir),
        param_overrides=param_overrides,
        run_phenix=run_phenix,
    )
    rows = phosphate_pair_rows(Path(str(result["pdb_aligned"])))
    summary = summarize_rows(rows)
    return summary, rows, result


class _Evaluator:
    """Cache X-disp evaluations, report progress, and honour cancellation."""

    def __init__(
        self,
        *,
        seq: str,
        params_file: Optional[Union[str, Path]],
        fixed_overrides: Dict[str, float],
        measure: Optional[Callable[[float, bool], Dict[str, float]]] = None,
        progress: Optional[Callable[[str], None]] = None,
        should_cancel: Optional[Callable[[], bool]] = None,
    ) -> None:
        self.seq = seq
        self.params_file = params_file
        self.fixed_overrides = fixed_overrides
        self.measure = measure
        self.progress = progress
        self.should_cancel = should_cancel
        self.cache: Dict[Tuple[str, float], Dict[str, object]] = {}
        self.records: List[Dict[str, object]] = []
        self.build_counts: Dict[str, int] = {MODE_DSSR: 0, MODE_PHENIX: 0}

    def emit(self, message: str) -> None:
        if self.progress is not None:
            self.progress(message)

    def _check_cancel(self) -> None:
        if self.should_cancel is not None and self.should_cancel():
            raise OpposingPhosphateCancelled("Search cancelled.")

    def __call__(self, x_disp: float, *, run_phenix: bool, stage: str) -> Dict[str, object]:
        x_value = quantize_x(x_disp)
        mode = MODE_PHENIX if run_phenix else MODE_DSSR
        key = (mode, x_value)
        cached = self.cache.get(key)
        if cached is not None:
            record = dict(cached)
            record["stage"] = stage
            return record

        self._check_cancel()
        if self.measure is not None:
            summary = self.measure(x_value, run_phenix)
        else:
            with tempfile.TemporaryDirectory(prefix=f"oppo_phos_{mode}_{x_label(x_value)}_") as tmp:
                summary, _rows, _result = build_and_measure(
                    seq=self.seq,
                    x_disp=x_value,
                    output_dir=Path(tmp),
                    run_phenix=run_phenix,
                    params_file=self.params_file,
                    fixed_overrides=self.fixed_overrides,
                )
        self.build_counts[mode] += 1

        record: Dict[str, object] = dict(summary)
        record.update(
            {
                "mode": mode,
                "stage": stage,
                "selected": False,
                "x_disp": x_value,
                "fixed_parameter_overrides": overrides_text(self.fixed_overrides),
            }
        )
        self.cache[key] = dict(record)
        self.records.append(record)
        self.emit(
            f"  [{mode}] X-disp {x_value:.4f} -> mean signed cross {float(record['mean_signed_cross_xy']):+.6f}, "
            f"RMS axis-line distance {float(record['rms_axis_line_distance_A']):.6f} A"
        )
        return record

    def best(self, mode: str) -> Dict[str, object]:
        candidates = [row for row in self.records if row["mode"] == mode]
        if not candidates:
            raise OpposingPhosphateError(f"No {mode} evaluations were recorded.")
        return min(candidates, key=_score)


# ----------------------------------------------------------------------
# Search
# ----------------------------------------------------------------------
def _scan_values(x_min: float, x_max: float, x_step: float) -> List[float]:
    if x_step <= 0:
        raise OpposingPhosphateError("X-disp scan step must be positive.")
    if x_max <= x_min:
        raise OpposingPhosphateError("X-disp scan max must be larger than min.")
    values: List[float] = []
    current = quantize_x(x_min)
    limit = quantize_x(x_max)
    while current <= limit + 1e-9:
        values.append(current)
        current = quantize_x(current + x_step)
    if values[-1] < limit:
        values.append(limit)
    return values


def _grid_sweep(
    evaluator: _Evaluator,
    *,
    center: float,
    halfwidth: float,
    run_phenix: bool,
    stage: str,
) -> None:
    """Evaluate every four-decimal X-disp within ``halfwidth`` of ``center``."""
    steps = int(round(halfwidth / X_DISP_QUANTUM))
    start = int(round(quantize_x(center) / X_DISP_QUANTUM))
    for offset in range(-steps, steps + 1):
        evaluator((start + offset) * X_DISP_QUANTUM, run_phenix=run_phenix, stage=stage)


def _cross_slope(records: Sequence[Dict[str, object]]) -> Optional[float]:
    """Least-squares d(mean signed cross)/d(X-disp) over the supplied records."""
    points = [(float(row["x_disp"]), float(row["mean_signed_cross_xy"])) for row in records]
    unique = {x: y for x, y in points}
    if len(unique) < 2:
        return None
    xs = list(unique)
    ys = [unique[x] for x in xs]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denominator = sum((x - mean_x) ** 2 for x in xs)
    if denominator <= 0:
        return None
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denominator
    if abs(slope) < 1e-9:
        return None
    return slope


def _dssr_search(
    evaluator: _Evaluator,
    *,
    x_min: float,
    x_max: float,
    x_step: float,
    grid_halfwidth: float,
) -> Dict[str, object]:
    evaluator.emit("Raw DSSR coarse scan...")
    bracket: Optional[Tuple[Dict[str, object], Dict[str, object]]] = None
    previous: Optional[Dict[str, object]] = None
    for x_value in _scan_values(x_min, x_max, x_step):
        record = evaluator(x_value, run_phenix=False, stage="coarse_grid")
        if previous is not None:
            f0 = float(previous["mean_signed_cross_xy"])
            f1 = float(record["mean_signed_cross_xy"])
            if f0 == 0.0 or f0 * f1 <= 0.0:
                bracket = (previous, record)
                break
        previous = record

    if bracket is None:
        best = evaluator.best(MODE_DSSR)
        raise OpposingPhosphateError(
            "The raw DSSR X-disp scan never crossed the target. "
            f"The closest scanned value was X-disp={float(best['x_disp']):.4f} A with "
            f"RMS axis-line distance={float(best['rms_axis_line_distance_A']):.6f} A. "
            "Widen the scan range and try again."
        )

    lo = float(bracket[0]["x_disp"])
    hi = float(bracket[1]["x_disp"])
    flo = float(bracket[0]["mean_signed_cross_xy"])
    evaluator.emit(f"Bracketed between X-disp {lo:.4f} and {hi:.4f}; bisecting...")
    for iteration in range(BISECTION_STEPS):
        mid = quantize_x((lo + hi) / 2.0)
        record = evaluator(mid, run_phenix=False, stage=f"bisection_{iteration:02d}")
        fm = float(record["mean_signed_cross_xy"])
        if flo * fm > 0:
            lo = mid
            flo = fm
        else:
            hi = mid

    center = quantize_x((lo + hi) / 2.0)
    evaluator.emit(f"Sweeping four-decimal values within {grid_halfwidth:.4f} A of X-disp {center:.4f}...")
    _grid_sweep(
        evaluator,
        center=center,
        halfwidth=grid_halfwidth,
        run_phenix=False,
        stage="four_decimal_grid",
    )
    return evaluator.best(MODE_DSSR)


def _phenix_search(
    evaluator: _Evaluator,
    *,
    start_x: float,
    slope: Optional[float],
    grid_halfwidth: float,
    max_builds: int,
) -> Dict[str, object]:
    """Refine the root for the Phenix + phosphate-regularization workflow.

    Minimization shifts the root by roughly 0.3 A from the raw DSSR answer, far
    outside any narrow window around it, so this seeds a Newton step with the
    slope measured during the raw sweep and then runs a secant iteration. The
    response is close to linear in X-disp, so this normally reaches the
    four-decimal root in about three builds.
    """
    evaluator.emit("Refining with Phenix minimization + phosphate regularization...")

    def evaluate(x_value: float, stage: str) -> Tuple[float, float]:
        if evaluator.build_counts[MODE_PHENIX] >= max_builds:
            raise OpposingPhosphateError(
                f"Stopped after {max_builds} Phenix builds without converging. "
                "Raise the Phenix build cap or check the fixed parameter values."
            )
        record = evaluator(x_value, run_phenix=True, stage=stage)
        return float(record["x_disp"]), float(record["mean_signed_cross_xy"])

    x0, f0 = evaluate(start_x, "phenix_seed")
    if f0 == 0.0:
        center = x0
    else:
        if slope is None:
            raise OpposingPhosphateError(
                "Could not estimate the X-disp response slope from the raw DSSR sweep."
            )
        x1, f1 = evaluate(max(0.0, x0 - f0 / slope), "phenix_newton")

        for iteration in range(SECANT_STEPS):
            if f1 == f0 or abs(x1 - x0) < X_DISP_QUANTUM / 2.0:
                break
            step = f1 * (x1 - x0) / (f1 - f0)
            nxt = quantize_x(max(0.0, x1 - step))
            if nxt == x1:
                break
            x0, f0 = x1, f1
            x1, f1 = evaluate(nxt, f"phenix_secant_{iteration:02d}")
        center = x1

    # The sweep costs one build per point, so trim it to whatever the cap still
    # allows. Cached points from the secant phase are free and often land inside
    # the window, but budgeting for the worst case keeps the cap honest.
    wanted_steps = int(round(grid_halfwidth / X_DISP_QUANTUM))
    remaining = max_builds - evaluator.build_counts[MODE_PHENIX]
    allowed_steps = max(0, min(wanted_steps, (remaining - 1) // 2))
    if allowed_steps < wanted_steps:
        evaluator.emit(
            f"  Phenix build cap of {max_builds} reached; narrowing the final sweep to "
            f"{allowed_steps * X_DISP_QUANTUM:.4f} A."
        )
    if remaining > 0:
        evaluator.emit(
            f"Sweeping four-decimal values within {allowed_steps * X_DISP_QUANTUM:.4f} A of X-disp {center:.4f}..."
        )
        _grid_sweep(
            evaluator,
            center=center,
            halfwidth=allowed_steps * X_DISP_QUANTUM,
            run_phenix=True,
            stage="four_decimal_grid",
        )
    return evaluator.best(MODE_PHENIX)


def find_opposing_phosphate_xdisp(
    seq: str,
    *,
    params_file: Optional[Union[str, Path]] = None,
    fixed_overrides: Optional[Dict[str, float]] = None,
    refine_with_phenix: bool = False,
    dssr_x_min: float = DEFAULT_DSSR_X_MIN,
    dssr_x_max: float = DEFAULT_DSSR_X_MAX,
    dssr_x_step: float = DEFAULT_DSSR_X_STEP,
    dssr_grid_halfwidth: float = DEFAULT_DSSR_GRID_HALFWIDTH,
    phenix_grid_halfwidth: float = DEFAULT_PHENIX_GRID_HALFWIDTH,
    max_phenix_builds: int = DEFAULT_MAX_PHENIX_BUILDS,
    progress: Optional[Callable[[str], None]] = None,
    should_cancel: Optional[Callable[[], bool]] = None,
    measure: Optional[Callable[[float, bool], Dict[str, float]]] = None,
) -> Dict[str, object]:
    """Search for the X-disp that puts opposing phosphates across the helix axis.

    The raw DSSR search always runs; it is fast and it supplies both the seed
    and the response slope used by the optional Phenix refinement. Pass
    ``measure`` to score X-disp values without building real models.
    """
    seq = (seq or "").strip()
    if not seq:
        raise OpposingPhosphateError("Sequence is blank.")
    clean_overrides = validate_fixed_overrides(fixed_overrides)
    if refine_with_phenix and measure is None and not params_file:
        raise OpposingPhosphateError(
            "A params file is required for the Phenix + phosphate-regularization refinement."
        )

    evaluator = _Evaluator(
        seq=seq,
        params_file=params_file,
        fixed_overrides=clean_overrides,
        measure=measure,
        progress=progress,
        should_cancel=should_cancel,
    )
    evaluator.emit(f"Sequence: {seq}")
    evaluator.emit(f"Fixed non-X parameters: {overrides_text(clean_overrides)}")

    dssr_best = _dssr_search(
        evaluator,
        x_min=dssr_x_min,
        x_max=dssr_x_max,
        x_step=dssr_x_step,
        grid_halfwidth=dssr_grid_halfwidth,
    )
    dssr_x = float(dssr_best["x_disp"])
    evaluator.emit(f"Raw DSSR X-disp = {dssr_x:.4f} A")

    phenix_best: Optional[Dict[str, object]] = None
    if refine_with_phenix:
        grid_records = [
            row
            for row in evaluator.records
            if row["mode"] == MODE_DSSR and row["stage"] == "four_decimal_grid"
        ]
        phenix_best = _phenix_search(
            evaluator,
            start_x=dssr_x,
            slope=_cross_slope(grid_records),
            grid_halfwidth=phenix_grid_halfwidth,
            max_builds=max_phenix_builds,
        )
        evaluator.emit(f"Phenix + phosphate-regularized X-disp = {float(phenix_best['x_disp']):.4f} A")

    recommended = phenix_best if phenix_best is not None else dssr_best
    for row in evaluator.records:
        row["selected"] = row is dssr_best or (phenix_best is not None and row is phenix_best)

    return {
        "seq": seq,
        "fixed_overrides": clean_overrides,
        "dssr": dict(dssr_best),
        "phenix": dict(phenix_best) if phenix_best is not None else None,
        "recommended_mode": str(recommended["mode"]),
        "recommended_x_disp": float(recommended["x_disp"]),
        "records": [dict(row) for row in evaluator.records],
        "build_counts": dict(evaluator.build_counts),
    }


def format_result_report(result: Dict[str, object]) -> str:
    """Human-readable summary of a completed search."""
    lines = [
        "=== Opposing phosphate X-disp search ===",
        f"Sequence                : {result['seq']}",
        f"Fixed non-X parameters  : {overrides_text(result.get('fixed_overrides'))}",  # type: ignore[arg-type]
    ]

    def block(title: str, record: Optional[Dict[str, object]]) -> None:
        if not record:
            return
        lines.extend(
            [
                "",
                title,
                f"  X-disp                     : {float(record['x_disp']):.4f} A",
                f"  RMS axis-line distance     : {float(record['rms_axis_line_distance_A']):.6f} A",
                f"  Mean around-axis angle     : {float(record['mean_radial_angle_deg']):.6f} deg",
                f"  Around-axis angle range    : {float(record['min_radial_angle_deg']):.6f}"
                f" to {float(record['max_radial_angle_deg']):.6f} deg",
                f"  Mean signed xy cross       : {float(record['mean_signed_cross_xy']):+.6f}",
            ]
        )

    block("Raw DSSR rebuild only:", result.get("dssr"))  # type: ignore[arg-type]
    block("Phenix minimization + phosphate regularization:", result.get("phenix"))  # type: ignore[arg-type]

    counts = result.get("build_counts") or {}
    lines.extend(
        [
            "",
            f"Models built            : {counts.get(MODE_DSSR, 0)} raw DSSR, {counts.get(MODE_PHENIX, 0)} Phenix",
            f"Recommended X-disp      : {float(result['recommended_x_disp']):.4f} A "
            f"({'Phenix + regularized' if result['recommended_mode'] == MODE_PHENIX else 'raw DSSR'})",
        ]
    )
    return "\n".join(lines)


def _parse_override_items(items: Sequence[str]) -> Dict[str, float]:
    raw: Dict[str, object] = {}
    for item in items:
        if "=" not in item:
            raise OpposingPhosphateError(f"Parameter override must use key=value form: {item!r}")
        key, raw_value = item.split("=", 1)
        raw[key.strip()] = raw_value.strip()
    return validate_fixed_overrides(raw)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find the B-DNA X-disp that places opposing phosphate P atoms across the helix axis."
        )
    )
    parser.add_argument("--seq", default="A31", help="Compact bnp_na sequence. Default: A31.")
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Fixed non-X DSSR parameter override, repeatable, for example --param h-Twist=30. "
            "Do not use this for X-disp."
        ),
    )
    parser.add_argument(
        "--phenix",
        action="store_true",
        help="Also refine with Phenix minimization plus phosphate regularization. This is much slower.",
    )
    parser.add_argument(
        "--params-file",
        default=None,
        help="Params file for phenix.geometry_minimization. Default: bnp_na_lib/min_P_C5.params.",
    )
    parser.add_argument("--dssr-x-min", type=float, default=DEFAULT_DSSR_X_MIN, help="Coarse scan minimum X-disp.")
    parser.add_argument("--dssr-x-max", type=float, default=DEFAULT_DSSR_X_MAX, help="Coarse scan maximum X-disp.")
    parser.add_argument("--dssr-x-step", type=float, default=DEFAULT_DSSR_X_STEP, help="Coarse scan X-disp step.")
    parser.add_argument(
        "--max-phenix-builds",
        type=int,
        default=DEFAULT_MAX_PHENIX_BUILDS,
        help="Safety cap on Phenix builds. Default: %(default)s.",
    )
    parser.add_argument("--quiet", action="store_true", help="Only print the final report.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    params_file = args.params_file
    if params_file is None:
        params_file = Path(__file__).resolve().parent / "min_P_C5.params"
    try:
        result = find_opposing_phosphate_xdisp(
            args.seq,
            params_file=params_file,
            fixed_overrides=_parse_override_items(args.param),
            refine_with_phenix=args.phenix,
            dssr_x_min=args.dssr_x_min,
            dssr_x_max=args.dssr_x_max,
            dssr_x_step=args.dssr_x_step,
            max_phenix_builds=args.max_phenix_builds,
            progress=None if args.quiet else lambda message: print(message, flush=True),
        )
    except OpposingPhosphateError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(format_result_report(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
