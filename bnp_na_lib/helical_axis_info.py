#!/usr/bin/env python3
"""Get DSSR helical-axis information for selected PDB chains."""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple, Union

from align2z import command_to_text, parse_first_axis_points, run_dssr_more_axis


Vector3 = Tuple[float, float, float]


class HelicalAxisInfoError(Exception):
    """Raised when selected-chain DSSR axis extraction fails."""


@dataclass
class HelicalAxisInfo:
    input_pdb: Path
    chains: Tuple[str, str]
    selected_pdb: Path
    dssr_output: Path
    dssr_command: str
    dssr_log: str
    selected_atom_counts: Dict[str, int]
    start_point: Vector3
    end_point: Vector3
    axis_vector: Vector3
    axis_length: float
    unit_vector: Vector3
    reference_vector: Optional[Vector3] = None
    reference_unit_vector: Optional[Vector3] = None
    angle_degrees: Optional[float] = None
    helix_bp_count: Optional[int] = None
    estimated_full_helix_length: Optional[float] = None
    bild_output: Optional[Path] = None


def _fmt(value: float) -> str:
    return f"{float(value):.6g}"


def _point_text(point: Vector3) -> str:
    return " ".join(_fmt(value) for value in point)


def parse_chain_ids(text: str) -> Tuple[str, str]:
    """Parse two one-character PDB chain IDs from comma/space text."""
    cleaned = text.replace(",", " ").replace(";", " ").strip()
    parts = cleaned.split()
    if len(parts) == 1 and len(parts[0]) == 2:
        parts = [parts[0][0], parts[0][1]]
    if len(parts) != 2:
        raise ValueError("Provide exactly two PDB chain IDs, for example: C D")
    chain_a, chain_b = parts[0], parts[1]
    if len(chain_a) != 1 or len(chain_b) != 1:
        raise ValueError("PDB chain IDs must be one character each.")
    if chain_a == chain_b:
        raise ValueError("The two chain IDs must be different.")
    return chain_a, chain_b


def parse_vector(text: str) -> Vector3:
    """Parse an xyz vector from comma- or whitespace-separated text."""
    parts = text.replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("Vector must have three numbers, for example: 0 0 1")
    try:
        return float(parts[0]), float(parts[1]), float(parts[2])
    except Exception as exc:
        raise ValueError(f"Could not parse vector values: {text!r}") from exc


def _sub(a: Vector3, b: Vector3) -> Vector3:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def _add(a: Vector3, b: Vector3) -> Vector3:
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def _scale(v: Vector3, factor: float) -> Vector3:
    return v[0] * factor, v[1] * factor, v[2] * factor


def _dot(a: Vector3, b: Vector3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(v: Vector3) -> float:
    return math.sqrt(_dot(v, v))


def _unit(v: Vector3, label: str) -> Vector3:
    length = _norm(v)
    if length <= 1e-10:
        raise HelicalAxisInfoError(f"{label} has near-zero length.")
    return v[0] / length, v[1] / length, v[2] / length


def _angle_degrees(a_unit: Vector3, b_unit: Vector3) -> float:
    dot = max(-1.0, min(1.0, _dot(a_unit, b_unit)))
    return math.degrees(math.acos(dot))


def estimate_full_helix_length(axis_distance: float, helix_bp_count: int) -> float:
    """Estimate full helix length from DSSR endpoint distance and base-pair count."""
    bp_count = int(helix_bp_count)
    if bp_count <= 1:
        raise HelicalAxisInfoError("Helix length in bp must be greater than 1 for length estimation.")
    return float(axis_distance) / float(bp_count - 1) * float(bp_count)


def available_chain_atom_counts(pdb_path: Union[str, Path]) -> Dict[str, int]:
    """Return ATOM/HETATM counts by chain ID."""
    counts: Dict[str, int] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain_id = line[21] if len(line) > 21 else " "
            counts[chain_id] = counts.get(chain_id, 0) + 1
    return counts


def _chain_label(chain_id: str) -> str:
    return chain_id if chain_id.strip() else "<blank>"


def _chain_file_label(chains: Tuple[str, str]) -> str:
    safe = []
    for chain_id in chains:
        if chain_id.strip():
            safe.append("".join(ch for ch in chain_id if ch.isalnum() or ch in ("_", "-")) or "chain")
        else:
            safe.append("blank")
    return "".join(safe)


def _canonical_chain_pair(chains: Tuple[str, str]) -> Tuple[str, str]:
    return tuple(sorted(chains))  # type: ignore[return-value]


def _reverse_axis_for_chain_order(chains: Tuple[str, str]) -> bool:
    """Use chain order to assign a deterministic axis direction."""
    return chains != _canonical_chain_pair(chains)


def write_selected_chains_pdb(
    input_pdb: Union[str, Path],
    chains: Tuple[str, str],
    out_pdb: Union[str, Path],
) -> Tuple[Path, Dict[str, int]]:
    """Write a temporary PDB containing only the selected chains."""
    input_path = Path(input_pdb).expanduser()
    output_path = Path(out_pdb).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    selected = set(chains)
    lines_by_chain = {chain_id: [] for chain_id in chains}
    counts = {chain_id: 0 for chain_id in chains}
    with open(input_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21] if len(line) > 21 else " "
                if chain_id in selected:
                    lines_by_chain[chain_id].append(line if line.endswith("\n") else line + "\n")
                    counts[chain_id] += 1
            elif line.startswith("TER"):
                chain_id = line[21] if len(line) > 21 else " "
                if chain_id in selected:
                    lines_by_chain[chain_id].append(line if line.endswith("\n") else line + "\n")

    missing = [chain_id for chain_id, count in counts.items() if count == 0]
    if missing:
        available = available_chain_atom_counts(input_path)
        available_text = ", ".join(
            f"{_chain_label(chain_id)}({count})" for chain_id, count in sorted(available.items())
        )
        if not available_text:
            available_text = "none"
        missing_text = ", ".join(_chain_label(chain_id) for chain_id in missing)
        raise HelicalAxisInfoError(
            f"No atoms found for chain ID(s): {missing_text}. Available chains: {available_text}."
        )

    selected_lines = []
    for chain_id in chains:
        selected_lines.extend(lines_by_chain[chain_id])
    selected_lines.append("END\n")
    output_path.write_text("".join(selected_lines), encoding="utf-8")
    return output_path, counts


def default_output_paths(
    input_pdb: Union[str, Path],
    chains: Tuple[str, str],
    workdir: Union[str, Path],
) -> Tuple[Path, Path]:
    """Return selected-chain PDB and DSSR .out paths."""
    input_path = Path(input_pdb).expanduser()
    workdir_path = Path(workdir).expanduser()
    label = _chain_file_label(chains)
    selected_pdb = workdir_path / f"{input_path.stem}_chains_{label}_dssr_input.pdb"
    dssr_out = workdir_path / f"{input_path.stem}_chains_{label}_dssr_more.out"
    return selected_pdb, dssr_out


def write_axis_bild(
    out_path: Union[str, Path],
    info: HelicalAxisInfo,
    *,
    axis_radius: float = 1.0,
    sphere_radius: float = 1.25,
    vector_radius: float = 0.7,
    draw_reference_vector: bool = True,
    reference_vector_length: Optional[float] = None,
) -> Path:
    """Write a Chimera/ChimeraX BILD drawing for the DSSR helical axis."""
    out = Path(out_path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    axis_head = max(axis_radius * 2.5, 1e-6)
    vector_head = max(vector_radius * 2.5, 1e-6)
    lines = [
        ".comment Generated by helical_axis_info.py",
        f".comment Input PDB: {info.input_pdb}",
        f".comment Selected chains: {info.chains[0]} {info.chains[1]}",
        f".comment Axis start point: {_point_text(info.start_point)}",
        f".comment Axis end point: {_point_text(info.end_point)}",
        f".comment Axis unit vector: {_point_text(info.unit_vector)}",
    ]
    if info.reference_vector is not None and info.angle_degrees is not None:
        lines.extend(
            [
                f".comment Reference vector: {_point_text(info.reference_vector)}",
                f".comment Angle to reference vector, degrees: {_fmt(info.angle_degrees)}",
            ]
        )
    lines.extend(
        [
            ".comment DSSR helical axis",
            ".color 0 0.65 1",
            f".arrow {_point_text(info.start_point)} {_point_text(info.end_point)} "
            f"{_fmt(axis_radius)} {_fmt(axis_head)} 0.85",
            ".comment Axis start/end markers",
            ".color 0 1 0",
            f".sphere {_point_text(info.start_point)} {_fmt(sphere_radius)}",
            ".color 1 0 0",
            f".sphere {_point_text(info.end_point)} {_fmt(sphere_radius)}",
        ]
    )
    if info.reference_unit_vector is not None and draw_reference_vector:
        ref_length = info.axis_length if reference_vector_length is None else float(reference_vector_length)
        if ref_length <= 0:
            raise HelicalAxisInfoError("Reference-vector BILD length must be positive.")
        reference_end = _add(info.start_point, _scale(info.reference_unit_vector, ref_length))
        lines.extend(
            [
                ".comment Reference vector drawn from the DSSR axis start point",
                ".color 1 0.55 0",
                f".arrow {_point_text(info.start_point)} {_point_text(reference_end)} "
                f"{_fmt(vector_radius)} {_fmt(vector_head)} 0.85",
            ]
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def get_helical_axis_info(
    input_pdb: Union[str, Path],
    chains: Union[str, Sequence[str]],
    *,
    reference_vector: Optional[Union[str, Sequence[float]]] = None,
    workdir: Optional[Union[str, Path]] = None,
    bild_output: Optional[Union[str, Path]] = None,
    draw_reference_vector: bool = True,
    reference_vector_length: Optional[float] = None,
    helix_bp_count: Optional[int] = None,
) -> HelicalAxisInfo:
    """Run DSSR --more on selected chains and return axis metadata."""
    input_path = Path(input_pdb).expanduser()
    if not input_path.is_absolute():
        input_path = input_path.resolve()
    if not input_path.exists():
        raise HelicalAxisInfoError(f"Input PDB not found: {input_path}")

    if isinstance(chains, str):
        chain_pair = parse_chain_ids(chains)
    else:
        chain_list = list(chains)
        if len(chain_list) != 2:
            raise ValueError("Provide exactly two PDB chain IDs.")
        chain_pair = parse_chain_ids(" ".join(str(chain_id) for chain_id in chain_list))

    if reference_vector is None:
        ref_vec = None
    elif isinstance(reference_vector, str):
        ref_vec = parse_vector(reference_vector)
    else:
        ref_values = tuple(float(value) for value in reference_vector)
        if len(ref_values) != 3:
            raise ValueError("Reference vector must have three numbers.")
        ref_vec = ref_values  # type: ignore[assignment]

    workdir_path = Path(workdir).expanduser() if workdir is not None else input_path.parent
    if not workdir_path.is_absolute():
        workdir_path = workdir_path.resolve()
    workdir_path.mkdir(parents=True, exist_ok=True)

    dssr_chain_pair = _canonical_chain_pair(chain_pair)
    selected_pdb, dssr_out = default_output_paths(input_path, chain_pair, workdir_path)
    selected_pdb, dssr_counts = write_selected_chains_pdb(input_path, dssr_chain_pair, selected_pdb)
    counts = {chain_id: dssr_counts[chain_id] for chain_id in chain_pair}

    ok, dssr_log, cmd = run_dssr_more_axis(selected_pdb, dssr_out, cwd=workdir_path)
    if not ok or not dssr_out.exists():
        raise HelicalAxisInfoError(
            "DSSR --more failed while extracting the selected-chain helical axis.\n"
            f"Command: {command_to_text(cmd)}\n"
            f"Expected output file: {dssr_out}\n"
            f"Output:\n{dssr_log}"
        )

    p1_raw, p2_raw = parse_first_axis_points(dssr_out)
    start = float(p1_raw[0]), float(p1_raw[1]), float(p1_raw[2])
    end = float(p2_raw[0]), float(p2_raw[1]), float(p2_raw[2])
    if _reverse_axis_for_chain_order(chain_pair):
        start, end = end, start
    axis_vector = _sub(end, start)
    axis_length = _norm(axis_vector)
    unit = _unit(axis_vector, "DSSR helical axis")

    if ref_vec is None:
        ref_unit = None
        angle = None
    else:
        ref_unit = _unit(ref_vec, "Reference vector")
        angle = _angle_degrees(unit, ref_unit)

    if helix_bp_count is None:
        estimated_full_helix_length = None
    else:
        estimated_full_helix_length = estimate_full_helix_length(axis_length, int(helix_bp_count))

    info = HelicalAxisInfo(
        input_pdb=input_path,
        chains=chain_pair,
        selected_pdb=selected_pdb,
        dssr_output=dssr_out,
        dssr_command=command_to_text(cmd),
        dssr_log=dssr_log,
        selected_atom_counts=counts,
        start_point=start,
        end_point=end,
        axis_vector=axis_vector,
        axis_length=axis_length,
        unit_vector=unit,
        reference_vector=ref_vec,
        reference_unit_vector=ref_unit,
        angle_degrees=angle,
        helix_bp_count=helix_bp_count,
        estimated_full_helix_length=estimated_full_helix_length,
    )
    if bild_output:
        info.bild_output = write_axis_bild(
            bild_output,
            info,
            draw_reference_vector=draw_reference_vector,
            reference_vector_length=reference_vector_length,
        )
    return info


def format_axis_info_report(info: HelicalAxisInfo) -> str:
    """Return a readable report for GUI logs and CLI output."""
    chain_counts = ", ".join(
        f"{_chain_label(chain_id)}={count}" for chain_id, count in info.selected_atom_counts.items()
    )
    lines = [
        "=== DSSR helical-axis info ===",
        f"Input PDB: {info.input_pdb}",
        f"Selected chains: {info.chains[0]} {info.chains[1]}",
        f"Selected-chain PDB: {info.selected_pdb}",
        f"Selected atom counts: {chain_counts}",
        f"DSSR --more output: {info.dssr_output}",
        f"DSSR command: {info.dssr_command}",
        f"Axis start point (A): {_point_text(info.start_point)}",
        f"Axis end point (A): {_point_text(info.end_point)}",
        f"Axis vector (A): {_point_text(info.axis_vector)}",
        f"Start-to-end distance (A): {_fmt(info.axis_length)}",
        "Distance note: DSSR point-one/point-two span the start and end base-pair centers, "
        "so this distance is approximately one base-pair step shorter than the full helix length.",
        f"Axis unit vector: {_point_text(info.unit_vector)}",
    ]
    if info.helix_bp_count is not None and info.estimated_full_helix_length is not None:
        lines.extend(
            [
                f"Helix length input (bp): {info.helix_bp_count}",
                f"Estimated full helix length (A): {_fmt(info.estimated_full_helix_length)}",
                "Estimate formula: start-to-end distance / (bp - 1) * bp",
            ]
        )
    if info.reference_vector is not None and info.reference_unit_vector is not None:
        lines.extend(
            [
                f"Reference vector: {_point_text(info.reference_vector)}",
                f"Reference unit vector: {_point_text(info.reference_unit_vector)}",
                f"Angle between axis and reference vector (deg): {_fmt(info.angle_degrees or 0.0)}",
            ]
        )
    else:
        lines.append("Reference vector: not provided")
    if info.bild_output is not None:
        lines.append(f"BILD output: {info.bild_output}")
    if info.dssr_log:
        lines.append(f"DSSR output:\n{info.dssr_log}")
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Get DSSR helical-axis information for two selected PDB chains.")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file.")
    parser.add_argument("--chains", required=True, help="Two chain IDs, for example: 'C D' or 'C,D'.")
    parser.add_argument(
        "--vector",
        default=None,
        help="Optional reference vector for angle reporting, for example: '0 0 1'.",
    )
    parser.add_argument("--workdir", default=None, help="Folder for selected-chain PDB and DSSR .out files.")
    parser.add_argument("--bild", default=None, help="Optional output Chimera/ChimeraX .bild file for the axis.")
    parser.add_argument(
        "--helix-bp",
        type=int,
        default=None,
        help="Optional helix length in base pairs for full-length estimate. Must be greater than 1.",
    )
    parser.add_argument(
        "--no-reference-vector-bild",
        action="store_true",
        help="Do not draw the reference vector in the optional .bild output.",
    )
    parser.add_argument(
        "--reference-vector-length",
        type=float,
        default=None,
        help="Reference-vector drawing length in Angstrom. Default: DSSR axis length.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    info = get_helical_axis_info(
        args.input,
        args.chains,
        reference_vector=args.vector,
        workdir=args.workdir,
        bild_output=args.bild,
        draw_reference_vector=not args.no_reference_vector_bild,
        reference_vector_length=args.reference_vector_length,
        helix_bp_count=args.helix_bp,
    )
    print(format_axis_info_report(info))


if __name__ == "__main__":
    main()
