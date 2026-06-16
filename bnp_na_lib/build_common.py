"""Common helpers for bnp_na build modules."""
from __future__ import annotations

import os
import re
import shutil
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union


class PipelineError(Exception):
    """Raised when a build pipeline fails. Carries accumulated log text."""

    def __init__(self, message: str, log_text: str = ""):
        super().__init__(message)
        self.log_text = log_text


PARAM_KEYS = [
    "Shear",
    "Stretch",
    "Stagger",
    "Buckle",
    "Propeller",
    "Opening",
    "X-disp",
    "Y-disp",
    "h-Rise",
    "Incl.",
    "Tip",
    "h-Twist",
]

HELICAL_STEP_KEYS = ["X-disp", "Y-disp", "h-Rise", "Incl.", "Tip", "h-Twist"]
SENTINEL_LAST6 = set(HELICAL_STEP_KEYS)

# Defaults supplied by Di. Tables are written with four digits after the decimal.
DEFAULT_PARAMS: Dict[str, List[float]] = {
    "B-DNA": [
        0.0,
        -0.15,
        0.09,
        0.5,
        -11.4,
        0.6,
        0.05,
        0.02,
        3.4,
        2.1,
        0.0,
        34.2857,
    ],
    "A-DNA": [
        0.00014375,
        -0.14475,
        0.06379375,
        0.0003375,
        -10.515775,
        -1.81695,
        -4.46163125,
        7.5e-05,
        2.5466,
        22.646025,
        9.375e-05,
        32.72727273,
    ],
    "A-RNA": [
        0.01373125,
        -0.0848125,
        0.0125875,
        -0.00438125,
        -2.07649375,
        -1.66756875,
        -4.05126875,
        0.0677625,
        2.8120125,
        15.51483125,
        0.7866125,
        32.72727273,
    ],
}

DNA_PAIR_TAGS = {"A": "A-T", "T": "T-A", "C": "C-G", "G": "G-C"}
RNA_PAIR_TAGS = {"A": "A-U", "U": "U-A", "C": "C-G", "G": "G-C"}


def sanitize_basename(name: str) -> str:
    if not name:
        return ""
    safe = []
    for ch in name:
        safe.append(ch if (ch.isalnum() or ch in "._-") else "_")
    value = "".join(safe).strip(" .")
    if value in {"", ".", ".."}:
        return ""
    return value.replace("/", "_").replace("\\", "_")


def strip_pdb_ext(name: str) -> str:
    base = (name or "").strip().strip(" .")
    if base.lower().endswith(".pdb"):
        base = base[:-4]
    return sanitize_basename(base)


def sequence_alphabet(na_type: str) -> str:
    if na_type == "A-RNA":
        return "AUCG"
    return "ATCG"


def expand_sequence(seq: str, alphabet: str = "ATCG") -> str:
    """Expand compact syntax such as A10T5C2G.

    Parameters
    ----------
    seq:
        Input sequence. Whitespace is ignored.
    alphabet:
        Allowed residue letters, e.g. ATCG for DNA or AUCG for RNA.
    """
    allowed = set(alphabet.upper())
    s = re.sub(r"\s+", "", (seq or "").upper())
    if not s:
        return ""
    if ":" in s:
        raise ValueError("Two-strand 'lead:follow' syntax is not supported yet.")
    if all(ch in allowed for ch in s):
        return s

    out: List[str] = []
    i = 0
    while i < len(s):
        base = s[i]
        if base not in allowed:
            allowed_text = "/".join(sorted(allowed))
            raise ValueError(f"Invalid character '{base}'. Only {allowed_text} allowed.")
        i += 1
        j = i
        while j < len(s) and s[j].isdigit():
            j += 1
        count = int(s[i:j]) if j > i else 1
        if count <= 0:
            raise ValueError("Compact sequence counts must be positive integers.")
        out.append(base * count)
        i = j
    return "".join(out)


def fmt_value(value: object) -> str:
    try:
        numeric = Decimal(str(value))
        if numeric == Decimal("999999"):
            return "999999"
        return str(numeric.quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP))
    except (InvalidOperation, ValueError):
        return str(value)


def _column_header() -> str:
    return "#bp      " + " ".join(f"{name:>10s}" for name in PARAM_KEYS)


def _format_table_row(pair_tag: str, values: Sequence[object]) -> str:
    return f"{pair_tag:<7s}" + " ".join(f"{fmt_value(value):>10s}" for value in values)


def get_default_params(na_type: str) -> List[float]:
    if na_type not in DEFAULT_PARAMS:
        raise ValueError(f"No default parameters defined for {na_type}.")
    return list(DEFAULT_PARAMS[na_type])


def pair_tag_for_base(base: str, na_type: str) -> str:
    if na_type == "A-RNA":
        mapping = RNA_PAIR_TAGS
    else:
        mapping = DNA_PAIR_TAGS
    try:
        return mapping[base]
    except KeyError as exc:
        raise ValueError(f"Invalid base {base!r} for {na_type}.") from exc


def write_helical_table(
    seq_in: str,
    out_path: Union[str, Path],
    na_type: str,
    param_overrides: Optional[Dict[str, float]] = None,
) -> Tuple[Path, str]:
    """Write a DSSR rebuild table using local base-pair helical parameters.

    The first six values are local base-pair parameters. The final six are local
    helical parameters, and are set to 999999 on the final base-pair row, as DSSR
    rebuild expects for inter-base-pair quantities.
    """
    alphabet = sequence_alphabet(na_type)
    seq = expand_sequence(seq_in, alphabet=alphabet)
    if not seq:
        raise ValueError("Empty sequence.")

    defaults = get_default_params(na_type)
    overrides = dict(param_overrides or {})
    bad_keys = set(overrides) - set(PARAM_KEYS)
    if bad_keys:
        raise ValueError(f"Unsupported parameter override key(s): {sorted(bad_keys)}")

    out = Path(out_path)
    lines = [f"# {len(seq)} (no. of base pairs)", _column_header()]
    for idx, base in enumerate(seq):
        values: List[object] = list(defaults)
        for key, value in overrides.items():
            values[PARAM_KEYS.index(key)] = float(value)
        if idx == len(seq) - 1:
            for key in HELICAL_STEP_KEYS:
                values[PARAM_KEYS.index(key)] = 999999
        lines.append(_format_table_row(pair_tag_for_base(base, na_type), values))

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out, seq


def which_or(path_guess: str) -> Optional[str]:
    hit = shutil.which(path_guess)
    if hit:
        return hit
    if Path(path_guess).is_file():
        return path_guess
    return None


def _path_arg(path: Union[str, Path], cwd: Union[str, Path]) -> str:
    p = Path(path)
    cwd_path = Path(cwd).resolve()
    try:
        if p.resolve().parent == cwd_path:
            return p.name
    except Exception:
        pass
    return str(p)


def run_dssr_rebuild(
    par_file: Union[str, Path],
    pdb_out: Union[str, Path],
    cwd: Union[str, Path],
    backbone: str,
) -> Tuple[bool, str, List[str]]:
    exe = which_or("x3dna-dssr") or "/usr/local/bin/x3dna-dssr"
    if not which_or(exe):
        return False, "x3dna-dssr not found in PATH or at /usr/local/bin/x3dna-dssr", []

    cmd = [
        exe,
        "rebuild",
        f"--backbone={backbone}",
        "--par-type=heli",
        f"--par-file={_path_arg(par_file, cwd)}",
        f"-o={_path_arg(pdb_out, cwd)}",
    ]
    try:
        cp = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
        output = (cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")
        return cp.returncode == 0, output.strip(), cmd
    except Exception as exc:
        return False, str(exc), cmd


def run_dssr_fiber_z(repeat: int, pdb_out: Union[str, Path], cwd: Union[str, Path]) -> Tuple[bool, str, List[str]]:
    exe = which_or("x3dna-dssr") or "/usr/local/bin/x3dna-dssr"
    if not which_or(exe):
        return False, "x3dna-dssr not found in PATH or at /usr/local/bin/x3dna-dssr", []

    cmd = [exe, "fiber", "--model=Z-DNA", f"--repeat={int(repeat)}", f"-o={_path_arg(pdb_out, cwd)}"]
    try:
        cp = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
        output = (cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")
        return cp.returncode == 0, output.strip(), cmd
    except Exception as exc:
        return False, str(exc), cmd


def command_to_text(cmd: Iterable[str]) -> str:
    return " ".join(str(part) for part in cmd)


def stage_params_to_output_dir(params_file: Union[str, Path], out_dir: Union[str, Path]) -> Path:
    params_path = Path(params_file)
    if not params_path.is_absolute():
        params_path = Path.cwd() / params_path
    if not params_path.exists():
        raise FileNotFoundError(f"Phenix params not found: {params_path}")

    target_dir = Path(out_dir)
    try:
        if params_path.resolve().parent == target_dir.resolve():
            return params_path
    except Exception:
        pass

    target = target_dir / params_path.name
    shutil.copy2(params_path, target)
    return target


def run_phenix_minimization(pdb_path: Union[str, Path], params_path: Union[str, Path]) -> Tuple[bool, str, List[str]]:
    pdb = Path(pdb_path)
    params = Path(params_path) if params_path else None
    workdir = str(pdb.parent)
    pdb_name = pdb.name
    params_name = params.name if params else ""

    exe = shutil.which("phenix.geometry_minimization")
    if exe:
        cmd = [exe, pdb_name] + ([params_name] if params_name else [])
        try:
            cp = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
            output = (cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")
            return cp.returncode == 0, output.strip(), cmd
        except Exception as exc:
            return False, str(exc), cmd

    phenix_env = os.environ.get("PHENIX_ENV")
    if not phenix_env:
        candidate = "/Applications/phenix-1.21.2-5419/phenix_env.sh"
        if Path(candidate).is_file():
            phenix_env = candidate

    if phenix_env and Path(phenix_env).is_file():
        cmd_text = f"source '{phenix_env}'; phenix.geometry_minimization '{pdb_name}'"
        if params_name:
            cmd_text += f" '{params_name}'"
        cmd = ["/bin/zsh", "-lc", cmd_text]
        try:
            cp = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
            output = (cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")
            return cp.returncode == 0, output.strip(), cmd
        except Exception as exc:
            return False, str(exc), cmd

    return False, "phenix.geometry_minimization not found in PATH and PHENIX_ENV not set/invalid.", []


def expected_phenix_minimized_path(pdb_path: Union[str, Path]) -> Path:
    path = Path(pdb_path)
    suffix = path.suffix if path.suffix else ".pdb"
    return path.with_name(path.stem + "_minimized" + suffix)


def check_dssr_installation(timeout_seconds: int = 8) -> Dict[str, object]:
    """Locate x3dna-dssr and return a short installation/version report.

    The exact version flag supported by DSSR can vary by installation, so this
    helper tries a small set of safe commands and returns the first informative
    output it receives.
    """
    exe = which_or("x3dna-dssr")
    fallback = Path("/usr/local/bin/x3dna-dssr")
    if exe is None and fallback.is_file():
        exe = str(fallback)
    if exe is None:
        return {
            "installed": False,
            "executable": None,
            "command": None,
            "returncode": None,
            "output": "x3dna-dssr was not found in PATH or at /usr/local/bin/x3dna-dssr.",
        }

    attempts = [[exe, "--version"], [exe, "-v"], [exe, "--help"]]
    last_error = ""
    for cmd in attempts:
        try:
            cp = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_seconds)
            output = ((cp.stdout or "") + ("\n" + cp.stderr if cp.stderr else "")).strip()
            if output or cp.returncode == 0:
                return {
                    "installed": True,
                    "executable": exe,
                    "command": command_to_text(cmd),
                    "returncode": cp.returncode,
                    "output": output or "Command completed without text output.",
                }
            last_error = f"Command returned {cp.returncode} without output: {command_to_text(cmd)}"
        except subprocess.TimeoutExpired:
            last_error = f"Timed out while running: {command_to_text(cmd)}"
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"

    return {
        "installed": True,
        "executable": exe,
        "command": None,
        "returncode": None,
        "output": last_error or "x3dna-dssr executable was found, but no version information could be obtained.",
    }
