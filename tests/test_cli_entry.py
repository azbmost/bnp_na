from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Optional, Sequence


APP_PATH = Path(__file__).resolve().parents[1] / "bnp_na.py"


def run_app(args: Sequence[str], *, without_tkinter: bool = False) -> subprocess.CompletedProcess:
    """Run bnp_na.py, optionally with `import tkinter` forced to fail.

    The shim stands in for a machine with no Tk installed, where the GUI cannot
    be imported but --help and --version still have to answer.
    """
    env = None
    stack: Optional[tempfile.TemporaryDirectory] = None
    try:
        if without_tkinter:
            import os

            stack = tempfile.TemporaryDirectory()
            (Path(stack.name) / "tkinter.py").write_text(
                'raise ImportError("no Tk on this machine")\n', encoding="utf-8"
            )
            env = dict(os.environ)
            env["PYTHONPATH"] = stack.name + os.pathsep + env.get("PYTHONPATH", "")
        return subprocess.run(
            [sys.executable, str(APP_PATH), *args],
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
    finally:
        if stack is not None:
            stack.cleanup()


class HelpOptionTests(unittest.TestCase):
    def test_long_help_prints_usage(self) -> None:
        result = run_app(["--help"])
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.startswith("usage: bnp_na.py"))
        self.assertIn("-h, --help", result.stdout)
        self.assertIn("-v, --version", result.stdout)

    def test_short_help_matches_long_help(self) -> None:
        self.assertEqual(run_app(["-h"]).stdout, run_app(["--help"]).stdout)

    def test_help_does_not_open_the_gui(self) -> None:
        """The whole point: --help must answer instead of launching a window."""
        result = run_app(["--help"], without_tkinter=True)
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.startswith("usage: bnp_na.py"))

    def test_help_names_the_helper_cli_tools(self) -> None:
        self.assertIn("bnp_na_lib/", run_app(["--help"]).stdout)

    def test_help_reports_the_current_version(self) -> None:
        version = run_app(["--version"]).stdout.split()[-1]
        self.assertIn(version, run_app(["--help"]).stdout)


class VersionOptionTests(unittest.TestCase):
    def test_long_version(self) -> None:
        result = run_app(["--version"])
        self.assertEqual(result.returncode, 0)
        self.assertRegex(result.stdout.strip(), r"^bnp_na V\d+\.\d+$")

    def test_short_version_matches_long_version(self) -> None:
        self.assertEqual(run_app(["-v"]).stdout, run_app(["--version"]).stdout)

    def test_version_works_without_tkinter(self) -> None:
        result = run_app(["--version"], without_tkinter=True)
        self.assertEqual(result.returncode, 0)
        self.assertRegex(result.stdout.strip(), r"^bnp_na V\d+\.\d+$")

    def test_help_wins_when_both_are_given(self) -> None:
        result = run_app(["--version", "--help"])
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.startswith("usage: bnp_na.py"))


class UnrecognizedArgumentTests(unittest.TestCase):
    def test_unknown_flag_is_a_usage_error(self) -> None:
        """It must not fall through and open the GUI, which is what it used to do."""
        result = run_app(["--nope"], without_tkinter=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized argument: --nope", result.stderr)
        self.assertIn("--help", result.stderr)
        self.assertEqual(result.stdout, "")

    def test_stray_positional_is_a_usage_error(self) -> None:
        result = run_app(["helix.pdb"], without_tkinter=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized argument: helix.pdb", result.stderr)

    def test_first_unrecognized_argument_is_named(self) -> None:
        result = run_app(["--one", "--two"], without_tkinter=True)
        self.assertIn("unrecognized argument: --one", result.stderr)


class ArgumentDispatchTests(unittest.TestCase):
    """Direct checks on the helper, without paying for a subprocess."""

    @staticmethod
    def _load_dispatcher():
        # bnp_na imports tkinter at module scope, so read the helper out of the
        # pre-import block rather than importing the module itself.
        source = APP_PATH.read_text(encoding="utf-8")
        start = source.index("HELP_TEXT = ")
        end = source.index('if __name__ == "__main__":', start)
        namespace: dict = {"sys": sys, "List": list, "APP_NAME": "bnp_na", "__version__": "V0.0"}
        exec(compile(source[start:end], str(APP_PATH), "exec"), namespace)
        return namespace["_answer_cli_request"]

    def test_no_arguments_returns_and_lets_the_gui_start(self) -> None:
        self.assertIsNone(self._load_dispatcher()([]))

    def test_each_recognized_flag_exits_zero(self) -> None:
        dispatch = self._load_dispatcher()
        for flag in ("-h", "--help", "-v", "--version"):
            with self.subTest(flag=flag), self.assertRaises(SystemExit) as ctx:
                dispatch([flag])
            self.assertEqual(ctx.exception.code, 0)

    def test_unknown_argument_exits_two(self) -> None:
        with self.assertRaises(SystemExit) as ctx:
            self._load_dispatcher()(["--bogus"])
        self.assertEqual(ctx.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
