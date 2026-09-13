from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bnp_na_lib.build_common import PipelineError
from bnp_na_lib.build_triplex import (
    DEFAULT_PARAMS_FILE,
    TMP_DIR_NAME,
    build_triplex_from_duplex,
    resolve_triplex_paths,
)


class OutputRoutingTests(unittest.TestCase):
    """Where each stage writes, for the post-processing options."""

    def test_without_post_processing_conversion_writes_to_the_requested_path(self) -> None:
        final, tmp_dir, convert_target = resolve_triplex_paths("/tmp/out/my.pdb", post_process=False)

        self.assertEqual(final, Path("/tmp/out/my.pdb"))
        self.assertIsNone(tmp_dir)
        self.assertEqual(convert_target, final)

    def test_with_post_processing_conversion_is_diverted_into_the_tmp_folder(self) -> None:
        final, tmp_dir, convert_target = resolve_triplex_paths("/tmp/out/my.pdb", post_process=True)

        self.assertEqual(final, Path("/tmp/out/my.pdb"))
        self.assertEqual(tmp_dir, Path("/tmp/out") / TMP_DIR_NAME)
        self.assertEqual(convert_target, tmp_dir / "my_triplex_raw.pdb")

    def test_tmp_folder_sits_beside_the_requested_output(self) -> None:
        _final, tmp_dir, _target = resolve_triplex_paths("/data/models/t.pdb", post_process=True)
        self.assertEqual(tmp_dir.parent, Path("/data/models"))
        self.assertEqual(tmp_dir.name, "triplex_gen_tmp")

    def test_missing_suffix_still_produces_a_pdb_intermediate(self) -> None:
        _final, _tmp, convert_target = resolve_triplex_paths("/tmp/out/noext", post_process=True)
        self.assertEqual(convert_target.name, "noext_triplex_raw.pdb")

    def test_requested_path_is_never_the_conversion_target_when_post_processing(self) -> None:
        """The named output must end up holding the final model, not the raw conversion."""
        final, _tmp, convert_target = resolve_triplex_paths("/tmp/out/my.pdb", post_process=True)
        self.assertNotEqual(final, convert_target)


class OptionDefaultTests(unittest.TestCase):
    """These run far enough to hit the option checks and stop before any DSSR work."""

    @staticmethod
    def _dummy_duplex(directory: str) -> Path:
        path = Path(directory) / "duplex.pdb"
        path.write_text("END\n", encoding="utf-8")
        return path

    def _call(self, duplex: Path, out: Path, **kwargs):
        return build_triplex_from_duplex(
            duplex,
            out,
            strand_i_chain="A",
            residue_range=(2, 11),
            mode="antiparallel",
            **kwargs,
        )

    def test_bundled_params_file_exists_to_back_the_default(self) -> None:
        self.assertTrue(DEFAULT_PARAMS_FILE.is_file(), DEFAULT_PARAMS_FILE)

    def test_phenix_without_a_params_file_falls_back_to_the_bundled_one(self) -> None:
        """It must reach the conversion, not stop on a missing-params complaint."""
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            with self.assertRaises(PipelineError) as ctx:
                self._call(duplex, Path(tmp) / "out.pdb", run_phenix=True)
        self.assertNotIn("params file", str(ctx.exception))

    def test_blank_params_file_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            with self.assertRaises(PipelineError) as ctx:
                self._call(duplex, Path(tmp) / "out.pdb", run_phenix=True, params_file="  ")
        self.assertIn("no params file was specified", str(ctx.exception))

    def test_missing_params_file_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            with self.assertRaises(PipelineError) as ctx:
                self._call(
                    duplex,
                    Path(tmp) / "out.pdb",
                    run_phenix=True,
                    params_file=Path(tmp) / "nope.params",
                )
        self.assertIn("Params file not found", str(ctx.exception))

    def test_regularization_alone_does_not_need_a_params_file(self) -> None:
        """It must fail later, on the real conversion, not on the params check."""
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            with self.assertRaises(PipelineError) as ctx:
                self._call(
                    duplex,
                    Path(tmp) / "out.pdb",
                    run_phenix=False,
                    run_regularize_phosphates=True,
                )
        self.assertNotIn("params file", str(ctx.exception))

    def test_missing_duplex_is_reported_before_anything_else(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(PipelineError) as ctx:
                self._call(Path(tmp) / "nope.pdb", Path(tmp) / "out.pdb", run_phenix=True)
        self.assertIn("not found", str(ctx.exception))

    def test_blank_strand_i_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            with self.assertRaises(PipelineError) as ctx:
                build_triplex_from_duplex(
                    duplex,
                    Path(tmp) / "out.pdb",
                    strand_i_chain="  ",
                    residue_range=(2, 11),
                    mode="antiparallel",
                )
        self.assertIn("strand I", str(ctx.exception))

    def test_defaults_turn_both_stages_on(self) -> None:
        """A caller that passes no options gets minimization and regularization."""
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            out = Path(tmp) / "out.pdb"
            with self.assertRaises(PipelineError) as ctx:
                self._call(duplex, out)
            # It got past the params gate and set up the tmp folder before the
            # dummy duplex failed conversion, so post-processing was requested.
            self.assertNotIn("params file", str(ctx.exception))
            self.assertTrue((out.parent / TMP_DIR_NAME).is_dir())

    def test_both_stages_off_writes_straight_to_the_requested_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            out = Path(tmp) / "out.pdb"
            with self.assertRaises(PipelineError):
                self._call(duplex, out, run_phenix=False, run_regularize_phosphates=False)
            self.assertFalse((out.parent / TMP_DIR_NAME).exists())

    def test_regularization_follows_minimization_when_unset(self) -> None:
        """run_regularize_phosphates=None tracks run_phenix, so off means no tmp folder."""
        with tempfile.TemporaryDirectory() as tmp:
            duplex = self._dummy_duplex(tmp)
            out = Path(tmp) / "out.pdb"
            with self.assertRaises(PipelineError):
                self._call(duplex, out, run_phenix=False)
            self.assertFalse((out.parent / TMP_DIR_NAME).exists())


if __name__ == "__main__":
    unittest.main()
