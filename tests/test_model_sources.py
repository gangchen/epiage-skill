"""Regression oracles from original coefficients and actual output semantics.

DNAmTL fixture: Lu et al. 2019 Supplementary Data 1, worksheet
DNAmLTL_addgene_addmQTL, rows 7-147. Provenance: references/model-audit.md.
No personal methylation data or network access is needed.
"""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "epigenetic-clocks/scripts/compute_clocks.py"
spec = importlib.util.spec_from_file_location("source_audited_clocks", SCRIPT)
clocks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clocks)


class ModelSources(unittest.TestCase):
    def test_mccartney_weights_and_raw_scores_match_original_supplement(self):
        source = pd.read_csv(ROOT / "tests/fixtures/McCartney2018_original.tsv", sep="\t")
        anchors = json.loads((ROOT / "tests/fixtures/McCartney2018_expected.json").read_text())
        for key, expected in anchors.items():
            with self.subTest(model=key):
                original = source[source.model.eq(key)].set_index("CpG")["coefficient"].sort_index()
                bundled = pd.read_csv(Path(clocks.DATA) / clocks.CLOCKS[key]["file"], index_col=0).iloc[:, 0].sort_index()
                self.assertEqual(list(original.index), list(bundled.index))
                self.assertEqual(len(original), expected["n_cpgs"])
                # Source tables carry more precision than the historical CSVs.
                np.testing.assert_allclose(bundled, original, rtol=0, atol=5e-8)
                data = pd.DataFrame({"half": 0.5}, index=original.index)
                actual = clocks.predict_linear(data, clocks.CLOCKS[key])["half"]
                self.assertAlmostEqual(actual, expected["expected_weighted_sum"], places=6)

    def test_unreconstructed_models_return_no_scores_while_other_models_continue(self):
        keys = ["cvd", "depression", "dnamtl"]
        features = sorted({c for key in keys for c in clocks.model_cpgs(clocks.CLOCKS[key])})
        data = pd.DataFrame({"synthetic": 0.5}, index=features)
        for key in ("cvd", "depression"):
            with self.subTest(model=key), self.assertRaises(ValueError):
                clocks.predict_linear(data, clocks.CLOCKS[key])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "betas.csv"
            data.to_csv(path, index_label="CpG")
            result = subprocess.run([sys.executable, str(SCRIPT), "--input", str(path),
                                     "--clocks", *keys], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = {r["clock"]: r for r in json.loads(result.stdout.split("JSON:", 1)[1])}
        for key in ("cvd", "depression"):
            self.assertEqual(rows[key]["status"], "unavailable")
            self.assertIsNone(rows[key]["value"])
            self.assertEqual(rows[key]["n_imputed"], 0)
            self.assertTrue(rows[key]["status_reason"])
            self.assertIn(key, result.stderr)
        self.assertEqual(rows["dnamtl"]["status"], "ok")
        self.assertEqual(rows["dnamtl"]["value"], 8.42)

    def test_dnamtl_all_weights_and_intercept_match_original_supplement(self):
        source = pd.read_csv(ROOT / "tests/fixtures/DNAmTL_Lu2019_original.csv", index_col=0)
        bundled = pd.read_csv(Path(clocks.DATA) / "DNAmTL.csv", index_col=0)
        self.assertEqual(len(source), 141)
        self.assertTrue(source.index.is_unique and bundled.index.is_unique)
        pd.testing.assert_frame_equal(bundled.sort_index(), source.sort_index(), check_exact=True)

    def test_dnamtl_source_anchored_predictions_and_feature_alignment(self):
        source = pd.read_csv(ROOT / "tests/fixtures/DNAmTL_Lu2019_original.csv", index_col=0).iloc[:, 0]
        weights = source.drop("intercept")
        data = pd.DataFrame({"zero": 0.0, "half": 0.5, "one": 1.0,
                             "gradient": np.linspace(0.05, 0.95, len(weights))}, index=weights.index)
        expected = pd.Series({"zero": 7.924780053, "half": 8.417606178,
                              "one": 8.910432303,
                              "gradient": 7.924780053 + weights.dot(data["gradient"])})
        # Reversing input rows must not change which CpG each coefficient uses.
        actual = clocks.predict_linear(data.iloc[::-1], clocks.CLOCKS["dnamtl"])
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    def test_cli_age_and_non_age_units_and_acceleration(self):
        keys = ["garagnani", "bocklandt", "dnamtl", "stocz", "zhang", "epitoc1", "dunedinpoam"]
        features = sorted({c for key in keys for c in clocks.model_cpgs(clocks.CLOCKS[key])})
        data = pd.DataFrame({"synthetic": 0.5}, index=features)
        data.loc["cg16867657", "synthetic"] = 0.3
        data.loc["cg09809672", "synthetic"] = 0.7
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "betas.csv"
            data.to_csv(path, index_label="CpG")
            result = subprocess.run([sys.executable, str(SCRIPT), "--input", str(path),
                                     "--age", "45", "--clocks", *keys],
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = {r["clock"]: r for r in json.loads(result.stdout.split("JSON:", 1)[1])}
        for key, expected in (("garagnani", 0.3), ("bocklandt", 0.7)):
            self.assertEqual(rows[key]["unit"], "beta")
            self.assertEqual(rows[key]["value"], expected)
        for key in set(keys) - {"stocz"}:
            self.assertEqual(rows[key]["accel"], "", key)
            self.assertEqual(rows[key]["n_missing"], 0, key)
        self.assertEqual(rows["dnamtl"]["unit"], "kb")
        self.assertEqual(rows["zhang"]["unit"], "score")
        self.assertEqual(rows["epitoc1"]["unit"], "score")
        self.assertEqual(rows["dunedinpoam"]["unit"], "years/year")
        self.assertEqual(rows["stocz"]["unit"], "years")
        self.assertAlmostEqual(rows["stocz"]["accel"], rows["stocz"]["value"] - 45, places=2)


if __name__ == "__main__":
    unittest.main()
