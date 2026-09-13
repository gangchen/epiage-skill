import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "epigenetic-clocks/scripts/compute_clocks.py"
spec = importlib.util.spec_from_file_location("clocks", SCRIPT)
clocks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clocks)


class ClockInputs(unittest.TestCase):
    def load(self, text, name="input.csv"):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / name
            path.write_text(text)
            return clocks.load_betas(path)

    def test_named_single_sample_and_masks(self):
        data = self.load("CpG,person A\ncg1,0.3\ncg2,NA\n")
        self.assertEqual(list(data.columns), ["person A"])
        self.assertTrue(pd.isna(data.loc["cg2", "person A"]))
        legacy = self.load("CpG_site,Beta_value\ncg1,0.3\n")
        self.assertEqual(list(legacy.columns), ["Sample"])

    def test_missing_resource_inventory_is_offline_and_complete(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(clocks, "DATA", directory):
            with self.assertRaises(ValueError) as error:
                clocks.check_clock_resources(["horvath", "dunedinpace"])
            for name in ["Horvath1.csv", "DunedinPACE.csv", "DunedinPACE_Gold_Means.csv", "blood_panel.npz"]:
                self.assertIn(name, str(error.exception))

    def test_tsv_and_gzip(self):
        data = self.load("CpG\tA\tB\ncg1\t0.2\t0.8\n", "input.tsv")
        self.assertEqual(data.loc["cg1", "B"], 0.8)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "betas.csv.gz"
            data.to_csv(path, compression="gzip")
            pd.testing.assert_frame_equal(data, clocks.load_betas(path))

    def test_reject_invalid_values_and_ambiguous_identifiers(self):
        invalid = ["CpG,A\ncg1,1.1\n", "CpG,A\ncg1,-0.1\n", "CpG,A\ncg1,inf\n",
                   "CpG,A\ncg1,oops\n", "CpG,A\ncg1,0.2\ncg1,0.3\n",
                   "CpG,A,A\ncg1,0.2,0.3\n", "CpG,A\n,0.2\n", "CpG,A\ncg1,NA\n",
                   "gene,A\nTP53,0.3\n", "CpG,A\ncg000001_BC11,0.3\n"]
        for text in invalid:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.load(text)

    def test_impute_masked_and_absent_independently(self):
        data = pd.DataFrame({"A": [0.2, np.nan], "B": [0.8, 0.9]}, index=["cg1", "cg2"])
        ref = SimpleNamespace(mode="global-median", sd=None,
                              median=pd.Series({"cg2": 0.4, "cg3": 0.6}))
        filled, missing, lowconf = clocks.impute_missing(data, ["cg1", "cg2", "cg3", "cg4"], ref)
        self.assertEqual(missing, {"A": ["cg2", "cg3", "cg4"], "B": ["cg3", "cg4"]})
        self.assertEqual(filled.loc["cg2", "A"], 0.4)
        self.assertEqual(filled.loc["cg2", "B"], 0.9)
        self.assertTrue(filled.loc["cg4"].isna().all())
        self.assertEqual(lowconf, {"A": 2, "B": 1})

    def test_dunedin_batch_matches_independent_samples(self):
        gold = pd.read_csv(Path(clocks.DATA) / "DunedinPACE_Gold_Means.csv", index_col=0)["mean"]
        data = pd.DataFrame({"A": gold, "B": gold})
        data.iloc[:13, 0] = np.nan
        data.iloc[20:23, 1] = np.nan
        data.iloc[24:30, 1] = 0.7
        values, n, missing = clocks.predict_dunedin(data)
        self.assertEqual(missing.to_dict(), {"A": 13, "B": 3})
        self.assertEqual(n, len(gold))
        for sample in data:
            single, _, _ = clocks.predict_dunedin(data[[sample]])
            self.assertAlmostEqual(values[sample], single[sample], places=12)

    def test_cli_per_sample_coverage_and_unavailable_model(self):
        feature = clocks.model_cpgs(clocks.CLOCKS["bocklandt"])[0]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "betas.csv"
            path.write_text(f"CpG,A,B\n{feature},0.3,NA\ncg00000000,0.4,0.4\n")
            result = subprocess.run([sys.executable, str(SCRIPT), "--input", str(path),
                                     "--clocks", "bocklandt", "education"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = json.loads(result.stdout.split("JSON:", 1)[1])
            model = {r["sample"]: r for r in rows if r["clock"] == "bocklandt"}
            self.assertEqual(model["A"]["n_missing"], 0)
            self.assertEqual(model["B"]["n_missing"], 1)
            self.assertEqual(model["B"]["n_imputed"], 1)
            self.assertTrue(all(r["status"] == "ok" for r in model.values()))
            unavailable = [r for r in rows if r["clock"] == "education"]
            self.assertTrue(all(r["status"] == "unavailable" and r["value"] is None for r in unavailable))
            self.assertTrue(all(r["n_unresolved"] > 0 for r in unavailable))


if __name__ == "__main__":
    unittest.main()
