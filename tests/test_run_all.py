import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import run_all

class OneCommandTests(unittest.TestCase):
    def test_default_has_all_steps(self):
        self.assertEqual(len(run_all.selected_steps()), 8)

    def test_skip_dependencies(self):
        names = [x[0] for x in run_all.selected_steps(skip_phylogeny=True, skip_structure=True)]
        self.assertEqual(names, ["dataset", "genus_audit", "extended", "paper_validation", "manuscript"])

    def test_missing_phylogeny_still_runs_structure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def mocked_run(args, cwd, check):
                script = Path(args[1]).name
                outputs = {
                    "run_analysis.py": "dataset_results/summary.json",
                    "audit_phylogeny.py": "dataset_results/audit/summary.json",
                    "extended_research.py": "dataset_results/extended/summary.json",
                    "structure_distances.py": "dataset_results/structures/summary.json",
                }
                dest = root / outputs[script]
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_text("{}", encoding="utf-8")
                return SimpleNamespace(returncode=0)
            with patch.object(run_all, "ROOT", root), patch.object(run_all.shutil, "which", return_value=None):
                with patch.object(run_all.subprocess, "run", side_effect=mocked_run) as call:
                    exitcode = run_all.execute()
            self.assertEqual(exitcode, 1)
            self.assertEqual(call.call_count, 4)
            import json
            manifest = json.loads((root / "dataset_results/full_run_manifest.json").read_text())
            statuses = {step["step"]: step["status"] for step in manifest["steps"]}
            self.assertEqual(statuses["phylogeny"], "missing_dependency")
            self.assertEqual(statuses["tree_support"], "blocked")
            self.assertEqual(statuses["structure"], "success")
            self.assertEqual(statuses["paper_validation"], "blocked")
            self.assertEqual(statuses["manuscript"], "blocked")

    def test_pipeline_order(self):
        names = [x[0] for x in run_all.selected_steps()]
        self.assertLess(names.index("phylogeny"), names.index("tree_support"))

if __name__ == "__main__":
    unittest.main()
