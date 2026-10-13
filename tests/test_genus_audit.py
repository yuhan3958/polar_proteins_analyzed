"""Deterministic regression tests for the genus-audit calculations."""
import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from scripts import audit_phylogeny as audit


class GenusAuditTests(unittest.TestCase):
    def test_mean(self):
        self.assertEqual(audit.mean([2, 4, 6]), 4)
        self.assertIsNone(audit.mean([]))

    def test_fasta_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "one.fasta"
            fasta.write_text(">A\nACDG\n", encoding="utf-8")
            self.assertEqual(audit.read_sequence(fasta), "ACDG")
            fasta.write_text(">A\nACD?\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                audit.read_sequence(fasta)
            fasta.write_text(">A\nACD\n>B\nACD\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                audit.read_sequence(fasta)

    def test_end_to_end_small_dataset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "data").mkdir()
            rows = [
                ("c1", "cold", "Alpha sp1", "GGGG"),
                ("c2", "cold", "Alpha sp2", "GGGG"),
                ("c3", "cold", "Beta sp1", "PPPP"),
                ("m1", "meso", "Gamma sp1", "RRRR"),
            ]
            fields = ["accession", "group", "organism", "fasta_path",
                      "sequence_length", "growth_optimum_c",
                      "temperature_source_url", "enzyme_family", "thermal_class"]
            with (root / "data" / "metadata.csv").open("w", encoding="utf-8", newline="") as output:
                writer = csv.DictWriter(output, fieldnames=fields)
                writer.writeheader()
                for accession, group, species, sequence in rows:
                    filename = f"{accession}.fasta"
                    (root / filename).write_text(f">{accession}\n{sequence}\n", encoding="utf-8")
                    writer.writerow(dict(accession=accession, group=group,
                                         organism=species, fasta_path=filename,
                                         sequence_length=len(sequence),
                                         growth_optimum_c="", temperature_source_url="",
                                         enzyme_family="GH2", thermal_class=""))
            with patch.object(audit, "ROOT", root), patch.object(audit, "INPUT", root / "data" / "metadata.csv"), patch.object(audit, "OUTPUT", root / "results"):
                audit.audit()
                with (root / "results" / "group_sensitivity.csv").open(encoding="utf-8-sig") as output:
                    g = next(r for r in csv.DictReader(output) if r["amino_acid"] == "G")
                # Naive cold G%: (100 + 100 + 0)/3; equal-genus: (100 + 0)/2.
                self.assertAlmostEqual(float(g["cold_protein_mean"]), 200 / 3)
                self.assertAlmostEqual(float(g["cold_equal_genus_mean"]), 50)
                with (root / "results" / "leave_one_genus_out_contrasts.csv").open(encoding="utf-8-sig") as output:
                    contrast = next(r for r in csv.DictReader(output)
                                    if r["amino_acid"] == "G"
                                    and r["excluded_group"] == "cold"
                                    and r["excluded_genus"] == "Alpha")
                self.assertAlmostEqual(float(contrast["baseline_pp"]), 50)
                self.assertAlmostEqual(float(contrast["cold_minus_meso_pp"]), 0)
                self.assertAlmostEqual(float(contrast["change_from_baseline_pp"]), -50)


if __name__ == "__main__":
    unittest.main()
