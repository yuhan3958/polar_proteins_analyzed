import csv
import json
import tempfile
import unittest
from pathlib import Path
from scripts.validate_manuscript import check

class ManuscriptValidationTests(unittest.TestCase):
    def test_consistent_synthetic_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            data=root/"dataset_results"
            for dirname in ("audit","phylogeny","structures"):
                (data/dirname).mkdir(parents=True)
            (root/"data").mkdir()
            def write(rel,obj):
                path=data/rel
                path.write_text(json.dumps(obj),encoding="utf-8")
            with (root/"data/metadata.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["accession","group","growth_optimum_c"])
                w.writeheader()
                w.writerows([{"accession":"a","group":"cold","growth_optimum_c":""},
                             {"accession":"b","group":"meso","growth_optimum_c":""}])
            write("summary.json",{"group_counts":{"cold":1,"meso":1},"total_pairs":1})
            with (data/"audit/group_sensitivity.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["amino_acid","cold_protein_mean","meso_protein_mean",
                    "protein_weighted_difference_pp","cold_equal_genus_mean",
                    "meso_equal_genus_mean","equal_genus_difference_pp"])
                w.writeheader()
                for aa in "GPR":
                    w.writerow(dict(amino_acid=aa,cold_protein_mean=3,meso_protein_mean=2,
                                    protein_weighted_difference_pp=1,cold_equal_genus_mean=4,
                                    meso_equal_genus_mean=2,equal_genus_difference_pp=2))
            with (data/"audit/genus_counts.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["group","genus","n"])
                w.writeheader()
                w.writerows([{"group":"cold","genus":"Alpha","n":1},
                             {"group":"meso","genus":"Beta","n":1}])
            write("phylogeny/support_audit.json",
                  {"leaf_count":2,"nodes_with_support":0,"internal_nodes":0,
                   "local_support_below_0_70":0,"local_support_below_0_90":0})
            write("structures/summary.json",{"n_distances":1})
            with (data/"structures/distances.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["ca_distance_angstrom","nearest_heavy_atom_angstrom"])
                w.writeheader()
                w.writerow({"ca_distance_angstrom":4,"nearest_heavy_atom_angstrom":2})
            self.assertTrue(check(root)["passed"])
            write("summary.json",{"group_counts":{"cold":1,"meso":1},"total_pairs":2})
            with self.assertRaises(ValueError):
                check(root)

if __name__ == "__main__":
    unittest.main()
