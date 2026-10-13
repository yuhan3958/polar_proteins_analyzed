"""Pre-publication consistency checks for GH2 exploratory manuscript.

Produces a machine-readable audit and rejects contradictory/stale results.
This checks computational consistency, not biological source authenticity.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "dataset_results"
OUTPUT = DATA / "manuscript_validation.json"

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))

def check(root=ROOT):
    data = root / "dataset_results"
    metadata = read_csv(root / "data" / "metadata.csv")
    summary = json.loads((data / "summary.json").read_text(encoding="utf-8"))
    group = read_csv(data / "audit" / "group_sensitivity.csv")
    genus = read_csv(data / "audit" / "genus_counts.csv")
    support = json.loads((data / "phylogeny" / "support_audit.json").read_text(encoding="utf-8"))
    structure = json.loads((data / "structures" / "summary.json").read_text(encoding="utf-8"))
    distances = read_csv(data / "structures" / "distances.csv")
    errors, warnings = [], []
    counts = {g:sum(x["group"] == g for x in metadata) for g in ("cold", "meso")}
    if counts != summary["group_counts"]:
        errors.append("Metadata group counts contradict dataset summary")
    n = len(metadata)
    if n * (n - 1) // 2 != summary["total_pairs"]:
        errors.append("All-pair count inconsistent with n choose 2")
    if len({x["accession"] for x in metadata}) != n:
        errors.append("Duplicate accession in metadata")
    if support["leaf_count"] != n:
        errors.append("Tree leaf count contradicts metadata")
    if support["nodes_with_support"] > support["internal_nodes"]:
        errors.append("Impossible tree support count")
    if support["local_support_below_0_70"] > support["local_support_below_0_90"]:
        errors.append("Support threshold counts non-monotonic")
    if len(group) != 3 or {x["amino_acid"] for x in group} != {"G", "P", "R"}:
        errors.append("Expected exactly G/P/R comparison rows")
    for row in group:
        cold, meso = float(row["cold_protein_mean"]), float(row["meso_protein_mean"])
        diff = float(row["protein_weighted_difference_pp"])
        balanced = float(row["equal_genus_difference_pp"])
        genus_delta = float(row["cold_equal_genus_mean"]) - float(row["meso_equal_genus_mean"])
        if not math.isclose(cold - meso, diff, abs_tol=1e-8):
            errors.append("Protein-weighted contrast mismatch: " + row["amino_acid"])
        if not math.isclose(genus_delta, balanced, abs_tol=1e-8):
            errors.append("Equal-genus contrast mismatch: " + row["amino_acid"])
    genus_sizes = {}
    for g in ("cold", "meso"):
        subset = [x for x in genus if x["group"] == g]
        genus_sizes[g] = len(subset)
        if sum(int(x["n"]) for x in subset) != counts[g]:
            errors.append("Genus sample counts inconsistent: " + g)
    if genus_sizes["cold"] < 3:
        warnings.append("Cold group has very few genera; avoid independence claims")
    if any(not x.get("growth_optimum_c","").strip() for x in metadata):
        warnings.append("Growth temperature incompletely recorded; no continuous-temperature regression")
    if len(distances) != structure["n_distances"]:
        errors.append("PDB distance rows differ from structure summary")
    if any(float(row["ca_distance_angstrom"]) <= 0 or
           float(row["nearest_heavy_atom_angstrom"]) <= 0 or
           float(row["nearest_heavy_atom_angstrom"]) > float(row["ca_distance_angstrom"])
           for row in distances):
        errors.append("Impossible or nonpositive atom distances")
    report = {
        "passed": not errors,
        "errors": errors,
        "warnings": warnings,
        "metadata_n": n,
        "counts": counts,
        "genus_counts": genus_sizes,
        "pair_count": summary["total_pairs"],
        "distance_count": len(distances),
        "scope": "Checks internal consistency; original papers and experimental conclusions require independent scientific review."
    }
    output = data / "manuscript_validation.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if errors:
        raise ValueError("Manuscript consistency validation failed: " + "; ".join(errors))
    return report

if __name__ == "__main__":
    check()
