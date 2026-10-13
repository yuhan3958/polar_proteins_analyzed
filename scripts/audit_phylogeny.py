"""Audit thermal-group labels and genus imbalance in the GH2 dataset.

Run: python scripts/audit_phylogeny.py
Outputs CSV/JSON in dataset_results/audit/.
Uses only the Python standard library; this is genus-stratified sensitivity
analysis, NOT a phylogenetic comparative model or proof of independence.
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "metadata.csv"
OUTPUT = ROOT / "dataset_results" / "audit"
AA = "ACDEFGHIKLMNPQRSTVWY"
FOCUS = "GPR"


def read_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = [x for x in lines if x.startswith(">")]
    if len(headers) != 1:
        raise ValueError(f"{path}: expected exactly one FASTA record")
    sequence = "".join(x.strip() for x in lines if x and not x.startswith(">")).upper()
    if not sequence or set(sequence) - set(AA):
        raise ValueError(f"{path}: invalid amino-acid sequence")
    return sequence


def mean(values):
    return sum(values) / len(values) if values else None


def save_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def audit():
    with INPUT.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("No metadata rows")
    OUTPUT.mkdir(parents=True, exist_ok=True)

    observations, flags = [], []
    seen_accessions = set()
    for row in rows:
        accession = row["accession"]
        group = row["group"]
        if not accession or accession in seen_accessions:
            raise ValueError(f"Duplicate or empty accession: {accession!r}")
        seen_accessions.add(accession)
        if group not in {"cold", "meso"}:
            raise ValueError(f"{accession}: unexpected group {group!r}")
        organism = row["organism"].strip()
        genus = organism.split()[0] if organism else "UNKNOWN"
        fasta = Path(row["fasta_path"])
        if not fasta.is_absolute():
            fasta = ROOT / fasta
        sequence = read_sequence(fasta)
        if len(sequence) != int(row["sequence_length"]):
            raise ValueError(f"{accession}: inconsistent sequence length")
        freq = Counter(sequence)
        obs = {"accession": accession, "group": group, "genus": genus,
               "organism": organism, "length": len(sequence)}
        for aa in FOCUS:
            obs[aa] = 100 * freq[aa] / len(sequence)
        observations.append(obs)

        if not row.get("growth_optimum_c", "").strip():
            flags.append({"accession": accession, "flag": "growth_optimum_unknown",
                          "detail": "Growth optimum unavailable; culture temperature is not equivalent"})
        if not row.get("temperature_source_url", "").strip():
            flags.append({"accession": accession, "flag": "thermal_source_missing",
                          "detail": "Group assignment has no temperature source URL"})
        if row.get("enzyme_family", "").strip() != "GH2":
            flags.append({"accession": accession, "flag": "family_not_GH2",
                          "detail": row.get("enzyme_family", "")})
        if "psychrotolerant" in row.get("thermal_class", "").lower() or "psychrotrophic" in row.get("thermal_class", "").lower():
            flags.append({"accession": accession, "flag": "cold_tolerant_not_necessarily_psychrophilic",
                          "detail": row.get("thermal_class", "")})
        if genus == "UNKNOWN":
            flags.append({"accession": accession, "flag": "genus_unknown", "detail": organism})

    groups = defaultdict(list)
    for obs in observations:
        groups[obs["group"]].append(obs)

    genus_rows = []
    for group in ("cold", "meso"):
        by_genus = defaultdict(list)
        for obs in groups[group]:
            by_genus[obs["genus"]].append(obs)
        for genus, subset in sorted(by_genus.items()):
            genus_rows.append({"group": group, "genus": genus, "n": len(subset),
                               **{aa: mean([x[aa] for x in subset]) for aa in FOCUS}})
        if not by_genus:
            raise ValueError(f"No observations in group {group}")

    result_rows = []
    by_group_genus = defaultdict(list)
    for row in genus_rows:
        by_group_genus[row["group"]].append(row)
    for aa in FOCUS:
        cold_naive = mean([x[aa] for x in groups["cold"]])
        meso_naive = mean([x[aa] for x in groups["meso"]])
        cold_balanced = mean([x[aa] for x in by_group_genus["cold"]])
        meso_balanced = mean([x[aa] for x in by_group_genus["meso"]])
        result_rows.append({
            "amino_acid": aa,
            "cold_protein_mean": cold_naive,
            "meso_protein_mean": meso_naive,
            "protein_weighted_difference_pp": cold_naive - meso_naive,
            "cold_equal_genus_mean": cold_balanced,
            "meso_equal_genus_mean": meso_balanced,
            "equal_genus_difference_pp": cold_balanced - meso_balanced,
        })

    loo = []
    for group in ("cold", "meso"):
        for genus_row in by_group_genus[group]:
            remaining = [x for x in by_group_genus[group] if x["genus"] != genus_row["genus"]]
            if not remaining:
                continue
            for aa in FOCUS:
                loo.append({"group": group, "excluded_genus": genus_row["genus"],
                            "amino_acid": aa, "remaining_genus_count": len(remaining),
                            "equal_genus_mean": mean([x[aa] for x in remaining])})

    # Repeat full cold-minus-meso contrasts after excluding each genus in turn.
    # This is a stability diagnostic, not an independent-sample confidence interval.
    loo_contrasts = []
    for group in ("cold", "meso"):
        for genus_row in by_group_genus[group]:
            retained = {g: [x for x in by_group_genus[g]
                             if not (g == group and x["genus"] == genus_row["genus"])]
                        for g in ("cold", "meso")}
            if not retained["cold"] or not retained["meso"]:
                continue
            for aa in FOCUS:
                difference = (mean([x[aa] for x in retained["cold"]])
                              - mean([x[aa] for x in retained["meso"]]))
                baseline = (mean([x[aa] for x in by_group_genus["cold"]])
                            - mean([x[aa] for x in by_group_genus["meso"]]))
                loo_contrasts.append({
                    "excluded_group": group,
                    "excluded_genus": genus_row["genus"],
                    "amino_acid": aa,
                    "remaining_cold_genera": len(retained["cold"]),
                    "remaining_meso_genera": len(retained["meso"]),
                    "cold_minus_meso_pp": difference,
                    "baseline_pp": baseline,
                    "change_from_baseline_pp": difference - baseline,
                    "sign_reversal": (difference * baseline < 0),
                })

    save_csv(OUTPUT / "leave_one_genus_out_contrasts.csv", loo_contrasts,
             ["excluded_group", "excluded_genus", "amino_acid",
              "remaining_cold_genera", "remaining_meso_genera",
              "cold_minus_meso_pp", "baseline_pp",
              "change_from_baseline_pp", "sign_reversal"])
    save_csv(OUTPUT / "genus_counts.csv", genus_rows, ["group", "genus", "n", *FOCUS])
    save_csv(OUTPUT / "group_sensitivity.csv", result_rows,
             ["amino_acid", "cold_protein_mean", "meso_protein_mean",
              "protein_weighted_difference_pp", "cold_equal_genus_mean",
              "meso_equal_genus_mean", "equal_genus_difference_pp"])
    save_csv(OUTPUT / "leave_one_genus_out.csv", loo,
             ["group", "excluded_genus", "amino_acid", "remaining_genus_count", "equal_genus_mean"])
    save_csv(OUTPUT / "classification_flags.csv", flags, ["accession", "flag", "detail"])
    summary = {"n_proteins": len(observations),
               "group_sizes": {g: len(groups[g]) for g in ("cold", "meso")},
               "n_genera": {g: len(by_group_genus[g]) for g in ("cold", "meso")},
               "n_classification_flags": len(flags),
               "warning": ("Equal-genus weights are a sensitivity check, not phylogenetic "
                           "correction. Growth optimum missingness and group-definition "
                           "differences prohibit causal claims.")}
    (OUTPUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = [
        "# GH2 taxon-sampling sensitivity report",
        "",
        f"Sequences: {len(observations)}; cold genera: {len(by_group_genus['cold'])}; "
        f"mesophilic genera: {len(by_group_genus['meso'])}.",
        "",
        "## Cold minus mesophilic amino-acid percentages (percentage points)",
        "",
        "| Residue | Protein-weighted | Equal-genus-weighted | LOO contrast range | Sign reversals |",
        "|---|---:|---:|---:|---:|",
    ]
    for result in result_rows:
        aa = result["amino_acid"]
        group_exclusions = [x for x in loo_contrasts if x["amino_acid"] == aa]
        vals = [x["cold_minus_meso_pp"] for x in group_exclusions]
        flip_count = sum(x["sign_reversal"] for x in group_exclusions)
        interval = f"{min(vals):+.3f} to {max(vals):+.3f}" if vals else "N/A"
        report.append(
            f"| {aa} | {result['protein_weighted_difference_pp']:+.3f} "
            f"| {result['equal_genus_difference_pp']:+.3f} "
            f"| {interval} | {flip_count} / {len(vals)} |"
        )
    report.extend([
        "",
        "## Methodological qualifications",
        "",
        "- Equal-genus weighting and leave-one-genus-out are sampling-sensitivity checks, "
        "**not** phylogenetic correction.",
        "- These are descriptive contrasts without p-values or causal attribution.",
        "- One genus is omitted from a single group at a time; remaining genera receive "
        "equal weight within their group.",
        "- Growth optimum and enzyme activity optimum differ. Missing growth-optimum "
        "metadata precludes temperature-response regression.",
        "- LOO ranges are **not** confidence intervals. The number of alignments "
        "is not the number of independent observations.",
        "",
    ])
    (OUTPUT / "research_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    for row in result_rows:
        print(f"{row['amino_acid']}: protein Δ={row['protein_weighted_difference_pp']:+.3f} pp; "
              f"equal-genus Δ={row['equal_genus_difference_pp']:+.3f} pp")


if __name__ == "__main__":
    audit()
