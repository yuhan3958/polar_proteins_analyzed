"""Extended GH2 analyses: thermal evidence audit, catalytic-site mapping,
and exploratory pairwise-distance clustering.

python scripts/extended_research.py

This does NOT reconstruct a statistically validated species/protein phylogeny.
For that use a quality-controlled multiple sequence alignment, an evolutionary
substitution model, bootstraps and relevant outgroups.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from Bio import Align, Phylo, SeqIO
from Bio.Align import substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "dataset_results" / "extended"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
# E. coli LacZ residues established by literature. UniProt P00722.
SITES = {388: "Arg accessory", 461: "Glu acid/base",
         502: "Met accessory", 503: "Tyr accessory",
         537: "Glu nucleophile"}
PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC92765/"
PDB_EVIDENCE = "https://pubmed.ncbi.nlm.nih.gov/11732897/"


def read_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{path}: exactly one FASTA record needed")
    seq = str(records[0].seq).upper()
    if not seq or set(seq) - STANDARD:
        raise ValueError(f"{path}: noncanonical protein sequence")
    return seq


def save_csv(path, rows, headers):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def aligner():
    x = Align.PairwiseAligner(mode="global")
    x.substitution_matrix = substitution_matrices.load("BLOSUM62")
    x.open_gap_score = -10
    x.extend_gap_score = -0.5
    return x


def comparable_identity(a, b, engine):
    aln = engine.align(a, b)[0]
    paired = [(int(i), int(j)) for i, j in zip(*aln.indices) if i >= 0 and j >= 0]
    return (sum(a[i] == b[j] for i, j in paired) / len(paired)) if paired else 0.


def map_reference_sites(reference, target, engine, flank=15):
    aln = engine.align(reference, target)[0]
    ref_to_target = {int(i) + 1: int(j) + 1
                     for i, j in zip(*aln.indices) if i >= 0 and j >= 0}
    rows = []
    for ref_pos, role in SITES.items():
        if ref_pos > len(reference):
            raise ValueError("Reference is shorter than annotated catalytic site")
        target_pos = ref_to_target.get(ref_pos)
        residue = target[target_pos - 1] if target_pos else "-"
        start = max(1, ref_pos - flank)
        end = min(len(reference), ref_pos + flank)
        mapped = [ref_to_target[p] for p in range(start, end + 1)
                  if p in ref_to_target]
        target_context = [target[p - 1] for p in mapped]
        rows.append({
            "reference_position": ref_pos,
            "role": role,
            "reference_residue": reference[ref_pos - 1],
            "target_position": target_pos or "",
            "target_residue": residue,
            "conserved": residue == reference[ref_pos - 1],
            "reference_window_G_percent": 100 * reference[start - 1:end].count("G") / (end - start + 1),
            "target_mapped_window_G_percent": (
                100 * target_context.count("G") / len(target_context)
                if target_context else ""),
            "mapped_window_pairs": len(target_context),
        })
    return rows


def audit_thermal_metadata(meta):
    rows = []
    for x in meta:
        optimum = x.get("growth_optimum_c", "").strip()
        culture = x.get("culture_temperature_c", "").strip()
        enzyme = x.get("enzyme_activity_optimum_c", "").strip()
        source = x.get("temperature_source_url", "").strip()
        group = x["group"]
        classification = x.get("thermal_class", "")
        # 'Known' means recorded in the table, NOT independently source-verified.
        rows.append({
            "accession": x["accession"], "group": group,
            "organism": x["organism"], "thermal_class": classification,
            "growth_optimum_c": optimum, "culture_temperature_c": culture,
            "enzyme_activity_optimum_c": enzyme,
            "source_url": source,
            "growth_optimum_recorded": bool(optimum),
            "source_link_recorded": bool(source),
            "cold_tolerant_label": (group == "cold" and
                any(t in classification.lower() for t in ("psychrotolerant", "psychrotrophic"))),
            "needs_independent_source_review": not bool(optimum),
        })
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with (ROOT / "data" / "metadata.csv").open(encoding="utf-8-sig", newline="") as file:
        meta = list(csv.DictReader(file))
    sequences = {}
    for x in meta:
        filename = Path(x["fasta_path"])
        path = filename if filename.is_absolute() else ROOT / filename
        sequences[x["accession"]] = read_fasta(path)
    reference = sequences["P00722"]
    engine = aligner()
    sites = []
    for x in meta:
        accession = x["accession"]
        for row in map_reference_sites(reference, sequences[accession], engine):
            sites.append({"accession": accession, "group": x["group"], **row})
    save_csv(OUT / "active_site_mapping.csv", sites,
             ["accession", "group", "reference_position", "role",
              "reference_residue", "target_position", "target_residue",
              "conserved", "reference_window_G_percent",
              "target_mapped_window_G_percent", "mapped_window_pairs"])
    audit = audit_thermal_metadata(meta)
    save_csv(OUT / "thermal_metadata_audit.csv", audit, list(audit[0]))
    # Exploratory neighbor joining of 1 - pairwise non-gap identity;
    # it is not a site-consistent MSA-based evolutionary phylogeny.
    names = list(sequences)
    lower = []
    for i, item in enumerate(names):
        lower.append([round(1 - comparable_identity(
            sequences[item], sequences[names[j]], engine), 8)
            if i != j else 0.0 for j in range(i + 1)])
    matrix = DistanceMatrix(names=names, matrix=lower)
    tree = DistanceTreeConstructor().nj(matrix)
    Phylo.write(tree, OUT / "identity_distance_NJ_EXPLORATORY.nwk", "newick")
    group_stats = {}
    for group in ["cold", "meso"]:
        subset = [r for r in audit if r["group"] == group]
        group_stats[group] = {
            "n": len(subset),
            "growth_optimum_recorded": sum(r["growth_optimum_recorded"] for r in subset),
            "thermal_source_url_recorded": sum(r["source_link_recorded"] for r in subset),
            "cold_tolerant_label": sum(r["cold_tolerant_label"] for r in subset),
        }
    stats = {
        "group_metadata": group_stats,
        "active_site_rows": len(sites),
        "active_sites_reference": SITES,
        "method": "BLOSUM62 global alignment, open_gap=-10, extend_gap=-0.5",
        "site_source": PAPER,
        "structure_source": PDB_EVIDENCE,
        "warning": "Distances are 1 - pairwise aligned nongap identity. NJ topology is exploratory, NOT a reliable phylogeny, without an MSA, substitution model and bootstrap. Mapping is sequence-based, not measured 3D proximity. Metadata link recorded does not mean source verified.",
    }
    (OUT / "summary.json").write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n",
                                      encoding="utf-8")
    report = [
        "# Extended GH2 investigation",
        "",
        "## Thermal metadata completeness",
        "",
        "| Group | Proteins | Growth optima recorded | Temperature source URLs |",
        "|---|---:|---:|---:|",
    ]
    for group, s in group_stats.items():
        report.append(f"| {group} | {s['n']} | {s['growth_optimum_recorded']} | {s['thermal_source_url_recorded']} |")
    report += [
        "",
        "A missing optimum is not repaired by substituting an enzyme activity optimum or culture temperature.",
        "Entries with source URLs still need independent reading and strain matching.",
        "",
        "## Catalytic site alignment",
        "",
        "Residues are mapped to the E. coli P00722 reference sequence (Arg388, Glu461, Met502, Tyr503, Glu537).",
        "The local G percentages are for aligned sequence positions, not a structural radius.",
        "",
        "## Exploratory distance tree",
        "",
        "An NJ tree is written from pairwise (1 - identity) distances as a **diagnostic only**.",
        "It must not be presented as a validated evolutionary phylogeny; rigorous inference needs a vetted multiple sequence alignment and evolutionary model.",
        "",
        "## Caveats",
        "",
        "- The protein dataset remains taxonomically imbalanced.",
        "- No 3D structures or catalytic activity data were newly measured here.",
        "- The identity-derived distances may violate tree metric assumptions.",
        "- Source URLs and annotated residues require critical verification.",
    ]
    (OUT / "research_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(stats, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
