"""Generate a reproducible protein multiple sequence alignment and exploratory
approximate maximum-likelihood tree.

External command-line tools: MAFFT and FastTree.
Run: python scripts/build_phylogeny.py
These sequences are proteins, so -wag is a protein substitution model.
Do not interpret this gene tree as a species tree or a correction for temperature.
"""
from __future__ import annotations
import csv
import json
import shutil
import subprocess
from pathlib import Path
from Bio import AlignIO, Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "dataset_results" / "phylogeny"

def commands(mafft="mafft", fasttree="FastTree"):
    if not shutil.which(mafft) or not shutil.which(fasttree):
        raise RuntimeError("MAFFT and FastTree are required. Install system binaries first.")
    return mafft, fasttree

def write_input():
    with (ROOT / "data/metadata.csv").open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or len(set(x["accession"] for x in rows)) != len(rows):
        raise ValueError("Missing or duplicated accessions")
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "input.fasta"
    with path.open("w", encoding="ascii") as f:
        for row in rows:
            raw = Path(row["fasta_path"])
            raw = raw if raw.is_absolute() else ROOT / raw
            records = list(SeqIO.parse(raw, "fasta"))
            if len(records) != 1:
                raise ValueError(f"{raw}: expected one FASTA")
            seq = str(records[0].seq).upper()
            if len(seq) != int(row["sequence_length"]):
                raise ValueError(f"{row['accession']}: length mismatch")
            if not seq or set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
                raise ValueError(f"{row['accession']}: invalid sequence")
            f.write(f">{row['accession']}\n{seq}\n")
    return path, rows

def main():
    mafft, fasttree = commands()
    input_path, rows = write_input()
    alignment_path = OUT / "gh2_mafft.fasta"
    with alignment_path.open("w", encoding="ascii") as f:
        subprocess.run([mafft, "--auto", "--inputorder", str(input_path)],
                       stdout=f, stderr=subprocess.PIPE, text=True, check=True)
    alignment = AlignIO.read(alignment_path, "fasta")
    observed = {r.id for r in alignment}
    expected = {r["accession"] for r in rows}
    if observed != expected or len(alignment) != len(expected):
        raise ValueError("MAFFT output does not match expected accessions")
    tree_path = OUT / "gh2_fasttree.nwk"
    with tree_path.open("w", encoding="ascii") as f:
        subprocess.run([fasttree, "-wag", str(alignment_path)],
                       stdout=f, stderr=subprocess.PIPE, text=True, check=True)
    tree = Phylo.read(tree_path, "newick")
    if {x.name for x in tree.get_terminals()} != expected:
        raise ValueError("FastTree tips do not match sequences")
    versions = {}
    for cmd, tool in (([mafft, "--version"], "mafft"), ([fasttree, "-help"], "fasttree")):
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        versions[tool] = (result.stdout + result.stderr).splitlines()[:2]
    summary = {
        "protein_count": len(alignment),
        "alignment_columns": alignment.get_alignment_length(),
        "tools": versions,
        "method": "MAFFT --auto --inputorder; FastTree -wag",
        "limitations": [
            "Approximate ML protein gene tree, not the microbial species tree.",
            "No root or outgroup; no claim of molecular-clock direction.",
            "Close lineages can remain oversampled.",
            "FastTree local support values are not standard bootstrap replicates.",
            "Formal phylogenetically controlled thermal effects are not estimated."
        ]
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
