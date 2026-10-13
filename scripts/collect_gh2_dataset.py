"""Download and validate the explicitly selected GH2 protein dataset.

Uses public UniProt/NCBI APIs. Run with --offline to rebuild from cached files.
Thermal classifications are curated evidence, not inferred from sequence.
"""

import argparse
import concurrent.futures
import csv
import hashlib
import io
import json
import sys
import time
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
AA = set("ACDEFGHIKLMNPQRSTVWY")
DOMAINS = {"PF02836", "PF02837", "PF02929", "PF16353"}


def download(url, path, offline):
    if path.exists():
        return path.read_bytes()
    if offline:
        raise FileNotFoundError(f"Missing cached input: {path}")
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "GH2Research/1.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return data
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def acquire(row, offline):
    row = dict(row)
    accession = row["accession"]
    raw = ROOT / "data" / "raw" / row["database"].lower()
    if row["database"] == "UniProt":
        json_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        fasta_url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        entry = json.loads(download(json_url, raw / (accession + ".json"), offline))
        body = download(fasta_url, raw / (accession + ".fasta"), offline).decode("utf-8")
        if entry.get("primaryAccession") != accession or "sequence" not in entry:
            raise ValueError(f"{accession}: inactive or unexpected UniProt record")
        pfam = {r["id"] for r in entry.get("uniProtKBCrossReferences", []) if r["database"] == "Pfam"}
        if not DOMAINS <= pfam:
            raise ValueError(f"{accession}: missing LacZ-like domain annotations {DOMAINS - pfam}")
        names = entry.get("proteinDescription", {})
        row["protein_name"] = names.get("recommendedName", {}).get("fullName", {}).get("value", "")
        if not row["protein_name"]:
            row["protein_name"] = names.get("submissionNames", [{}])[0].get("fullName", {}).get("value", "")
        if "galactosidase" not in row["protein_name"].lower():
            raise ValueError(f"{accession}: function annotation mismatch")
        row["database_organism"] = entry["organism"]["scientificName"]
        if not row["database_organism"].startswith(row["organism"].split(" (strain")[0]):
            raise ValueError(f"{accession}: organism mismatch")
        row["taxon_id"] = entry["organism"]["taxonId"]
        row["domain_ids"] = ";".join(sorted(pfam))
        row["interpro_ids"] = ";".join(r["id"] for r in entry.get("uniProtKBCrossReferences", []) if r["database"] == "InterPro")
        row["reviewed_status"] = entry["entryType"]
        row["sequence_version"] = entry["entryAudit"]["sequenceVersion"]
        row["protein_existence"] = entry.get("proteinExistence", "")
        row["family_evidence"] = "UniProt GH2 annotation and four LacZ-like Pfam domains"
        row["family_source_url"] = json_url
        row["phylogenetic_lineage"] = ";".join(entry["organism"].get("lineage", []))
        expected = entry["sequence"]["value"]
    else:
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
        fasta_url = base + f"&id={accession}&rettype=fasta&retmode=text"
        body = download(fasta_url, raw / (accession + ".fasta"), offline).decode("utf-8")
        gp_url = base + f"&id={accession}&rettype=gp&retmode=text"
        gp_body = download(gp_url, raw / (accession + ".gp"), offline).decode("utf-8")
        gp = SeqIO.read(io.StringIO(gp_body), "genbank")
        expected = str(gp.seq)
        if "20B" not in gp.annotations["organism"] or gp.id != accession:
            raise ValueError(f"{accession}: NCBI organism/accession mismatch")
        row["database_organism"] = gp.annotations["organism"]
        row["protein_name"] = gp.description.split(" [")[0]
        row["phylogenetic_lineage"] = ";".join(gp.annotations.get("taxonomy", []))
        row["family_evidence"] = "CAZy GH2_10 experimentally characterized beta-galactosidase; NCBI original sequence"
        row["family_source_url"] = "https://www.cazy.org/GH2_10_characterized.html"
        row["reviewed_status"] = "NCBI GenBank original submission"
        row["sequence_version"] = accession.split(".")[-1]
    records = list(SeqIO.parse(io.StringIO(body), "fasta"))
    if len(records) != 1:
        raise ValueError(f"{accession}: expected exactly one FASTA record")
    sequence = str(records[0].seq).upper()
    if not sequence or set(sequence) - AA:
        raise ValueError(f"{accession}: empty/nonstandard sequence")
    if sequence != expected or not 950 <= len(sequence) <= 1150:
        raise ValueError(f"{accession}: FASTA/metadata discrepancy or unexpected length")
    row["sequence_length"] = len(sequence)
    row["sequence_sha256"] = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    row["sequence_url"] = fasta_url
    row["enzyme_family"] = "GH2"
    row["validation_status"] = "sequence_family_checked;thermal_evidence_curated"
    row["fasta_path"] = f"data/validated/{row['group']}/{accession}.fasta"
    row["accessed_on"] = datetime.now(timezone.utc).date().isoformat()
    row["growth_optimum_c"] = row.get("growth_optimum_c", "")
    row["enzyme_activity_optimum_c"] = row.get("enzyme_activity_optimum_c", "")
    row["culture_temperature_c"] = row.get("culture_temperature_c", "")
    target = ROOT / row["fasta_path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(body, encoding="utf-8")
    print(f"Validated {row['group']} {accession}: {len(sequence)} aa", flush=True)
    return row, body


def save_csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "data" / "selection.json").read_text(encoding="utf-8"))
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            acquired = list(pool.map(lambda row: acquire(row, args.offline), manifest["records"]))
        rows = [row for row, body in acquired]
        counts = Counter(row["group"] for row in rows)
        if counts != {"cold": 10, "meso": 10}:
            raise ValueError(f"Incorrect group counts: {counts}")
        if len({r["sequence_sha256"] for r in rows}) != 20:
            raise ValueError("Duplicate protein sequences in selected dataset")
        for group in ["cold", "meso", "all"]:
            joined = "".join(body.rstrip() + "\n" for row, body in acquired if group == "all" or row["group"] == group)
            (ROOT / "data" / "validated" / (group + ".fasta")).write_text(joined, encoding="utf-8")
        save_csv(ROOT / "data" / "metadata.csv", rows)
        save_csv(ROOT / "data" / "exclusions.csv", manifest["exclusions"])
        report = {"dataset_version": manifest["dataset_version"], "validated_at_utc": datetime.now(timezone.utc).isoformat(), "group_counts": dict(counts), "unique_sequences": 20, "nonstandard_residues": 0, "length_min": min(r["sequence_length"] for r in rows), "length_max": max(r["sequence_length"] for r in rows), "warnings": ["Cold group includes psychrophiles and psychrotolerants.", "Ten proteins per group do not mean ten independent species in the cold group.", "Some mesophilic classifications rely on culture-collection or species-level evidence; exact growth optima remain unknown.", "A1SWB8 and Q6LL68 enzyme activities are homology annotations, not direct cold-activity assays.", "Phylogeny, high-pressure adaptation and sample selection confound temperature comparisons."]}
        (ROOT / "data" / "dataset_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps(report, ensure_ascii=False, indent=2))
    except Exception as error:
        (ROOT / "data" / "download_error.json").write_text(json.dumps({"error": str(error), "time": datetime.now(timezone.utc).isoformat()}, ensure_ascii=False, indent=2), encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
