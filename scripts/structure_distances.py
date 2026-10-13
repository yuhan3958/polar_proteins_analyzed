"""Measure distances among mapped catalytic-site residues in experimental PDB structures.

Run: python scripts/structure_distances.py
Internet required on first run to download RCSB PDB CIF files.
Outputs residue mappings, heavy-atom minimum distances, CA distances,
and a provenance report; NOT protein flexibility or dynamics.
"""
from __future__ import annotations
import csv
import json
import math
import urllib.request
from pathlib import Path
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "dataset_results" / "structures"
SITES = (388, 461, 502, 503, 537)
MODELS = [("P00722", "1DP0"), ("Q8KRF6", "1YQ2")]
RCSB = "https://files.rcsb.org/download/{pdb}.cif"


def fetch_structure(pdb):
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f"{pdb}.cif"
    if not dest.exists():
        request = urllib.request.Request(RCSB.format(pdb=pdb), headers={"User-Agent":"Polar-GH2-research/1.0"})
        with urllib.request.urlopen(request, timeout=45) as response:
            dest.write_bytes(response.read())
    return dest


def read_fasta(accession):
    group = "meso" if accession == "P00722" else "cold"
    path = ROOT / "data" / "validated" / group / (accession + ".fasta")
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected a single sequence")
    return str(records[0].seq).upper()


def map_positions(a, b):
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    pair = aligner.align(a, b)[0]
    coords = [(int(i), int(j)) for i, j in zip(*pair.indices) if i >= 0 and j >= 0]
    matches = sum(a[i] == b[j] for i, j in coords)
    return {i + 1: j + 1 for i, j in coords}, matches / len(coords) if coords else 0.0


def parse_chain(cif):
    d = MMCIF2Dict(str(cif))
    keys = ["_atom_site.label_asym_id", "_atom_site.label_seq_id",
            "_atom_site.label_comp_id", "_atom_site.label_atom_id",
            "_atom_site.type_symbol", "_atom_site.Cartn_x",
            "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "_atom_site.pdbx_PDB_model_num"]
    arrays = [d[k] for k in keys]
    sites = {}
    first_model = arrays[-1][0]
    for chain, number, residue, atom, element, x, y, z, model in zip(*arrays):
        if chain != "A" or number in {".", "?"} or model != first_model:
            continue
        if element.upper() in {"H", "D"}:
            continue
        position = int(number)
        item = sites.setdefault(position, {"residue":seq1(residue, custom_map={"UNK":"X"}), "atoms":{}})
        item["atoms"].setdefault(atom, (float(x), float(y), float(z)))
    if not sites:
        raise ValueError("No chain A atom records found")
    positions = sorted(pos for pos, site in sites.items() if "CA" in site["atoms"])
    sequence = "".join(sites[pos]["residue"] for pos in positions)
    return sites, positions, sequence


def distance(a, b):
    return math.dist(a, b)


def pair_distances(site_a, site_b):
    ca_a, ca_b = site_a["atoms"].get("CA"), site_b["atoms"].get("CA")
    if ca_a is None or ca_b is None:
        raise ValueError("Missing CA")
    minimum = min(distance(a, b) for a in site_a["atoms"].values() for b in site_b["atoms"].values())
    return round(distance(ca_a, ca_b), 3), round(minimum, 3)


def analyze():
    OUT.mkdir(parents=True, exist_ok=True)
    reference = read_fasta("P00722")
    all_pairs, maps, summaries = [], [], []
    for accession, pdb in MODELS:
        raw_seq = read_fasta(accession)
        ref_map, reference_identity = map_positions(reference, raw_seq)
        sites, positions, pdb_seq = parse_chain(fetch_structure(pdb))
        target_map, structure_identity = map_positions(raw_seq, pdb_seq)
        selected = {}
        for ref_pos in SITES:
            target_pos = ref_map.get(ref_pos)
            structure_index = target_map.get(target_pos) if target_pos is not None else None
            label_pos = positions[structure_index - 1] if structure_index is not None else None
            matched = label_pos is not None and raw_seq[target_pos - 1] == sites[label_pos]["residue"]
            if matched and label_pos is not None:
                selected[ref_pos] = sites[label_pos]
            maps.append({"accession":accession,"pdb":pdb,
                         "reference_position":ref_pos,
                         "reference_residue":reference[ref_pos - 1],
                         "target_position":target_pos or "",
                         "pdb_label_seq_id":label_pos or "",
                         "pdb_residue":sites[label_pos]["residue"] if label_pos is not None else "",
                         "residue_matched":matched})
        for ix, a in enumerate(SITES):
            for b in SITES[ix + 1:]:
                if a not in selected or b not in selected:
                    continue
                ca, heavy = pair_distances(selected[a], selected[b])
                all_pairs.append({"accession":accession,"pdb":pdb,
                                  "ref_position_a":a,"ref_position_b":b,
                                  "ca_distance_angstrom":ca,
                                  "nearest_heavy_atom_angstrom":heavy})
        summaries.append({"accession":accession,"pdb":pdb,
                          "reference_to_target_identity":reference_identity,
                          "fasta_to_structure_identity":structure_identity,
                          "mapped_residues":len(selected),
                          "structure_source":f"https://www.rcsb.org/structure/{pdb}"})
    with (OUT / "site_mapping.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        w=csv.DictWriter(stream, fieldnames=list(maps[0]));w.writeheader();w.writerows(maps)
    with (OUT / "distances.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        w=csv.DictWriter(stream,fieldnames=["accession","pdb","ref_position_a","ref_position_b","ca_distance_angstrom","nearest_heavy_atom_angstrom"])
        w.writeheader();w.writerows(all_pairs)
    summary={"structures":summaries, "n_distances":len(all_pairs),
             "limitations":["Mapping based on sequence alignment and first model/label chain A",
                            "CA and heavy-atom distances are static geometry, not flexibility",
                            "Different oligomerization states and crystallization conditions limit comparisons",
                            "1YQ2 contains an engineered mutation and represents C2-2-1 isoenzyme"]}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    analyze()
