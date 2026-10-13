"""Audit FastTree local support and taxon distribution without mistaking it for bootstrap."""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path
from Bio import Phylo

ROOT=Path(__file__).resolve().parents[1]
TREE=ROOT/"dataset_results"/"phylogeny"/"gh2_fasttree.nwk"
OUT=ROOT/"dataset_results"/"phylogeny"/"support_audit.json"

def audit(tree_path=TREE,metadata_path=ROOT/"data"/"metadata.csv"):
    tree=Phylo.read(tree_path,"newick")
    with metadata_path.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    taxa={r["accession"]:r for r in rows}
    tips=[x.name for x in tree.get_terminals()]
    if len(set(tips))!=len(tips) or set(tips)!=set(taxa):
        raise ValueError("Tree leaves differ from curated metadata")
    supports=[float(x.confidence) for x in tree.get_nonterminals()
              if x.confidence is not None]
    counts=Counter((r["group"],r["organism"].split()[0]) for r in rows)
    result={
       "leaf_count":len(tips), "internal_nodes":len(tree.get_nonterminals()),
       "nodes_with_support":len(supports),
       "local_support_below_0_70":sum(v<0.70 for v in supports),
       "local_support_below_0_90":sum(v<0.90 for v in supports),
       "local_support_min":min(supports) if supports else None,
       "local_support_max":max(supports) if supports else None,
       "genus_counts":{f"{group}:{genus}":n for (group,genus),n in sorted(counts.items())},
       "caveat":"FastTree SH-like local support values are not bootstrap replicates; tree is a protein gene tree, not an independent species phylogeny."
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,indent=2))
