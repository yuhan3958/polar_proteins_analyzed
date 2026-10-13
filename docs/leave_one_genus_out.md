# Leave-one-genus-out sensitivity analysis

Run `python scripts/audit_phylogeny.py` and inspect `dataset_results/audit/leave_one_genus_out_contrasts.csv`.

For each taxonomic genus, remove its contribution from one group and recalculate the **full cold-minus-mesophilic** contrast after giving each remaining genus equal weight within its group. Columns show the new contrast, baseline, change, and whether the sign reverses.

This differs from `leave_one_genus_out.csv`, which contains only the within-group average after exclusion. If all genera in one group would be excluded, the operation is skipped.

Interpret large shifts, particularly after excluding *Arthrobacter* from the cold group, as **sensitivity to taxon sampling**, not as an estimate of independent-sample uncertainty or evidence for a temperature effect. Do not interpret these variants as bootstrap confidence intervals or treat 190 protein pairs as 190 independent observations.

The CI workflow runs three regression tests and regenerates the 20-protein audit results. More robust biological inference requires a suitable phylogenetic tree, matched homologous proteins and reliable microbial growth-temperature metadata.
