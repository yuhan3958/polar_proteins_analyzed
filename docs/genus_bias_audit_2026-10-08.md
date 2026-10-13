# Research note — genus imbalance and thermal-label audit (2026-10-08)

## Metadata-only audit (current master snapshot)
- 20 sequences: 10 cold, 10 mesophilic.
- Cold: Arthrobacter 6, Pseudoalteromonas 2, Psychromonas 1, Photobacterium 1 (4 genera).
- Meso: 10 distinct genera, 1 sequence each.
- Reported growth optimum missing: 8/10 cold and 10/10 meso.
- Cold classifications explicitly mentioning psychrotolerant or psychrotrophic: 6/10.

These counts are derived from the *metadata table*, not yet from a full rerun of protein analyses. They do not quantify a temperature effect.

## New script
`python scripts/audit_phylogeny.py` writes:
- `dataset_results/audit/genus_counts.csv`
- `dataset_results/audit/group_sensitivity.csv`: amino-acid G/P/R differences computed both per-protein and after averaging **within genera first**.
- `dataset_results/audit/leave_one_genus_out.csv`: within-group leave-one-genus-out means.
- `dataset_results/audit/classification_flags.csv`: missing optimum, thermal source or differing cold terminology.
- `dataset_results/audit/summary.json`

No new packages required.

## Interpretation rules
1. Equal-genus weighting is **not** a true phylogenetic correction; it is a diagnostic for sampling imbalance.
2. Growth temperatures, culture recommendations, and *enzyme* activity optima must never be conflated.
3. The 190 pairwise alignments generated from 20 sequences are not 190 independent specimens.
4. If the sign or magnitude of the G/P/R contrast changes after genus balancing, flag it as fragile and avoid inference of a universal thermal-adaptation mechanism.
5. Even if signs agree, sampling and taxonomy may still explain the patterns.

## Next experimental decision
- Inspect `classification_flags.csv` and validate temperature evidence for each record.
- Compare ordinary vs genus-balanced G/P/R results.
- Seek additional, taxonomically balanced GH2 sequences with documented growth temperature, or restrict conclusions to descriptive analyses.
- A credible phylogenetic comparison later requires homologous sequences, quality-checked multiple alignment and an explicit phylogenetic tree/model.

## Status
The code is proposed for review. Metadata counts have been inspected independently, but full execution of the new script against all repository FASTA files has **not** been verified in this session.


## Verified CI findings (2026-10-08)
GitHub Actions job `audit-tests` (run 37775986706) completed successfully: 3 unit tests, 20-sequence data audit, and artifact archival.

Values below are cold − meso differences in **percentage points**:

| Residue | Protein-weighted difference | Equal-genus-weighted difference |
|---|---:|---:|
| G | +1.359 | +0.033 |
| P | +0.116 | −0.565 |
| R | −0.284 | −1.433 |

**Scientific interpretation:** the G contrast nearly disappears after equal-genus weighting, P reverses direction, and R becomes more negative. This demonstrates that the descriptive results are sensitive to taxonomic sampling. It does **not** establish a causal temperature adaptation signal nor phylogenetically independent replicates. Use these results in the final discussion with the audit method and its limitations.

CI generated files can be downloaded from the run's `genus-audit-results` artifact.
