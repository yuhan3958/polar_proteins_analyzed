# Growth-temperature evidence and structural validation (2026-10-08)

## Source review (strain-specific statements only)

| Accession | Organism | Evidence inspected | Supported observation | Unsupported leap |
|---|---|---|---|---|
| P81650 | *Pseudoalteromonas haloplanktis* TAE79 | https://pmc.ncbi.nlm.nih.gov/articles/PMC92765/ | Antarctic psychrophilic isolate; published beta-galactosidase has high catalytic efficiency at low temperatures; enzyme's apparent activity optimum is approximately 45°C in a later review (https://pmc.ncbi.nlm.nih.gov/articles/PMC7832830/) | Enzyme activity optimum is **not** the organism's optimal growth temperature; no numeric growth optimum confirmed here |
| Q7WTU5 | *Arthrobacter* sp. SB | https://pmc.ncbi.nlm.nih.gov/articles/PMC193751/ | Cells grew from 0–25°C and failed to form colonies at 30°C on the reported media; enzyme activity optimum ~18°C | A single organism-level optimum was **not** measured by these statements |
| A1SWB8 | *Psychromonas ingrahamii* strain 37 | https://link.springer.com/article/10.1007/s12033-018-0142-z | Growth optimum 5°C; no growth above 15°C as reported | Does not establish beta-galactosidase activity optimum |
| Q8KRF6 | *Arthrobacter* sp. C2-2 | https://www.rcsb.org/structure/1YQ2 | Experimental 1.90 Å X-ray structure of cold-active GH2 beta-galactosidase isoenzyme C2-2-1 | The crystal structure is not a direct flexibility measurement, and cannot independently establish microbial growth optimum |

The remaining records retain their existing source links but have **not** been independently verified line-by-line in this audit. Blank growth optima must remain blank.

## Phylogeny confidence

The current MAFFT + FastTree protein tree is an exploratory **gene tree**, not a species phylogeny. The FastTree local branch supports are not standard nonparametric bootstrap replicates. Reliability requires:
1. inspecting alignment gaps/domain coverage and excluding badly aligned portions,
2. checking model fit and support via resampling or an independent inference method,
3. comparing with a species tree / taxonomic clades,
4. assessing root/outgroup selection and sensitivity to overrepresented *Arthrobacter* sequences.

The script `scripts/tree_support_audit.py` (if present) should summarize leaf coverage and FastTree local support without inventing bootstrap confidence.

## Experimental structure comparison

- Mesophilic reference: 1DP0, *E. coli* LacZ, X-ray resolution 1.70 Å; https://www.rcsb.org/structure/1DP0
- Cold-group structural example: 1YQ2, *Arthrobacter* sp. C2-2 isoenzyme C2-2-1, X-ray resolution 1.90 Å, reported engineered mutation; https://www.rcsb.org/structure/1YQ2

The script `scripts/structure_distances.py` maps reference catalytic sites Arg388, Glu461, Met502, Tyr503, Glu537 to PDB residue coordinates using sequence alignment and calculates CA and nearest heavy atom distances for each mapped pair. **P81650 has not thereby been assigned an experimentally determined structure.**

Distance differences in static crystals must not be interpreted as measured flexibility, catalytic rates, temperature dependence, or causal cold-adaptation mechanisms.

## Original publication
Hoyoux et al. (2001) https://doi.org/10.1128/AEM.67.4.1529-1535.2001; Coker et al. (2003) https://doi.org/10.1128/JB.185.18.5473-5482.2003; Skalova et al. (2005) https://doi.org/10.1016/j.jmb.2005.08.028
