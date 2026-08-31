# Kefir4all

Code and numerical source data for the Kefir4All citizen-science metagenomics study by Walsh et al.

The repository contains R scripts, Unix pipeline documentation, numerical source data, and rendered figures for the manuscript:

> **Adaptations and community changes in milk and water kefir microbiomes in response to environmental parameters as revealed by the Kefir4All Citizen Science Project.**

## Docker (zero install)

All 553 R packages are pre-installed in a container — `docker pull` then `docker run`. See **[docker/README.md](docker/README.md)** for command-line and RStudio instructions.

## Run bare-metal

Install packages with `Rscript install_R_packages.R`, then run any script below. See **[REPRODUCIBILITY.md](REPRODUCIBILITY.md)** for the full list.

```bash
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R   # Figure 2
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_stability.R   # Figure 3
Rscript scripts/r_scripts/07_metabolomics/07_metabolomics.R                 # Figure 4
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R       # Figure 5
Rscript scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes.R # Figure 6
Rscript scripts/r_scripts/06_strain_profiling/06_strain_profiling.R          # Figures 7–8
```

## Figure → script map

> Original-submission figures 7–8 are archived in `researchgate_legacy/` for traceability to the ResearchGate preprint.

| Revised manuscript | Script | What it generates |
|---|---|---|
| Figure 2 | `scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R` | Species composition heatmaps / bar charts |
| Figure 3 | `scripts/r_scripts/04_taxonomic_profiling/04_community_stability.R` | Community stability (Bray-Curtis, Mantel) |
| Figure 4 | `scripts/r_scripts/07_metabolomics/07_metabolomics.R` | VOC / volatilome analysis |
| Figure 5 | `scripts/r_scripts/04_taxonomic_profiling/04_community_types.R` | Community types |
| Figure 6 | `scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes.R` | Environmental microbes |
| Figures 7–8 | `scripts/r_scripts/06_strain_profiling/06_strain_profiling.R` | Strain profiling |
| Figure 8 (inStrain) | `scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R` | inStrain temporal cluster alluvial |
| Supplementary Fig. 5 | `scripts/r_scripts/04_taxonomic_profiling/04_supp_taxonomic_profiling.R` | ANI pairwise |
| Supplementary Fig. 6 | `scripts/r_scripts/04_taxonomic_profiling/04_supp_community_stability.R` | Community stability supp. |
| Supplementary Fig. 7 | `scripts/r_scripts/07_metabolomics/07_supp_metabolomics.R` | Metabolomics supp. |
| Supplementary Fig. 8 | `scripts/r_scripts/04_taxonomic_profiling/04_supp_community_types.R` | Community types supp. |
| Supplementary Fig. 9 | `scripts/r_scripts/06_strain_profiling/06_supp_strain_profiling.R` | dRep strain profiling supp. |
| Supplementary Table 1 | per-species secondary-cluster counts in `figures/Supplementary_Table_1.tsv` | — |
| Original-submission Figures 7–8 | `scripts/r_scripts/researchgate_legacy/Figure_7.R`, `Figure_8.R` | Archived for traceability |

## Data

Raw reads: **PRJEB77409** (European Nucleotide Archive).  
Numerical source data: `data/` (loaded directly by the figure scripts).  
Pre-rendered figures: `figures/`.

This repository contains **no individual citizen-scientist identifiers or re-identifiable metadata**.

## Reproducibility

Full instructions (Docker, bare-metal, script ordering): **[REPRODUCIBILITY.md](REPRODUCIBILITY.md)**.

## Citation

Please cite the manuscript. Update `CITATION.cff` with the DOI on acceptance.

## Licence

MIT — see `LICENSE`.
