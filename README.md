# Kefir4all

Code and numerical source data for the Kefir4All citizen-science metagenomics study by Walsh et al.

> **Adaptations and community changes in milk and water kefir microbiomes in response to environmental parameters as revealed by the Kefir4All Citizen Science Project.**

The repository contains:

- **R scripts** (`scripts/r_scripts/`) used to generate every original main-text figure of the manuscript.
- **Unix pipeline scripts** (`scripts/unix_scripts/`) documenting the upstream metagenomics workflow (Trim Galore, decontamination, MetaCache / MetaPhlAn, SUPER-FOCUS, dRep, inStrain, StrainPhlAn 4, etc.).
- **Numerical source data** (`data/`) that the figure scripts read from. Folder names match the original figure numbering used during initial submission; see `scripts/README.md` for a mapping between original-submission numbering and revised-manuscript numbering.
- **Rendered current-revision figures** (`figures/`) plus per-figure source-data TSVs.

## Reproducibility quick start

### R figures

Run any figure script directly from the repo root — each is self-contained:

```r
# Main-text figures
Rscript scripts/r_scripts/Figure_2.R
Rscript scripts/r_scripts/Figure_3.R
Rscript scripts/r_scripts/Figure_4.R
Rscript scripts/r_scripts/Figure_5.R
Rscript scripts/r_scripts/Figure_6.R

# Supplementary analyses (each script produces its own output)
Rscript scripts/r_scripts/04_supp_taxonomic_profiling.R
Rscript scripts/r_scripts/04_supp_community_stability.R
Rscript scripts/r_scripts/04_supp_community_types.R
Rscript scripts/r_scripts/06_supp_strain_profiling.R
Rscript scripts/r_scripts/07_supp_metabolomics.R
```

## Figure-numbering map (original submission to revised manuscript)

> **Original-submission numbering is retained in `researchgate_legacy/`** so that the publicly archived original-submission manuscript on ResearchGate remains traceable to the code that produced its figures. The mapping table below is the authoritative bridge between the two numbering systems.

| Original submission | Revised manuscript | Script(s) |
|---|---|---|
| Figure 2 | Figure 2 | `scripts/r_scripts/Figure_2.R` |
| Figure 3 | Figure 3 | `scripts/r_scripts/Figure_3.R` |
| Figure 4 | Figure 4 | `scripts/r_scripts/Figure_4.R` |
| Figure 5 | Figure 5 | `scripts/r_scripts/Figure_5.R` |
| Figure 6 | Figure 6 | `scripts/r_scripts/Figure_6.R` |
| Figure 7 (original strain overview) | moved to Supplementary Result S6 | `scripts/r_scripts/researchgate_legacy/Figure_7.R` |
| Figure 8 (original inStrain alluvial) | moved to Supplementary Result S6 | `scripts/r_scripts/researchgate_legacy/Figure_8.R` |
| Supplementary Fig. 5 (ANI pairwise) | Supplementary Fig. 5 | `scripts/r_scripts/04_supp_taxonomic_profiling.R` |
| Supplementary Fig. 6 (community stability) | Supplementary Fig. 6 | `scripts/r_scripts/04_supp_community_stability.R` |
| Supplementary Fig. 7 (metabolomics) | Supplementary Fig. 7 | `scripts/r_scripts/07_supp_metabolomics.R` |
| Supplementary Fig. 8 (community types) | Supplementary Fig. 8 | `scripts/r_scripts/04_supp_community_types.R` |
| Supplementary Fig. 9 (dRep strain profiling) | Supplementary Fig. 9 | `scripts/r_scripts/06_supp_strain_profiling.R` |
| n/a | Supplementary Table 1 | per-species secondary-cluster counts in `figures/Supplementary_Table_1.tsv` |

## Data

Raw sequencing reads are deposited at the European Nucleotide Archive under project accession **PRJEB77409**.

The `data/` folder in this repository contains the numerical source data for every figure in the manuscript and supplementary materials. The directory is flat and de-duplicated: each file is named descriptively and is loaded directly by name from the figure scripts in `scripts/`.

## Data protection

This repository contains **no individual citizen-scientist identifiers, no household-level metadata that could be re-identifying, and no participant raw response files**. Aggregated, anonymised metadata used as figure inputs (e.g. `kefir4all_sample_metadata_v2.csv`) is included where it appears as a column input to figure code; participant-identifying CSVs are explicitly excluded by `.gitignore`.

## Reproducibility

For exact reproduction of every figure (R version, Python version, package install, env-locking, container, Zenodo archival), see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Citation

If you use the code or data, please cite the manuscript. Citation metadata is provided in `CITATION.cff` for GitHub's "Cite this repository" button. Update the entry with a DOI on acceptance.

## Licence

MIT — see `LICENSE`.
