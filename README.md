# Kefir4all

Code and numerical source data for the Kefir4All citizen-science metagenomics study by Walsh et al.

> **Adaptations and community changes in milk and water kefir microbiomes in response to environmental parameters as revealed by the Kefir4All Citizen Science Project.**

The repository contains:

- **R scripts** (`scripts/r_scripts/`) used to generate every original main-text figure of the manuscript.
- **Unix pipeline scripts** (`scripts/unix_scripts/`) documenting the upstream metagenomics workflow (Trim Galore, decontamination, MetaCache / MetaPhlAn, SUPER-FOCUS, dRep, inStrain, StrainPhlAn 4, etc.).
- **Numerical source data** (`data/`) that the figure scripts read from. Folder names match the original figure numbering used during initial submission; see `scripts/README.md` for a mapping between original-submission numbering and revised-manuscript numbering.
- **Rendered current-revision figures** (`figures/`) plus per-figure source-data TSVs.

## Reproducibility with Docker

All 553 R packages are pre-installed in a container — no setup needed.

```bash
# Clone only the data and scripts (not the whole repo)
git clone --depth 1 --filter=blob:none --sparse https://github.com/LiamHWalsh/kefir4all.git
cd kefir4all
git sparse-checkout set data scripts
docker pull ghcr.io/liamhwalsh/kefir4all:4.4.2

# Run any script
docker run --rm \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2 \
  /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
```

This generates the figure at `output/04_taxonomic_profiling/Figure_2.jpeg`.

See [`docker/README.md`](docker/README.md) for RStudio mode and troubleshooting.

## Run bare-metal (manual install)

Run any figure script directly from the repo root — each is self-contained:

```bash
# Main-text figures
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R   # Figure 2
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_stability.R   # Figure 3
Rscript scripts/r_scripts/07_metabolomics/07_metabolomics.R                 # Figure 4
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R       # Figure 5
Rscript scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes.R # Figure 6
Rscript scripts/r_scripts/06_strain_profiling/06_strain_profiling.R          # Figures 7–8
Rscript scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R # Figure 8 (inStrain)

# Supplementary analyses
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_taxonomic_profiling.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_community_stability.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_community_types.R
Rscript scripts/r_scripts/05_functional_profiling/05_resistome.R
Rscript scripts/r_scripts/06_strain_profiling/06_supp_strain_profiling.R
Rscript scripts/r_scripts/06_strain_profiling/06_supp_strainphlan_ani.R
Rscript scripts/r_scripts/07_metabolomics/07_supp_metabolomics.R
```

## Figure-numbering map (original submission to revised manuscript)

> **Original-submission numbering is retained in `researchgate_legacy/`** so that the publicly archived original-submission manuscript on ResearchGate remains traceable to the code that produced its figures. The mapping table below is the authoritative bridge between the two numbering systems.

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

Raw sequencing reads are deposited at the European Nucleotide Archive under project accession **PRJEB77409**.

The `data/` folder in this repository contains the numerical source data for every figure in the manuscript and supplementary materials. The directory is flat and de-duplicated: each file is named descriptively and is loaded directly by name from the figure scripts in `scripts/`.

## Data protection

This repository contains **no individual citizen-scientist identifiers, no household-level metadata that could be re-identifying, and no participant raw response files**. Aggregated, anonymised metadata used as figure inputs (e.g. `kefir4all_sample_metadata_v2.csv`) is included where it appears as a column input to figure code; participant-identifying CSVs are explicitly excluded by `.gitignore`.

## Reproducibility

For exact reproduction of every figure (R version, Python version, package install, env-locking, Zenodo archival), see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Citation

If you use the code or data, please cite the manuscript. Citation metadata is provided in `CITATION.cff` for GitHub's "Cite this repository" button. Update the entry with a DOI on acceptance.

## Licence

MIT — see `LICENSE`.
