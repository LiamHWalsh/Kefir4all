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

Reviewer-revision analysis (not a manuscript figure):
`scripts/r_scripts/04_taxonomic_profiling/community_changes.R` writes shared laboratory
T0 reference distances, all five consecutive participant intervals, a participant-blocked
permutation omnibus (Kruskal-Wallis kept alongside for reference only) with Wilcoxon
follow-ups, a descriptive stage-distance summary with observed-vs-hypothesised peak, and
PERMANOVA, to `output/04_taxonomic_profiling/community_changes/` (see `METHODS.txt`
there), plus the figures below. Run it from the repository root:
`Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes.R`.
Volatile panels use Euclidean distances on measured concentrations; taxonomic and
functional panels use Bray-Curtis distances on row proportions. The same script renders
all ten panels (reference vs interval distance) into
`figures/community_changes_figure.{tiff,pdf}` and the peak-stage distribution into
`figures/community_changes_peak_stage.{tiff,pdf}`, using base R graphics only.
All figure scripts write TIFF at 1800 x 1200 px @ 300 dpi (6 x 4 in, LZW).

MetaCache's species-only profile (`04_metacache_total_species_profile_v2.csv`) keeps only
reads resolved to exactly species rank, discarding reads resolved more specifically
(subspecies/strain) and reads resolved only to genus or above.
`scripts/r_scripts/04_taxonomic_profiling/build_metacache_full_profile.R` rebuilds a full,
mixed-rank profile (`data/04_metacache_total_full_profile.csv`) from the raw per-sample
MetaCache classification reports instead, keeping every `rank::name` feature except
`unclassified`, and excluding samples below a 10,000-read depth floor or with conflicting
study metadata. It needs the raw reports, which are not distributed in this repository;
point `KEFIR4ALL_METACACHE_RAW_DIRS` (colon-separated directories) at them and run:
`Rscript scripts/r_scripts/04_taxonomic_profiling/build_metacache_full_profile.R`.
A full build audit (per-sample read depth, classified fraction, inclusion/exclusion reason,
and raw-file md5) is written to
`output/04_taxonomic_profiling/metacache_full_profile_build_audit.csv`.
`community_changes.R` and `community_changes_v2.R` both load this full profile rather than
the species-only file.

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

## Workspace hygiene (2026-09-19)

One-off probe/debug/test scripts, their outputs, figure-layout experiments,
docker experiment files, and superseded output runs that had accumulated in
the repository root (and in `scripts/r_scripts/04_taxonomic_profiling/`)
during the review were **deleted** (recoverable from the OneDrive recycle
bin). The 12 `v3_*` pipeline run logs remain in `v3_logs/`, and the latest
`community_stability_v3` output run is retained alongside the archived
layer-based approach.

**Unchanged**: all pipeline scripts (including the `v2`/`v3` development
versions and the final `04_community_stability_v3.R`), `data/` (including the
private functional data), `figures/`, the tracked Docker configuration, the
root docs (`README.md`, `REPRODUCIBILITY.md`, `AGENTS.md`, `CITATION.cff`,
`LICENSE`, `.gitignore`, `install_R_packages.R`), and `Pictures/` (personal
OneDrive content — deliberately untouched).
