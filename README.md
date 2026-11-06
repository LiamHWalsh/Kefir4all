# Kefir4All resubmission

Analysis code, figures, and public data supporting the Kefir4All manuscript
resubmission.

## Structure

- **`scripts/`** -- R analysis and figure-generation scripts, organized by
  pipeline stage (`03_mag_classification/`, `04_taxonomic_profiling/`,
  `05_functional_profiling/`, `06_strain_profiling/`, `07_metabolomics/`).
- **`figures/`** -- generated figure outputs, one folder per manuscript
  figure, in TIFF and PDF.
- **`data/`** -- public input data used by the scripts. Does not include
  participant-identifying citizen-scientist metadata.
- **`docker/`** -- the Docker image definition the scripts are run in.

## Script -> figure mapping

| Figure | Script |
|---|---|
| Figure 2 | `04_taxonomic_profiling.R` |
| Figure 3 | `community_changes_v3.R` |
| Figure 4 | `04_community_types.R` |
| Figure 5 | `04_environmental_microbes_v2.R` |
| Figure 6 | `06_figure6_method_comparison.R` |
| Figure 7 | `06_strain_profiling.R` |
| Figure 8 | `06_instrain_temporal_alluvial.R` |
| Extended Data Figure 1 | `04_community_stability_v2.R` |
| Extended Data Figure 3 | `04_supp_taxonomic_profiling.R` |
| Supplementary Figure 1 | `04_supp_superfocus_t0_boxplot.R` |
| Supplementary Figure 2 | `05_resistome.R` |
| Supplementary Figure 3 | `07_metabolomics.R` |
| Supplementary Figure 4 | `04_supp_mantel_gower.R` |
| Supplementary Figure 5 | `06_supp_instrain_pointrange.R` |
| Supplementary Figure S6B | `03_drep_cluster_breakdown.R` |
| Supplementary Figure S6E | `06_supp_strain_profiling.R` |
| Supplementary Table S1 | `03_mag_taxonomy.R` |

## Figure 3 statistics

`community_turnover_mixed_v1.R` tests whether within-participant taxonomic
turnover between consecutive visits differs by interval, per sample type,
without assuming where any change occurs. It fits turnover as a function of
a 5-level interval factor (no breakpoint), with a planned contrast comparing
the first interval to the mean of the later four.

## Running a script

Scripts run inside the project Docker image
(`ghcr.io/liamhwalsh/kefir4all:4.4.2-slim`), from this folder as the mounted
project root:

```
docker run --rm -v "$(pwd):/home/rstudio/kefir4all" -w /home/rstudio/kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2-slim \
  /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R
```
