# Reproducibility

This document explains how to reproduce the analyses and figures in the
manuscript from this repository.

## Docker — zero-install reproducibility (recommended)

The entire R 4.4.2 environment with all 553 packages is pre-built into a
single container. No local R or package installation is needed.

```bash
# 1. Pull the image (one-time)
docker pull ghcr.io/liamhwalsh/kefir4all:4.4.2

# 2. Start RStudio Server
docker run --rm -p 8787:8787 -e PASSWORD=kefir4all ghcr.io/liamhwalsh/kefir4all:4.4.2

# 3. Open http://localhost:8787 → login: rstudio / kefir4all
# 4. In RStudio's Terminal tab, run any script:
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
```

See [`docker/README.md`](docker/README.md) for:

- Running scripts from the command line without RStudio
- Running all scripts in batch mode
- Mounting your own data into the container
- Troubleshooting common issues

## Bare-metal install

- **R** version 4.2 or later (4.4 recommended). No other language is required.
- Install every package the scripts load — CRAN, Bioconductor, and two
  GitHub-only packages — with one command, from the repository root:

  ```
  Rscript install_R_packages.R
  ```

  On macOS and Windows this uses prebuilt CRAN binaries and takes a few
  minutes. On Ubuntu 22.04/24.04 the script automatically switches to a
  version-pinned Posit Package Manager **binary** snapshot, which keeps the
  install fast and shields it from future CRAN/R version mismatches.

Any script can then be run individually — open it in RStudio and Source it,
or run it from the command line (see below).

## Reproducing the analyses

Each script is standalone: run it from the repository root so that
`here::here()` resolves to the repo root. Scripts read their inputs from
`data/` and write their outputs to `output/<step>/` (directories are created
automatically). The pre-rendered manuscript figures live in `figures/`.
If a script reports a missing input, that dataset could not be redistributed
here; the corresponding figure is included pre-rendered under `figures/`.

### Canonical scripts, in recommended order

```
# Supplementary notes
Rscript scripts/r_scripts/supplementary_notes/Supplementary_Note_6_analysis.R
Rscript scripts/r_scripts/supplementary_notes/Supplementary_Note_7_analysis.R

# 03 — MAG classification
Rscript scripts/r_scripts/03_mag_classification/03_drep_cluster_breakdown.R
Rscript scripts/r_scripts/03_mag_classification/03_mag_novel_species_summary.R
Rscript scripts/r_scripts/03_mag_classification/03_mag_taxonomy.R

# 04 — Taxonomic profiling
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_stability.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_global_comparison.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_taxonomic_profiling.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_community_types.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_community_stability.R

# 05 — Functional profiling
Rscript scripts/r_scripts/05_functional_profiling/05_resistome.R

# 06 — Strain profiling
Rscript scripts/r_scripts/06_strain_profiling/06_strain_profiling.R
Rscript scripts/r_scripts/06_strain_profiling/06_supp_strain_profiling.R
Rscript scripts/r_scripts/06_strain_profiling/06_supp_strainphlan_ani.R
Rscript scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R

# 07 — Metabolomics
Rscript scripts/r_scripts/07_metabolomics/07_metabolomics.R
Rscript scripts/r_scripts/07_metabolomics/07_supp_metabolomics.R

# ENA submission helper
Rscript scripts/r_scripts/ena_submission/ena_submission_r.R
```

To run them all in one go from the repository root:

```
for s in $(find scripts/r_scripts -name '*.R' | sort); do Rscript "$s" || echo "FAILED: $s"; done
```

For example, `04_taxonomic_profiling.R` renders the main-text Figure 2 to
`output/04_taxonomic_profiling/Figure_2.jpeg`.

## Stricter reproducibility (locked package versions)

For exact, version-locked reproducibility we recommend pinning R
packages with `renv`:

```
install.packages("renv")
renv::init()
renv::snapshot()    # writes renv.lock
```

Anyone reproducing the analysis later runs `renv::restore()` to install
the exact package versions captured in `renv.lock`.

## Long-term archival

When the manuscript is accepted, we recommend:

1. Tag a release on this repository (`git tag v1.0`, `git push --tags`).
2. Archive that release on Zenodo via the GitHub--Zenodo integration so
   it is assigned a permanent DOI that can be cited from the manuscript.
3. Update `CITATION.cff` with the published DOI and the manuscript's
   journal reference.

