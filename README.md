# Kefir4all

Code and numerical source data for the Kefir4All citizen-science metagenomics study by Walsh et al.

The repository contains R scripts, Unix pipeline documentation, numerical source data, and rendered figures for the manuscript:

> **Adaptations and community changes in milk and water kefir microbiomes in response to environmental parameters as revealed by the Kefir4All Citizen Science Project.**

## Docker (zero install)

All 553 R packages are pre-installed in a container — `docker pull` then `docker run`. See **[docker/README.md](docker/README.md)** for command-line and RStudio instructions.

## Run bare-metal

Install packages with `Rscript install_R_packages.R`, then run any script below. See **[REPRODUCIBILITY.md](REPRODUCIBILITY.md)** for the full list.

```bash
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R    # Figure 2
Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes_v3.R     # Figure 3
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R       # Figure 4
Rscript scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes_v2.R # Figure 5
Rscript scripts/r_scripts/06_strain_profiling/06_strain_profiling.R          # Figure 7
Rscript scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R # Figure 8
```

Figure 6 currently has no producing script (see map below) — do not run anything expecting it to exist yet.

## Figure → script map

**Status as of 2026-11-04**, verified against the actual current manuscript
(`Kefir4All_Manuscript_TRACKED_SB.docx`) rather than carried forward from the
previous (stale) version of this table — several rows below previously
pointed at the wrong script or a figure number left over from an earlier
draft's numbering. See `scripts/README.md` for the fuller script-by-script
status table, including Extended Data figures and remaining gaps.

| Manuscript figure | Script | What it generates | Status |
|---|---|---|---|
| Figure 2 | `scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R` | Prevalent-species heatmaps | ✅ |
| Figure 3 | `scripts/r_scripts/04_taxonomic_profiling/community_changes_v3.R` | Taxonomic/functional/volatile-profile changes over fermentation | ✅ |
| Figure 4 | `scripts/r_scripts/04_taxonomic_profiling/04_community_types.R` | PCoA + community-type alluvial | ✅ |
| Figure 5 | `scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes_v2.R` | Diversity + environmental-taxa persistence | ✅ |
| Figure 6 | *(none)* | Method-comparison detection counts (dRep/inStrain/StrainPhlAn4) | ❌ no current script produces this |
| Figure 7 | `scripts/r_scripts/06_strain_profiling/06_strain_profiling.R` | Co-detection of secondary-cluster representatives | ✅ |
| Figure 8 | `scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R` | Temporal distribution of secondary clusters | ✅ |
| Extended Data Fig. 1 | `scripts/r_scripts/04_taxonomic_profiling/04_community_stability_v2.R` | Consecutive-visit Bray-Curtis/Gower dissimilarity | ⚠️ runs, but never writes its output file |
| Extended Data Fig. 2 | `scripts/r_scripts/06_strain_profiling/06_supp_strainphlan_ani.R` | StrainPhlAn SNP divergence + dRep ANI | ❌ blocked on a missing input file (`data/marker_names.csv`) |
| Extended Data Fig. 3 | `scripts/r_scripts/04_taxonomic_profiling/04_supp_taxonomic_profiling.R` | Pairwise MAG relatedness (ANI) | ✅ |
| Supplementary Fig. 1 | *(none)* | Functional (SUPER-FOCUS) Bray-Curtis vs T0, grain only | ❌ no current script produces this |
| Supplementary Fig. 2 | `scripts/r_scripts/05_functional_profiling/05_resistome.R` | AMR class prevalence | ✅ |
| Supplementary Fig. 4 | *(none)* | Bray-Curtis vs Gower Mantel test (hexbin) | ❌ no current script produces this |
| Supplementary Fig. 5 | *(none)* | InStrain SNS/SNV/r2 point-range | ❌ no current script produces this |
| Supplementary Table 1 | per-species secondary-cluster counts in `figures/Supplementary_Table_1.tsv` | — | — |
| Original-submission Figures 7–8 | `scripts/r_scripts/researchgate_legacy/Figure_7.R`, `Figure_8.R` | Archived for traceability | — |

Three figures (6, Supplementary 1, 4, 5) do not currently have a producing
script in this repo at all; legacy pre-repo candidate scripts exist outside
this repo and are being ported over as a separate piece of work.

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
