# Kefir4All resubmission -- clean snapshot

Assembled 2026-11-05 from `kefir4all_review_history/` as a single, de-duplicated
working copy: one version of each script (the most modern available), the
current canonical figure outputs, and public data only. This is a local-only
snapshot -- nothing here has been pushed anywhere, and it should not be until
the "NEEDS REVIEW" items below are resolved.

## How this folder was assembled

- **`scripts/`** -- merged from three generations that existed side by side in
  `kefir4all_review_history/` (`scripts/`, `scripts_reviewed/`,
  `scripts_reviewed_and_figure_improved/`). For each script, the most modern
  available version was kept; only 13 scripts had been through the full
  "figure_improved" review pass (see table below) -- everything else is
  carried over unchanged from the original/base generation. Internal path
  references (`here::here(...)` calls) were updated so every script still
  resolves correctly inside this single merged tree (old
  `scripts_reviewed_and_figure_improved` / `figures_ai_skill` path references
  were rewritten to `scripts` / `figures`). All touched scripts were
  syntax-checked (`parse()`) after the rewrite.
- **`figures/`** -- copied from `figures_ai_skill/`, the canonical, most
  current output location (per-manuscript-figure folders, with
  `save_figure()` provenance JSON where available).
- **`data/`** -- copied from `data/`, **excluding `data/private/`**
  (participant-identifying citizen-scientist metadata -- several scripts
  already fail gracefully with a clear message if it's absent, by design).
  **Do not add private data back into this folder if it is ever pushed
  publicly.**
- **`docker/`** -- copied from the last clean git commit (`HEAD`), not the
  working tree, because `kefir4all_review_history/docker/*` had uncommitted,
  unexplained edits that were never reviewed this pass.

## NEEDS REVIEW before this goes anywhere further

1. **Manuscript caption edit pending**: Supplementary Figure 4's Mantel
   result was corrected this pass (R=0.424, p=0.0001, n=437, superseding the
   previously captioned R=0.395 -- see `scripts/r_scripts/04_taxonomic_profiling/04_supp_mantel_gower.R`
   for the full root-cause writeup). The actual manuscript document
   (`Kefir4All_Manuscript_TRACKED_SB.docx`) caption has **not** been updated
   to match -- someone needs to do that by hand or ask for it explicitly.
2. **`07_metabolomics.R` / `07_supp_metabolomics.R`**: never been through any
   review pass. `07_supp_metabolomics.R` is confirmed to write misplaced
   strain-profiling content (`Figure_11.tiff`, `instrain_strain_clusters_cs_v3.tiff`)
   that duplicates Figure 8's analysis and has nothing to do with
   metabolomics -- needs a decision on scope before cleanup.
3. **`05_resistome.R`**: missing a `set.seed()` call -- flagged, not fixed.
4. **Figures 2, 5, 7, 8 and Supplementary 1/5**: carried an existing "matches
   manuscript" status from before this review pass, but were **not**
   independently re-verified against the manuscript images this time.
5. **Figure 5**: produces a top-level `Figure_5.tiff` plus separate
   `Figure_5A/B/C` and `Figure_5_combined` outputs -- never investigated
   whether all of these are intentional component outputs or redundant
   clutter (same question Figure 4 had, which was resolved by keeping only
   the final combined figure).
6. **Extended Data Figure 2** (`06_supp_strainphlan_ani.R`): confirmed
   unrecoverable with data currently in this repo -- needs the original
   StrainPhlAn4 `.mutation` files and a `marker_names.csv` lookup from
   elsewhere.
7. **Supplementary Figure S6A** (StrainPhlAn4 `.polymorphic` input) and
   **S6B/S6E** (`03_drep_cluster_breakdown.R`, `06_supp_strain_profiling.R`):
   S6A has no source data; S6B/S6E have outputs but were never compared
   against the manuscript's images (no embedded images exist to compare
   against in the current manuscript draft for these).
8. **Supplementary Table S1** (`03_mag_taxonomy.R`): not re-verified.
9. **`docker/` changes left out of this snapshot** (see above) were never
   reviewed -- check whether they were intentional before deciding whether to
   bring them in.
10. One known large file: `data/functional/superfocus_p3.csv` is 68MB --
    under GitHub's 100MB hard limit but above its 50MB soft-warning
    threshold; consider Git LFS if this ever gets pushed.

## Script -> figure mapping and status

| Figure | Script | Status |
|---|---|---|
| Figure 2 | `04_taxonomic_profiling.R` | Not re-checked this pass (carried-over match) |
| Figure 3 | `community_changes_v3.R` | Verified this pass -- trimmed to final output only, Docker-run confirmed |
| Figure 4 | `04_community_types.R` | Verified this pass -- dead legacy code removed, trimmed to final output only, Docker-run confirmed |
| Figure 5 | `04_environmental_microbes_v2.R` | Not re-checked this pass; see NEEDS REVIEW #5 |
| Figure 6 | `06_figure6_method_comparison.R` | Size discrepancy vs. reference investigated and explained (rendering-environment noise, not a bug); **not yet pushed to origin per prior note** |
| Figure 7 | `06_strain_profiling.R` | Mislabeling in README fixed; figure content itself not independently re-verified |
| Figure 8 | `06_instrain_temporal_alluvial.R` | Not re-checked this pass (carried-over match) |
| Extended Data Figure 1 | `04_community_stability_v2.R` | Verified this pass -- content scope confirmed (Candidate A), long-standing rendering bug root-caused and fixed, Docker-run confirmed |
| Extended Data Figure 2 | `06_supp_strainphlan_ani.R` | Confirmed blocked -- missing source data |
| Extended Data Figure 3 | `04_supp_taxonomic_profiling.R` | Confirmed free of dead code; not visually re-verified against the manuscript image |
| Supplementary Figure 1 | `04_supp_superfocus_t0_boxplot.R` | Not re-checked this pass (carried-over match) |
| Supplementary Figure 2 | `05_resistome.R` | Carried-over match; missing `set.seed()` (see NEEDS REVIEW #3) |
| Supplementary Figure 3 | `07_metabolomics.R` | Known mismatch (PCoA labels/legends differ); never reviewed |
| Supplementary Figure 4 | `04_supp_mantel_gower.R` | Verified and corrected this pass -- R=0.424 (see NEEDS REVIEW #1 for the manuscript caption) |
| Supplementary Figure 5 | `06_supp_instrain_pointrange.R` | Not re-checked this pass (carried-over match) |
| Supplementary Figure S6A | StrainPhlAn4 `.polymorphic` input | Confirmed blocked -- no source data found |
| Supplementary Figure S6B | `03_drep_cluster_breakdown.R` | Output exists, never compared to manuscript image |
| Supplementary Figure S6E | `06_supp_strain_profiling.R` | Output exists, never compared to manuscript image |
| Supplementary Table S1 | `03_mag_taxonomy.R` | Not re-verified |
| -- | `07_supp_metabolomics.R` | Confirmed to contain misplaced strain-profiling content (see NEEDS REVIEW #2) |

## Running a script

Scripts expect to run inside the project Docker image
(`ghcr.io/liamhwalsh/kefir4all:4.4.2-slim`), from this folder as the mounted
project root, e.g.:

```
docker run --rm -v "$(pwd):/home/rstudio/kefir4all" -w /home/rstudio/kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2-slim \
  /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R
```

Each script's own project-root anchor block expects a `scripts/` directory to
exist alongside it at the mounted root -- this folder's layout already
satisfies that.
