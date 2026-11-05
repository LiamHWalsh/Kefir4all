# scripts/

## Folder structure

```
scripts/r_scripts/
├── 03_mag_classification/     Pipeline step 3 — MAG binning, GTDB taxonomy, dRep clustering
├── 04_taxonomic_profiling/    Pipeline step 4 — MetaCache/MetaPhlAn, community analysis, metadata
├── 05_functional_profiling/   Pipeline step 5 — SUPER-FOCUS, resistome
├── 06_strain_profiling/       Pipeline step 6 — inStrain, StrainPhlAn4, strain clusters
├── 07_metabolomics/           Pipeline step 7 — VOC / volatilome analysis
├── supplementary_notes/       Supplementary Notes 6 and 7
├── ena_submission/            ENA submission helper
└── researchgate_legacy/       Original-submission Figure scripts (Figure_2–8.R), archived for citation traceability
```

## Script-to-results cross-reference

**Status as of 2026-11-04.** The table below was rebuilt by checking each
script's actual code/output against the current manuscript
(`Kefir4All_Manuscript_TRACKED_SB.docx`) rather than carried forward from
assumption — the previous version of this table had drifted out of date on
several rows (wrong script, wrong figure number, or a stale output filename
left over from an earlier numbering scheme) after the manuscript was
revised. Status column: ✅ confirmed against the manuscript's actual
caption text/content, ⚠️ exists but has a known problem, ❌ no script
currently produces this figure at all.

| Script | Figure | Status |
|--------|--------|--------|
| `04_taxonomic_profiling/04_taxonomic_profiling.R` | Figure 2 | ✅ |
| `04_taxonomic_profiling/community_changes_v3.R` | Figure 3 | ✅ |
| `04_taxonomic_profiling/04_community_types.R` | Figure 4 | ✅ |
| `04_taxonomic_profiling/04_environmental_microbes_v2.R` | Figure 5 | ✅ |
| `06_strain_profiling/06_figure6_method_comparison.R` | Figure 6 | ✅ ported from a pre-repo legacy script; verified. Current output is ~16.8KB larger than the `figures_for_review/Figure_6.tiff` reference, but TIFF header tags (dimensions 4800x3000, 8-bit RGB, LZW compression, 300dpi) are identical between the two -- the size gap is consistent with LZW-compression-ratio noise from a different rendering environment (R/Cairo/font-hinting version), not a content or data change. **Not yet pushed to origin — pending review.** |
| `06_strain_profiling/06_strain_profiling.R` | Figure 7 | ✅ output file is `Figure_7.tiff`; "Figure_9_model" is a stale stem from an older manuscript draft, confirmed via git history (commit ee124e1) to have always been this same panel, not a separate figure |
| `06_strain_profiling/06_instrain_temporal_alluvial.R` | Figure 8 | ✅ single panel. A second, dRep-based panel in `06_supp_strain_profiling.R` was deliberately dropped, not combined |
| `04_taxonomic_profiling/04_community_stability_v2.R` | — | ⚠️ Content (Bray-Curtis within-stage convergence) is not referenced anywhere in the current manuscript. Panel-title and rendering bugs were fixed, but whether to keep this script is still an open decision |
| `04_taxonomic_profiling/04_supp_taxonomic_profiling.R` | Extended Data Figure 3 | ✅ (output file still named `Supplementary_Fig_5.tiff` — stale name, not yet renamed) |
| `04_taxonomic_profiling/04_community_stability_v2.R` (downstream block) | Extended Data Figure 1 | ⚠️ Runs without error but never produces its output file (`Supplementary_Figure_2.tiff`, also a stale name) — cause not fully diagnosed |
| `06_strain_profiling/06_supp_strainphlan_ani.R` | Extended Data Figure 2 | ❌ **Confirmed unrecoverable as-is.** Blocked on two missing inputs: `data/marker_names.csv` (a clade-ID-to-species lookup, never existed in this repo) and, more fundamentally, the raw per-pair StrainPhlAn4 `.mutation` files the SNP-divergence values are computed from -- these don't exist anywhere in this repo, its git history, or the old pre-repo personal-workstation folders checked. Unlike Figure 6 (which only needed detection *counts*, recoverable from tree tip labels), this figure needs the actual SNP values, which cannot be reconstructed from anything else in the repo. Would require re-running StrainPhlAn4's mutation-calling step or locating the original output elsewhere. |
| `04_taxonomic_profiling/04_supp_superfocus_t0_boxplot.R` | Supplementary Figure 1 | ✅ Visually matches the embedded image and caption; review copy is identical to generated output |
| `05_functional_profiling/05_resistome.R` | Supplementary Figure 2 | ✅ Four-panel output and updated caption match the embedded image (output file: `resistome_figure.tiff`) |
| `07_metabolomics/07_metabolomics.R` | Supplementary Figure 3 | ⚠️ Freshly rendered after source edits; PCoA labels/legends still differ from the reference, so the review copy was not replaced |
| `04_taxonomic_profiling/04_supp_mantel_gower.R` | Supplementary Figure 4 | ✅ Two-panel render matches embedded image; caption corrected to current result (R=0.424, p=0.0001, n=437 -- supersedes an earlier 0.395 figure computed before `data/kefir4all_sample_metadata_v2.csv` was corrected to drop the known conflicting sample `TG_FG_WLC_S410`, same exclusion `community_changes.R`/`community_changes_v2.R` already apply) |
| `06_strain_profiling/06_supp_instrain_pointrange.R` | Supplementary Figure 5 | ✅ Visually matches embedded image |
| StrainPhlAn4 `.polymorphic` input | Supplementary Figure S6A | ❌ No matching source or input data found |
| `03_mag_classification/03_drep_cluster_breakdown.R` | Supplementary Figure S6B | ⚠️ Output exists; no embedded S6B image for visual comparison |
| `06_strain_profiling/06_supp_strain_profiling.R` | Supplementary Figure S6E | ⚠️ `Supplementary_Fig_9.tiff` contains dRep secondary-cluster distributions over time; no embedded image for visual comparison |
| `03_mag_classification/03_mag_taxonomy.R` | Supp Table S1 | ❓ not re-verified this pass |
| `04_taxonomic_profiling/04_supp_community_stability.R` | — | ❓ not re-verified; actual output file (`Figure_9_v2.tiff`) does not match the old claimed "Supp Fig 6" |
| `04_taxonomic_profiling/04_supp_community_types.R` | Legacy Supplementary Figure 8 | InStrain polymorphism/nucleotide-diversity analysis; not current S6A, which requires StrainPhlAn4 `.polymorphic` data |
| `07_metabolomics/07_supp_metabolomics.R` | ? | ❓ Output (`Figure_11.tiff`) appears to contain leftover/duplicated strain-profiling code unrelated to metabolomics — flagged for review, not confirmed |
| `supplementary_notes/Supplementary_Note_6_analysis.R` | Supplementary Note 6 | — |
| `supplementary_notes/Supplementary_Note_7_analysis.R` | Supplementary Note 7 | — |

## unix_scripts/

Pipeline scripts documenting the upstream metagenomics workflow (Trim Galore,
decontamination, MetaCache, MetaPhlAn, SUPER-FOCUS, dRep, inStrain,
StrainPhlAn 4, antimicrobial-resistance pipeline, ENA submission). These are
documentation rather than turn-key scripts; consult the manuscript Methods
section for full parameters.
