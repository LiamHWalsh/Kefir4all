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
| *(none)* | Figure 6 | ❌ No script combines dRep + inStrain + StrainPhlAn4 detection counts as the caption describes. A 2024 pre-repo script with this logic exists outside this repo and needs porting (old `Q:/` network paths, not repo-relative) |
| `06_strain_profiling/06_strain_profiling.R` | Figure 7 | ✅ |
| `06_strain_profiling/06_instrain_temporal_alluvial.R` | Figure 8 | ✅ single panel. A second, dRep-based panel in `06_supp_strain_profiling.R` was deliberately dropped, not combined |
| `04_taxonomic_profiling/04_community_stability_v2.R` | — | ⚠️ Content (Bray-Curtis within-stage convergence) is not referenced anywhere in the current manuscript. Panel-title and rendering bugs were fixed, but whether to keep this script is still an open decision |
| `04_taxonomic_profiling/04_supp_taxonomic_profiling.R` | Extended Data Figure 3 | ✅ (output file still named `Supplementary_Fig_5.tiff` — stale name, not yet renamed) |
| `04_taxonomic_profiling/04_community_stability_v2.R` (downstream block) | Extended Data Figure 1 | ⚠️ Runs without error but never produces its output file (`Supplementary_Figure_2.tiff`, also a stale name) — cause not fully diagnosed |
| `06_strain_profiling/06_supp_strainphlan_ani.R` | Extended Data Figure 2 | ❌ Blocked: requires `data/marker_names.csv`, which has never existed in this repo (an external file, not a committed-then-deleted one) |
| `05_functional_profiling/05_resistome.R` | Supplementary Figure 2 | ✅ (output file still named `resistome_figure.tiff`) |
| `04_taxonomic_profiling/04_supp_superfocus_t0_boxplot.R` | Supplementary Figure 1 | ✅ reuses `community_changes_v3.R`'s own output data; verified |
| `04_taxonomic_profiling/04_supp_mantel_gower.R` | Supplementary Figure 4 | ✅ ported from a pre-repo legacy script; verified (Mantel r=0.395, p=1e-4, n=437) |
| `06_strain_profiling/06_supp_instrain_pointrange.R` | Supplementary Figure 5 | ✅ ported from a pre-repo legacy script; verified |
| `03_mag_classification/03_mag_taxonomy.R` | Supp Table S1 | ❓ not re-verified this pass |
| `03_mag_classification/03_drep_cluster_breakdown.R` | Supp Fig S6 | ❓ not re-verified this pass |
| `04_taxonomic_profiling/04_supp_community_stability.R` | — | ❓ not re-verified; actual output file (`Figure_9_v2.tiff`) does not match the old claimed "Supp Fig 6" |
| `04_taxonomic_profiling/04_supp_community_types.R` | Supplementary Figure 3? | ❓ not re-verified this pass |
| `06_strain_profiling/06_supp_strain_profiling.R` | — | Produces a dRep-based "clusters per visit" panel; not currently used in any figure on its own (see Figure 8 note above) |
| `07_metabolomics/07_metabolomics.R` | ? | ❓ Needs re-verification — its standalone VOC output (`Figure_9.tiff`) may now be superseded by Figure 3's volatile-profile row, which is produced by `community_changes_v3.R` instead |
| `07_metabolomics/07_supp_metabolomics.R` | ? | ❓ Output (`Figure_11.tiff`) appears to contain leftover/duplicated strain-profiling code unrelated to metabolomics — flagged for review, not confirmed |
| `supplementary_notes/Supplementary_Note_6_analysis.R` | Supplementary Note 6 | — |
| `supplementary_notes/Supplementary_Note_7_analysis.R` | Supplementary Note 7 | — |

## unix_scripts/

Pipeline scripts documenting the upstream metagenomics workflow (Trim Galore,
decontamination, MetaCache, MetaPhlAn, SUPER-FOCUS, dRep, inStrain,
StrainPhlAn 4, antimicrobial-resistance pipeline, ENA submission). These are
documentation rather than turn-key scripts; consult the manuscript Methods
section for full parameters.
