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

| Script | Results section | Figure(s) |
|--------|----------------|-----------|
| `03_mag_classification/03_mag_taxonomy.R` | 3.1 (MAG taxonomy, novel species) | Supp Table S1, Supp S6 |
| `03_mag_classification/03_drep_cluster_breakdown.R` | 3.5 (strain sub-clusters) | Supp Fig S6 (dRep) |
| `04_taxonomic_profiling/04_taxonomic_profiling.R` | 3.1 (prevalent species, heatmaps) | Figure 2 |
| `04_taxonomic_profiling/04_community_stability.R` | 3.2 (Bray-Curtis temporal), 3.4 (Gower/Mantel) | Figure 3 |
| `04_taxonomic_profiling/04_community_types.R` | 3.3 (community states) | Figure 5 |
| `04_taxonomic_profiling/04_global_comparison.R` | 3.3 (pan-metagenome comparison) | Figure 5 |
| `04_taxonomic_profiling/04_environmental_microbes.R` | 3.4 (environmental acquisition, filtering) | Figure 6 |
| `04_taxonomic_profiling/04_supp_taxonomic_profiling.R` | Supp Result S1 | Supp Fig 5 |
| `04_taxonomic_profiling/04_supp_community_stability.R` | Supp Result S4 | Supp Fig 6 |
| `04_taxonomic_profiling/04_supp_community_types.R` | Supp Result S4 | Supp Fig 8 |
| `05_functional_profiling/05_resistome.R` | Supp Result S2 (resistome) | Supp Fig S2 |
| `06_strain_profiling/06_strain_profiling.R` | 3.5–3.7 (strain diversity, co-occurrence) | Figures 7–8 |
| `06_strain_profiling/06_instrain_temporal_alluvial.R` | 3.7 (inStrain temporal clusters) | Figure 8 |
| `06_strain_profiling/06_supp_strain_profiling.R` | Supp Result S6 (dRep temporal) | Supp Fig 9 |
| `06_strain_profiling/06_supp_strainphlan_ani.R` | Supp Result S5 (mutation rates, ANI) | Supp Fig S5 |
| `07_metabolomics/07_metabolomics.R` | 3.2 (VOC early changes) | Figure 4 |
| `07_metabolomics/07_supp_metabolomics.R` | Supp Result S3 | Supp Fig 7 |
| `supplementary_notes/Supplementary_Note_6_analysis.R` | Supplementary Note 6 | — |
| `supplementary_notes/Supplementary_Note_7_analysis.R` | Supplementary Note 7 | — |

## unix_scripts/

Pipeline scripts documenting the upstream metagenomics workflow (Trim Galore,
decontamination, MetaCache, MetaPhlAn, SUPER-FOCUS, dRep, inStrain,
StrainPhlAn 4, antimicrobial-resistance pipeline, ENA submission). These are
documentation rather than turn-key scripts; consult the manuscript Methods
section for full parameters.
