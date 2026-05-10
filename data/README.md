# data/

Numerical source data inputs to each figure. Folder names follow the **original-submission** figure numbering for backward compatibility with the R scripts in `scripts/r_scripts/`. The mapping to revised-manuscript figures is in the top-level README.

| Folder | Used by |
|---|---|
| `Figure_2_data/` | `scripts/r_scripts/Figure_2.R` |
| `Figure_3_data/` | `scripts/r_scripts/Figure_3.R` |
| `Figure_4_data/` | `scripts/r_scripts/Figure_4.R` |
| `Figure_5_data/` | `scripts/r_scripts/Figure_5.R` |
| `Figure_6_data/` | `scripts/r_scripts/Figure_6.R` |
| `Figure_7_data/` | `scripts/r_scripts/Figure_7.R` |
| `Figure_8_data/` | `scripts/r_scripts/Figure_8_legacy.R`, `scripts/python_scripts/render_figure_8.py`, `scripts/python_scripts/render_supplementary_fig_8.py` |
| `Figure_9_data/` | `scripts/r_scripts/Figure_9_legacy.R` |
| `Figure_10_data/` | `scripts/r_scripts/Figure_10_legacy.R`, `scripts/python_scripts/render_supplementary_fig_5.py`, `scripts/python_scripts/render_supplementary_fig_9.py` |

Larger files (e.g. `ndb.csv.gz`, `mag_metadata_v3.csv.gz`) are stored gzipped to keep clone size manageable.
