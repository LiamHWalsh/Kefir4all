# scripts/

| Folder | Purpose |
|---|---|
| `r_scripts/` | Original R scripts that generated the main-text figures of the initial submission. Filenames retain the **original** figure numbers; some of these became Supplementary Figures or were superseded in the revision. See `_legacy` suffixes and the figure-numbering map in the top-level README. |
| `python_scripts/` | Python scripts added during revision. Each regenerates one current-manuscript figure from real upstream output (no cropping). |
| `unix_scripts/` | Pipeline scripts documenting the upstream metagenomics workflow (Trim Galore, decontamination, MetaCache, MetaPhlAn, SUPER-FOCUS, dRep, inStrain, StrainPhlAn 4, antimicrobial-resistance pipeline, ENA submission). These are documentation rather than turn-key scripts; consult the manuscript Methods section for full parameters. |

## Running the Python figure scripts

Python ≥ 3.9, pandas, numpy, matplotlib. From the repo root:

```bash
python scripts/python_scripts/render_figure_8.py
python scripts/python_scripts/render_supplementary_fig_5.py
python scripts/python_scripts/render_supplementary_fig_8.py
python scripts/python_scripts/render_supplementary_fig_9.py
```

Each writes outputs to `figures/`. None modifies inputs.
