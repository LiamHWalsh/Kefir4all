# Reproducibility

This document explains how to reproduce every figure in the manuscript
from this repository.

## Software versions

The figures were rendered with:

- **R**, version 4.2 or later (4.4 recommended). Required packages are
  listed in [`install_R_packages.R`](install_R_packages.R) and can be
  installed in one command:
  ```
  Rscript install_R_packages.R
  ```
- **Python**, version 3.9 or later. Required packages are pinned in
  [`requirements.txt`](requirements.txt) at minimum versions:
  ```
  python -m venv .venv
  source .venv/bin/activate          # or .venv\Scripts\activate on Windows
  pip install -r requirements.txt
  ```

## Reproducing the figures

All figure scripts read **only** from `data/` and write **only** to
`figures/`. No script modifies the input data; running them is idempotent.

### Main-text figures

```
Rscript scripts/r_scripts/Figure_2.R
Rscript scripts/r_scripts/Figure_3.R
Rscript scripts/r_scripts/Figure_4.R
Rscript scripts/r_scripts/Figure_5.R
Rscript scripts/r_scripts/Figure_6.R
Rscript scripts/r_scripts/Figure_7.R
Rscript scripts/r_scripts/Figure_8.R          # revised, single inStrain panel
```

### Supplementary figures

```
Rscript scripts/r_scripts/Supplementary_Fig_5.R
Rscript scripts/r_scripts/Supplementary_Fig_6.R   # was original Figure 8
Rscript scripts/r_scripts/Supplementary_Fig_7.R   # was original Figure 9
Rscript scripts/r_scripts/Supplementary_Fig_8.R
Rscript scripts/r_scripts/Supplementary_Fig_9.R
```

Equivalent Python implementations are provided for the four
revision-era figures (Figure 8, Supplementary Figs. 5, 8, 9) under
`scripts/python_scripts/`. Both implementations read the same input
files and produce numerically equivalent outputs; choose whichever
language you prefer.

## Stricter reproducibility (locked package versions)

For exact, version-locked reproducibility we recommend pinning R and
Python packages with environment-management tools:

### R: renv

```
install.packages("renv")
renv::init()
renv::snapshot()    # writes renv.lock
```

Anyone reproducing the analysis later runs `renv::restore()` to install
the exact package versions captured in `renv.lock`.

### Python: pinned requirements

The `requirements.txt` file specifies minimum versions; for exact
reproducibility, pin to specific versions and capture them once your
environment is set:

```
pip freeze > requirements.lock.txt
```

### Container

For the strongest reproducibility guarantee, build a container that
includes both R and Python with the exact package versions used. A
minimal `Dockerfile` could base on `rocker/r-ver:4.4.0`, install
Python 3.11, and run the install scripts on top.

## Long-term archival

When the manuscript is accepted, we recommend:

1. Tag a release on this repository (`git tag v1.0`, `git push --tags`).
2. Archive that release on Zenodo via the GitHub--Zenodo integration so
   it is assigned a permanent DOI that can be cited from the manuscript.
3. Update `CITATION.cff` with the published DOI and the manuscript's
   journal reference.

## Sanity checks

Each script writes a `*_data.tsv` source-data table next to its PNG/PDF
in `figures/`. After running all scripts, you should see, in `figures/`:

- `Figure_8.png`, `Figure_8.pdf`, `Figure_8_data.tsv`
- `Supplementary_Fig_5.png/.pdf/_data.tsv`
- `Supplementary_Fig_8.png/.pdf/_data.tsv`
- `Supplementary_Fig_9.png/.pdf/_data.tsv`
- `Supplementary_Table_1.tsv`

Compare visual outputs to the figures in the manuscript and the
numerical TSVs to those reported in the supplementary information.
