# Reproducibility

This document explains how to reproduce every figure in the manuscript
from this repository.

## Software versions

The figures were rendered with **R**, version 4.2 or later (4.4 recommended).
Required packages are listed in `install_R_packages.R` and can be installed
in one command:

```
Rscript install_R_packages.R
```

## Reproducing the figures

All figure scripts read **only** from `data/` and write **only** to
`figures/`. No script modifies the input data.

### Main-text figures

```
Rscript scripts/r_scripts/Figure_2.R
Rscript scripts/r_scripts/Figure_3.R
Rscript scripts/r_scripts/Figure_4.R
Rscript scripts/r_scripts/Figure_5.R
Rscript scripts/r_scripts/Figure_6.R
Rscript scripts/r_scripts/Figure_7.R
Rscript scripts/r_scripts/Figure_8.R          # revised single inStrain panel
```

### Supplementary figures

```
Rscript scripts/r_scripts/Supplementary_Fig_5.R
Rscript scripts/r_scripts/Supplementary_Fig_6.R   # was original Figure 8
Rscript scripts/r_scripts/Supplementary_Fig_7.R   # was original Figure 9
Rscript scripts/r_scripts/Supplementary_Fig_8.R
Rscript scripts/r_scripts/Supplementary_Fig_9.R
```

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
