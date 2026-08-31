# Reproducibility

This document explains how to reproduce the analyses and figures in the
manuscript from this repository.

## Docker — zero-install reproducibility (recommended)

The entire R 4.4.2 environment with all 553 packages is pre-built into a
single container. No local R or package installation is needed.

### Command line (no RStudio)

```bash
git clone --depth 1 --filter=blob:none --sparse https://github.com/LiamHWalsh/kefir4all.git
cd kefir4all
git sparse-checkout set data scripts
docker pull ghcr.io/liamhwalsh/kefir4all:4.4.2

# Run any script
docker run --rm \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2 \
  /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
```

### With RStudio

```bash
docker run --rm -p 8787:8787 \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  -e PASSWORD=kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2
```

Open http://localhost:8787 in your browser and log in with:
- **Username:** `rstudio`
- **Password:** `kefir4all`  (set with `-e PASSWORD=` at `docker run`)

The `-w` flag makes RStudio open straight into the project folder (`data/` and `scripts/` visible in the Files pane). If not, click **Home** (`~`) in the Files pane, then **kefir4all**.

Open any `.R` file and click **Source**.

Full details in [`docker/README.md`](docker/README.md).

## Bare-metal install

- **R** 4.2+ (4.4 recommended).
- Install all packages from the repo root:

  ```bash
  Rscript install_R_packages.R
  ```

  This uses prebuilt CRAN binaries (macOS, Windows) or a version-pinned Posit Package Manager snapshot (Ubuntu).

## Reproducing the analyses

Each script is standalone: run it from the repository root (`here::here()` resolves to the repo root). Scripts read from `data/` and write to `output/<step>/`.

### Figure scripts

```bash
# Main text
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_stability.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_community_types.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_global_comparison.R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_environmental_microbes.R
Rscript scripts/r_scripts/06_strain_profiling/06_strain_profiling.R
Rscript scripts/r_scripts/06_strain_profiling/06_instrain_temporal_alluvial.R
Rscript scripts/r_scripts/07_metabolomics/07_metabolomics.R

# Supplementary
Rscript scripts/r_scripts/04_taxonomic_profiling/04_supp_*.R
Rscript scripts/r_scripts/05_functional_profiling/05_resistome.R
Rscript scripts/r_scripts/06_strain_profiling/06_supp_*.R
Rscript scripts/r_scripts/07_metabolomics/07_supp_metabolomics.R
```

