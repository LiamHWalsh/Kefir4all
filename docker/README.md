# Reproducibility with Docker

The entire R environment is pre-built into a container — R 4.4.2 and all 553 packages. No installation needed.

## Run scripts without RStudio (command line)

```bash
git clone https://github.com/LiamHWalsh/kefir4all.git
cd kefir4all
docker pull ghcr.io/liamhwalsh/kefir4all:4.4.2

# Run one script
docker run --rm -v $(pwd):/home/rstudio/kefir4all -w /home/rstudio/kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2 \
  /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R

# Figure appears at: output/04_taxonomic_profiling/Figure_2.jpeg
```

## Run with RStudio (interactive)

```bash
docker run --rm -p 8787:8787 -v $(pwd):/home/rstudio/kefir4all \
  -e PASSWORD=kefir4all ghcr.io/liamhwalsh/kefir4all:4.4.2
```

Open **http://localhost:8787** → login: `rstudio` / `kefir4all`. Open any `.R` file in the Files pane and click **Source**.

## What's inside

| Component | Version |
|---|---|
| R | 4.4.2 |
| RStudio Server | 2024.12.1 |
| CRAN packages | 500+ (2025-06-15 snapshot) |
| Bioconductor | 3.20 |

Full manifest: [`package-manifest.txt`](package-manifest.txt)

## Troubleshooting

| Problem | Fix |
|---|---|
| Port in use | `-p 8788:8787` |
| Permission errors | Add `-u $(id -u):$(id -g)` |
| Missing packages | Report to authors — all 553 are pre-installed |