# Reproducibility with Docker

The container has R 4.4.2 and all 553 packages — no installation needed.

## Run a script directly (no RStudio)

You need the repo's `data/` and `scripts/` directories. Everything else in
the repo (manuscript, logs, figures, Docker build files) is irrelevant.

```bash
# Clone only what you need
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

# Figure appears at: output/04_taxonomic_profiling/Figure_2.jpeg
```

## Run with RStudio (interactive)

```bash
docker run --rm -p 8787:8787 \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  -e PASSWORD=kefir4all \
  ghcr.io/liamhwalsh/kefir4all:4.4.2
```

Open **http://localhost:8787** — login with:
- **Username:** `rstudio`
- **Password:** `kefir4all`

Open any `.R` file in `scripts/r_scripts/` and click **Source**.

## What's inside

| Component | Version |
|---|---|
| R | 4.4.2 |
| RStudio Server | 2024.12.1 |
| CRAN + Bioconductor | 553 packages (2025-06-15 snapshot) |

Full package manifest: [`package-manifest.txt`](package-manifest.txt)

## Troubleshooting

| Problem | Fix |
|---|---|
| Port in use | `-p 8788:8787` |
| Permission errors | Add `-u $(id -u):$(id -g)` |
| Missing packages | Report to authors — all 553 are pre-installed |