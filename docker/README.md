# Run Kefir4All analyses — zero installation

Skip the 2-hour package install. The container has R 4.4.2, RStudio Server,
and all **553 packages** pre-installed and version-locked. You just pull,
run, and open a browser tab.

---

## Quick start (30 seconds)

```bash
# 1. Pull the image from GitHub Container Registry
docker pull ghcr.io/liamhwalsh/kefir4all:4.4.2

# 2. Run it
docker run --rm -p 8787:8787 -e PASSWORD=kefir4all ghcr.io/liamhwalsh/kefir4all:4.4.2
```

Open **http://localhost:8787** in your browser and log in:
- **Username:** `rstudio`
- **Password:** whatever you set with `PASSWORD=` (above: `kefir4all`)

That's it. RStudio opens and the repo scripts are fully executable — no
`install.packages()`, no compilation, no dependency hunting.

---

## What is inside the container

| Component | Version | Notes |
|---|---|---|
| R | 4.4.2 (2024-10-31) | Installed at `/opt/R/4.4.2` |
| RStudio Server | 2024.12.1+563 | Web IDE; no local install needed |
| CRAN packages | 2025-06-15 snapshot | 500+ packages, all pre-built binaries |
| Bioconductor | 3.20 | phyloseq, DESeq2, ggtree, ComplexHeatmap, DECIPHER, edgeR, ALDEx2, … |
| Key figure packages | pre-installed | ggplot2, ggpubr, ggstatsplot, ggalluvial, ggmsa, ggtreeExtra, vegan, caret, … |

Full package list: [`package-manifest.txt`](package-manifest.txt)
Exact R session: [`sessionInfo.txt`](sessionInfo.txt)

---

## Running the analysis scripts

### Option A — inside RStudio (recommended)

Open the **Terminal** tab in RStudio and run any script:

```bash
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
Rscript scripts/r_scripts/05_functional_profiling/05_resistome.R
# … any script from REPRODUCIBILITY.md
```

Or use the **Files** pane → navigate to `scripts/r_scripts/` →
double-click a `.R` file → click **Source** to run it.

### Option B — one-liner from the host

Mount the repo into the container and run scripts directly:

```bash
docker run --rm \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  kefir4all:4.4.2 \
  Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
```

### Option C — run everything at once

```bash
docker run --rm \
  -v $(pwd):/home/rstudio/kefir4all \
  -w /home/rstudio/kefir4all \
  kefir4all:4.4.2 \
  bash -c 'for s in $(find scripts/r_scripts -name "*.R" | sort); do Rscript "$s" || echo "FAILED: $s"; done'
```

---

## Mounting your own data

The container's home directory is `/home/rstudio`. To work with local files:

```bash
docker run --rm -p 8787:8787 \
  -e PASSWORD=kefir4all \
  -v /path/to/your/data:/home/rstudio/mydata \
  ghcr.io/liamhwalsh/kefir4all:4.4.2
```

Your data appears at `~/mydata` inside RStudio.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `contrib.url … repos not set` | Repos already pinned in the image; this shouldn't happen |
| Package not found | The 553 packages cover every script. Report the missing package name |
| Port 8787 already in use | Use a different port: `-p 8788:8787` |
| Permission denied writing output | Mounted volumes inherit your host UID; use `-u $(id -u):$(id -g)` |
| RStudio login fails | Check the `PASSWORD=` value you set at `docker run` |

---

## Stopping the container

- **RStudio window**: just close the browser tab; the container keeps running
- **To stop**: press `Ctrl+C` in the terminal where you ran `docker run`
- **`--rm` flag** (shown above) ensures the container is deleted on stop — no leftover state

---

## Image size

~8 GB uncompressed (~7.5 GB compressed). Typical for a full Bioconductor +
tidyverse + phylogenetics stack. The container eliminates all installation
time, so this is a one-time download.

---

## How the image was built (for maintainers)

The image is built from a validated, frozen R library tarball — no
compilation happens during `docker build`. Key build decisions:

1. **R 4.4.2** from Posit's prebuilt `.deb` (not distro R 4.3.x)
2. **CRAN snapshot 2025-06-15** via Posit Package Manager Linux binaries — 
   chosen because its `ggplot2` (3.5.2) satisfies both `ggtree` and 
   `ggalluvial` compatibility requirements
3. **Bioconductor 3.20** (the release paired with R 4.4)
4. **Deriv 4.1.3** pinned — later versions require R 4.5 C++ API

Full build recipe: [`Dockerfile`](Dockerfile)
Build artefacts: `r442.deb`, `rstudio.deb`, `R-library-v2.tar.gz` (all 
regenerable; see `.gitignore`)

### Rebuilding from scratch

```bash
# Download prerequisites (not committed to git)
wget https://cdn.posit.co/r/ubuntu-2404/pkgs/r-4.4.2_1_amd64.deb -O r442.deb
wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.12.1-563-amd64.deb -O rstudio.deb
# Generate R-library-v2.tar.gz from a validated R 4.4.2 environment:
#   cd /opt/R/4.4.2/lib/R
#   tar czf R-library-v2.tar.gz library/

# Build and push
docker build -t kefir4all:4.4.2 .
docker tag kefir4all:4.4.2 ghcr.io/liamhwalsh/kefir4all:4.4.2
docker push ghcr.io/liamhwalsh/kefir4all:4.4.2
```