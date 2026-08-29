# Kefir4All reproducible container — build & usage guide

This directory contains everything needed to reproduce the **exact R
environment** used to run the Kefir4All scripts — baked into a single
Docker image so that nobody ever has to install or compile packages.

---

## 1. What is inside the image

| Component | Version | Source |
|---|---|---|
| OS | Ubuntu 24.04 (noble) | `ubuntu:24.04` base |
| R | **4.4.2** (2024-10-31) | Posit prebuilt `.deb` (`/opt/R/4.4.2`) |
| RStudio Server | 2024.12.1+563 | Posit `.deb` |
| R packages | **550** (manifest: `package-manifest.txt`) | Tarball of the validated library |
| CRAN snapshot | 2025-06-15 (pinned in `Rprofile.site`) | Posit Package Manager binaries |
| Bioconductor | 3.20 (ggtree 3.14.0, phyloseq, ComplexHeatmap, DECIPHER, DESeq2, …) | Bioc release for R 4.4 |

Image size ≈ **7.5 GB** (uncompressed), typical for a Bioconductor-scale stack.
The library is **byte-identical** to the one validated on the build machine:
no compilation happens for the user, ever.

---

## 2. How the environment was built (and why)

1. **R 4.4.2 from a prebuilt Posit binary** (`cdn.posit.co/r/ubuntu-2404/pkgs/r-4.4.2_1_amd64.deb`)
   installed at `/opt/R/4.4.2` — avoids the distro R (4.3.x) entirely.
2. **CRAN packages from a pinned Posit Package Manager snapshot**
   (`…/cran/__linux__/noble/2025-06-15`) — prebuilt Linux binaries, so installs
   are fast *and* version-coherent. That date was chosen deliberately: its
   `ggplot2` (3.5.2) satisfies **both** `ggtree` (needs `check_linewidth`,
   `is_ggplot`) **and** `ggalluvial`-based scripts (need pre-4.0 internals).
   Mixing "current CRAN" with Bioc 3.20 breaks this pairing.
3. **Bioconductor 3.20 packages** via `BiocManager` (R 4.4 → Bioc 3.20 release).
4. **Version-sensitive pins** discovered the hard way:
   - `Deriv 4.1.3` — current CRAN `Deriv` (4.2.x) uses the R-4.5-only
     `R_ClosureFormals` C++ API and **will not compile on R 4.4**. Install
     `Deriv 4.1.3` *before* anything that pulls it (`car` → `doBy`/`msm` →
     `rstatix` → `ggpubr`/`ggstatsplot` → …).
   - `ggplot2 3.5.2`, `ggfun 0.1.8`, `ggforce 0.4.2`, `ggalluvial 0.12.5`,
     `ggnewscale 0.5.1` — the mutually-compatible set (§8 if you rebuild).
   - `libuv1-dev`, `libwebp-dev`, `libproj-dev`, `cargo` — system deps that
     only surface when source fallbacks trigger.
5. **Validated by execution**: 8 scripts run end-to-end; the remaining 13 stop
   by design (§6).

### Pitfalls hit during the build (do not re-learn these)

| Symptom | Root cause | Fix |
|---|---|---|
| `contrib.url … repos not set` | Scripts call `install.packages()` inline; `Rscript` has no default repo | `Rprofile.site` sets pinned repos (already in image) |
| `R_ClosureFormals … Error 1` | `Deriv` ≥ 4.1.6 needs R 4.5 API | Pin `Deriv 4.1.3` |
| `uv.h: No such file` | `fs` needs libuv | `apt install libuv1-dev` |
| `is_ggplot not exported` | ggplot2 too old (3.5.1) | ggplot2 ≥ 3.5.2 |
| `check_linewidth not found` | ggplot2 4.0.x removed it | ggplot2 ≤ 3.5.2 (hence the 2025-06-15 snapshot) |
| `object 'gg_par' not found` | ggalluvial/ggnewscale built for ggplot2 4.x | Keep them on the same snapshot |
| `cc1plus terminated` (OOM) | Heavy source compiles on small VMs | Use PPM binaries; `-j1` for TMB/glmmTMB |
| RStudio: `Path to R not specified` | R not at a standard path | `rsession-which-r=/opt/R/4.4.2/bin/R` in `/etc/rstudio/rserver.conf` |
| `library/4.4/…` nested after untar | Tar entries are `./4.4/<pkg>` | Extract with `--strip-components=2` |
| `useradd rstudio` fails in build | RStudio deb already creates that user | Idempotent user creation in Dockerfile |

---

## 3. Build artefacts

| File | Role |
|---|---|
| `Dockerfile` | The full image recipe (this directory) |
| `start.sh` | Entrypoint: sets the RStudio password from `$PASSWORD`, runs `rserver` in the foreground |
| `package-manifest.txt` | `package version` list of all 551 packages baked in |
| `sessionInfo.txt` | Full `sessionInfo()` of the validated build machine |
| *(not committed)* `R-library.tar.gz` | 848 MB gzip of the library — regenerate, don't store (§4) |
| *(not committed)* `r442.deb`, `rstudio.deb` | Re-downloadable installers (URLs in §4) |