# Dev Container for Kefir4All

## One-click setup for reviewers

Open this folder in a dev container (VS Code: Remote-Containers extension):

1. Open VS Code
2. `Dev Containers: Open Folder in Container` (or the blue status bar prompt)
3. VS Code will pull `ghcr.io/liamhwalsh/kefir4all:4.4.2-slim` and configure R
4. The workspace folder is mounted at `/workspaces/<folder-name>` (e.g. `/workspaces/kefir4all`)

## Authentication

The container image is **private** on GitHub Container Registry. Before first run:

```bash
docker login ghcr.io -u YourGitHubUsername -p ghp_yourPATtoken
```

The PAT needs `read:packages` scope. Generate one at:
https://github.com/settings/tokens

## Using RStudio

Port 8787 is forwarded automatically. Once inside the container, start RStudio:

```bash
/usr/local/bin/start.sh
```

Then open `http://localhost:8787` and log in with:
- Username: `rstudio`
- Password: `kefir4all`

## Running scripts

In the VS Code terminal (inside the container):

```bash
# R is available at /opt/R/4.4.2/bin/R
Rscript scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R
```

> Private data files in `data/private/` are required by some scripts.
> See `data/private/README.md` for details.

## Extensions

The following VS Code extensions are pre-installed:
- `Ikuyadeu.r` — R language tools
- `ms-vscode.vscode-R` — R language support
