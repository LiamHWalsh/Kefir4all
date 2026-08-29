#!/bin/bash
# start.sh — run RStudio Server in the foreground (container main process).
set -e
PW="${PASSWORD:-kefir4all}"
echo "rstudio:${PW}" | chpasswd || true

mkdir -p /var/run/rstudio-server /var/log/rstudio /var/lib/rstudio-server

/opt/R/4.4.2/bin/Rscript -e 'cat("R ready:", R.version.string, "| packages:", length(rownames(installed.packages())), "\n")'

exec /usr/lib/rstudio-server/bin/rserver --server-daemonize=0