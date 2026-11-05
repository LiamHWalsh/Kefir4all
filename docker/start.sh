#!/bin/bash
# start.sh — run RStudio Server in the foreground (container main process).
set -e
PW="${PASSWORD:-kefir4all}"

# Set password: try chpasswd, fall back to usermod
echo "rstudio:${PW}" | chpasswd 2>/dev/null || \
  usermod -p "$(openssl passwd -1 "$PW")" rstudio

# RStudio needs rsession-which-r AND the database.conf pointing at R
mkdir -p /etc/rstudio
cat > /etc/rstudio/rserver.conf <<'RSCFG'
rsession-which-r=/opt/R/4.4.2/bin/R
RSCFG
# Ensure RStudio's session process can find R
cat > /etc/rstudio/database.conf <<'DBCFG'
provider=sqlite
directory=/var/lib/rstudio-server
DBCFG
# Default to the project directory so data and scripts are immediately visible
cat > /etc/rstudio/rsession.conf <<'RSESSCFG'
session-default-working-dir=/home/rstudio/kefir4all
session-default-new-project-dir=/home/rstudio/kefir4all
RSESSCFG

mkdir -p /var/run/rstudio-server /var/log/rstudio /var/lib/rstudio-server

# Fix permissions on mounted project directory so RStudio user can read it.
# Host-mounted files may be owned by arbitrary UIDs.
if [ -d /home/rstudio/kefir4all ]; then
  chmod -R a+rX /home/rstudio/kefir4all 2>/dev/null || true
fi

/opt/R/4.4.2/bin/Rscript -e 'cat("R ready:", R.version.string, "| packages:", length(rownames(installed.packages())), "\n")'

exec /usr/lib/rstudio-server/bin/rserver --server-daemonize=0
