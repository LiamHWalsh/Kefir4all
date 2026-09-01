# Kefir4All reproducible computation image
# ------------------------------------------
# Builds a self-contained R 4.4.2 + Bioconductor-style RStudio Server
# environment with the EXACT 551-package library baked in, so users never
# install or compile anything.
#
#   docker build -t kefir4all:4.4.2 .
#   docker run --rm -p 8787:8787 -e PASSWORD=kefir4all kefir4all:4.4.2
#   # open http://localhost:8787 (login: user / your PASSWORD)
#
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# System libraries matching this Codespace's build environment.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl gcc g++ gfortran make pkg-config \
      libssl-dev libcurl4-openssl-dev libxml2-dev zlib1g-dev \
      libpng-dev libjpeg-dev libtiff-dev libfontconfig1-dev libfreetype6-dev \
      libharfbuzz-dev libfribidi-dev libgit2-dev libicu-dev libbz2-dev \
      liblzma-dev libpcre2-dev libblas-dev liblapack-dev \
      libuv1-dev libproj-dev proj-bin libmagick++-dev libcairo2-dev \
      libglpk-dev \
      locales sudo \
 && rm -rf /var/lib/apt/lists/* \
 && locale-gen en_US.UTF-8 && update-locale
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# R 4.4.2 (prebuilt Posit binary; installs to /opt/R/4.4.2)
COPY r442.deb /tmp/r442.deb
RUN dpkg -i /tmp/r442.deb || (apt-get update && apt-get install -f -y && dpkg -i /tmp/r442.deb); \
    rm -f /tmp/r442.deb && /opt/R/4.4.2/bin/R --version

# Persist the frozen repository snapshot for reproducibility.
RUN echo 'options(repos=c(POSIT="https://packagemanager.posit.co/cran/__linux__/noble/2025-06-15",CRAN="https://cloud.r-project.org"))' \
    >> /opt/R/4.4.2/lib/R/etc/Rprofile.site

# The exact working 551-package library (tar'd from the validated env).
# Tar entries look like ./4.4/<pkg>, so strip two components to land
# packages directly in R's library dir.
COPY R-library-v2.tar.gz /tmp/R-library.tar.gz
RUN mkdir -p /opt/R/4.4.2/lib/R/library \
    && tar -xzf /tmp/R-library.tar.gz -C /opt/R/4.4.2/lib/R/library --strip-components=2 \
    && rm -f /tmp/R-library.tar.gz \
    && ls /opt/R/4.4.2/lib/R/library | wc -l
# ---- Slim unused packages to reduce image size (~733 MB) ----
# Computes strong-dependency closure for the 75 analysis packages;
# removes the remaining ~126 packages. Reproducible.
COPY packages_required.txt /tmp/packages_required.txt
COPY slim_packages.R /tmp/slim_packages.R
RUN /opt/R/4.4.2/bin/Rscript /tmp/slim_packages.R /tmp/packages_required.txt \
 && rm /tmp/packages_required.txt /tmp/slim_packages.R

# RStudio Server (prebuilt, jammy build is binary-compatible with noble).
COPY rstudio.deb /tmp/rstudio.deb
RUN dpkg -i /tmp/rstudio.deb || (apt-get update && apt-get install -f -y && dpkg -i /tmp/rstudio.deb); \
    rm -f /tmp/rstudio.deb

# Point RStudio at the non-standard R location.
RUN mkdir -p /etc/rstudio \
    && echo "rsession-which-r=/opt/R/4.4.2/bin/R" > /etc/rstudio/rserver.conf

# Convenience user for RStudio login. The RStudio deb already creates
# 'rstudio'; make this step tolerant so it never fails a rebuild.
# Password is (re)set at container start via the PASSWORD env (start.sh).
RUN bash -lc '\
  id rstudio >/dev/null 2>&1 || useradd -m -s /bin/bash rstudio ; \
  adduser rstudio sudo >/dev/null 2>&1 || true ; \
  chsh -s /bin/bash rstudio >/dev/null 2>&1 || true ; \
  echo "rstudio:rstudio" | chpasswd || true ; \
  true'

COPY start.sh /usr/local/bin/start.sh
RUN chmod +x /usr/local/bin/start.sh

EXPOSE 8787
CMD ["/usr/local/bin/start.sh"]