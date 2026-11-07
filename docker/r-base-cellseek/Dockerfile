FROM rocker/r-ver:4.4.0

LABEL maintainer="OpenOmics"
LABEL description="Docker image for cell-seek's Seurat QC/integration reports (gex, cite, multi pipelines)"

# rocker/r-ver is Ubuntu 22.04 (jammy) based and ships with CRAN preconfigured
# to a dated Posit Package Manager snapshot matching R 4.4.0, so packages
# install as precompiled binaries.
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libgit2-dev \
    libglpk-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff-dev \
    libjpeg-dev \
    libgsl-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libcairo2-dev \
    libxt-dev \
    libwebp-dev \
    libnlopt-dev \
    libudunits2-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install pak (r-ver's CRAN repo is a PPM binary snapshot, so plain
# install.packages works; pak is then used for dependency resolution)
RUN R -e 'install.packages("pak")'

# Install CRAN + Bioconductor packages. pak adds the Bioconductor repos that
# match R 4.4 (3.19) automatically. r-ver's default CRAN snapshot (spring 2024)
# is too old: current Bioc 3.19 packages need matrixStats >= 1.4.1 (Sep 2024).
# Use a later jammy binary snapshot that still supports R 4.4.
ENV CRAN_SNAPSHOT="https://packagemanager.posit.co/cran/__linux__/jammy/2025-03-01"
RUN R -e ' \
  options(repos = c(CRAN = Sys.getenv("CRAN_SNAPSHOT"))); \
  options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
  pak::pkg_install(c( \
    "BiocManager", "devtools", "remotes", "optparse", "dplyr", "stringr", "ggplot2", \
    "gridExtra", "cluster", "base64enc", "xfun", "knitr", "rmarkdown", "rmdformats", \
    "kableExtra", "flextable", "DT", "pander", "SoupX", "Seurat", "harmony", \
    "SingleCellExperiment", "SummarizedExperiment", "scater", "scRNAseq", "SingleR", "DropletUtils" \
  ), ask = FALSE)'

# Fail the build early if any key package did not install
RUN R -e 'for (p in c("Seurat","harmony","SoupX","SingleR","DropletUtils","scater","scRNAseq","flextable","kableExtra","rmdformats")) library(p, character.only = TRUE)'

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
