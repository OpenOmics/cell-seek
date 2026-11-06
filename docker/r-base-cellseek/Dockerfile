FROM r-base:4.4.0

LABEL maintainer="OpenOmics"
LABEL description="Docker image for cell-seek's Seurat QC/integration reports (gex, cite, multi pipelines)"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libhdf5-serial-dev \
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
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install pak, the fast R package manager
RUN R -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))'

# Install CRAN + Bioconductor packages from pre-compiled binaries served by
# Posit Package Manager (PPM), targeting this image's OS (Debian trixie, the
# codename of "debian:testing" that r-base:4.4.0 is built on) and the
# Bioconductor release matched to R 4.4 (3.19). All packages below were
# confirmed to have PPM binaries for this OS/R/Bioc combination.
RUN R -e ' \
  options(repos = c( \
    CRAN = "https://packagemanager.posit.co/cran/__linux__/trixie/latest", \
    BIOC = "https://packagemanager.posit.co/bioconductor/__linux__/trixie/packages/3.19/bioc", \
    BIOCDATAEXP = "https://packagemanager.posit.co/bioconductor/__linux__/trixie/packages/3.19/data/experiment", \
    BIOCDATAANN = "https://packagemanager.posit.co/bioconductor/__linux__/trixie/packages/3.19/data/annotation" \
  )); \
  pak::pkg_install(c( \
    "BiocManager", "devtools", "remotes", "optparse", "dplyr", "stringr", "ggplot2", \
    "gridExtra", "cluster", "base64enc", "xfun", "knitr", "rmarkdown", "rmdformats", \
    "kableExtra", "flextable", "DT", "pander", "SoupX", "Seurat", \
    "SingleCellExperiment", "SummarizedExperiment", "scater", "scRNAseq", "SingleR", "DropletUtils" \
  ), ask = FALSE)'

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
