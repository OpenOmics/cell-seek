FROM rocker/r-ver:4.5.0

LABEL maintainer="OpenOmics"
LABEL description="Docker image for Signac scATAC-seq QC analysis"

# Install system dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    zlib1g-dev \
    ca-certificates \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Test internet connectivity before installing R packages
RUN curl -s https://cloud.r-project.org/src/contrib/PACKAGES | head -5 || \
    (echo "ERROR: Cannot access CRAN repository" && exit 1)

# Configure R to use CRAN mirror and install BiocManager
RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages('BiocManager', dependencies=TRUE)"

# Install Bioconductor packages
# RUN R -e "BiocManager::install(c( \
#     'Signac', \
#     'Seurat', \
#     'GenomicRanges', \
#     'rtracklayer', \
#     'SeuratObject', \
#     'BSgenome', \
#     'GenomeInfoDb', \
#     'IRanges', \
#     'S4Vectors', \
#     'BiocGenerics' \
#     ), ask=FALSE, update=TRUE)"

# # Install additional CRAN packages
# RUN R -e "install.packages(c( \
#     'ggplot2', \
#     'gridExtra', \
#     'dplyr', \
#     'stringr', \
#     'optparse', \
#     'cluster', \
#     'devtools', \
#     'hdf5r', \
#     'knitr', \
#     'rmarkdown', \
#     'flextable', \
#     'base64enc', \
#     'xfun', \
#     'tools' \
#     ), repos='https://cloud.r-project.org')"

# # Create scripts directory
# RUN mkdir -p /opt/scripts
# RUN export PATH="/opt/scripts:${PATH}"

# # Set working directory
# WORKDIR /data

# # Default command
# CMD ["/bin/bash"]
