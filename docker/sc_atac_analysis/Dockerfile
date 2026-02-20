FROM rocker/r-ver:4.5.0

LABEL maintainer="OpenOmics"
LABEL description="Docker image for Signac scATAC-seq QC analysis"

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
    pandoc \
    zlib1g-dev \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN R -e "install.packages(c('hdf5r', 'optparse', 'rtracklayer', 'dplyr', 'ggplot2', 'cluster', 'stringr', 'Seurat', 'BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('GenomeInfoDb', 'GenomicRanges', 'IRanges', 'Rsamtools', 'S4Vectors', 'BiocGenerics', 'SeuratObject'))"
RUN R -e "install.packages(c('remotes', 'XML', 'hdf5r', 'restfulr', 'devtools')); BiocManager::install(c('Rsamtools', 'rtracklayer', 'fastmatch', 'RcppRoll')); remotes::install_github('stuart-lab/signac', dependencies = FALSE)"

# Create scripts directory
RUN mkdir -p /opt/scripts
RUN export PATH="/opt/scripts:${PATH}"

ADD signacSampleQC.R signacSampleQCReport.Rmd signacMultiSampleQC.R signacMultiSampleQCreport.Rmd /opt/scripts/

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
