FROM rocker/r-ver:4.3.2

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
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN R -e "install.packages(c('pak', 'hdf5r'), repos='https://cloud.r-project.org')"

# Install Bioconductor packages
RUN R -e "pak::pkg_install(c('Signac', 'Seurat', 'GenomicRanges', \
                             'rtracklayer', 'SeuratObject', 'BSgenome', \
                             'GenomeInfoDb', 'IRanges', 'S4Vectors', \
                             'BiocGenerics', 'hdf5r',     'ggplot2', \
                             'gridExtra', 'dplyr', 'stringr', \
                             'optparse', 'cluster', 'devtools', \
                             'knitr', 'rmarkdown',  'flextable',  \
                             'base64enc', 'xfun', 'tools'))"

# Create scripts directory
RUN mkdir -p /opt/scripts
RUN export PATH="/opt/scripts:${PATH}"

ADD signacSampleQC.R signacSampleQCReport.Rmd signacMultiSampleQC.R signacMultiSampleQCreport.Rmd /opt/scripts/

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
