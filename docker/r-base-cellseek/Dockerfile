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

# Install R packages from CRAN
RUN R -e "install.packages(c('BiocManager', 'devtools', 'remotes', 'optparse', 'dplyr', 'stringr', 'ggplot2', 'gridExtra', 'cluster', 'base64enc', 'xfun', 'knitr', 'rmarkdown', 'rmdformats', 'kableExtra', 'flextable', 'DT', 'pander', 'SoupX', 'Seurat'), repos='https://cloud.r-project.org')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('SingleCellExperiment', 'SummarizedExperiment', 'scater', 'scRNAseq', 'SingleR', 'DropletUtils'))"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
