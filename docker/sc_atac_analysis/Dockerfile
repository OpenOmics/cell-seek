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
RUN R -e "install.packages(c('hdf5r', 'remotes', 'sparseMatrixStats'), repos='https://cloud.r-project.org')"

# Install Signac and Seurat from github
RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5')"
RUN R -e "setRepositories(ind=1:3); install.packages('Signac')"

# Create scripts directory
RUN mkdir -p /opt/scripts
RUN export PATH="/opt/scripts:${PATH}"

ADD signacSampleQC.R signacSampleQCReport.Rmd signacMultiSampleQC.R signacMultiSampleQCreport.Rmd /opt/scripts/

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
