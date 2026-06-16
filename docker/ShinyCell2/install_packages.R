#!/usr/bin/env Rscript

# Number of CPUs processes
# to use for parallelizing
# the install of N packages
use_ncpus <- max(parallel::detectCores() - 2, 2)

# CRAN packages,
# add any missing/required
# CRAN packages to the list
# directly below. This script
# will install them if they
# are not already installed.
# Install via:
#   install.packages('packageName', repos='http://cran.r-project.org')
cran_packages <- c(
  "optparse",
  "data.table",
  "Matrix",
  "hdf5r",
  "reticulate",
  "R.utils",
  "ggplot2",
  "gridExtra",
  "glue",
  "readr",
  "future",
  "RColorBrewer",
  "devtools",
  "Seurat",
  "remotes",
  "DT",
  "BiocManager",
  "shinyhelper",
  "argparse",
  "ggdendro",
  "hdf5r",
  "magrittr",
  "shiny",
  "shinyhelper",
  "ggpubr",
  "ggrepel"
)


# Install missing CRAN packages
install.packages(
  setdiff(
    cran_packages, rownames(installed.packages())
  ),
  Ncpus = use_ncpus,
  repos = "http://cran.r-project.org"
)

# Bioconductor packages,
# add any missing/required
# Bioconductor packages to the list
# directly below. This script will
# install them if they are missing.
# Install via:
#   BiocManager::install('packageName')
bioc_packages <- c(
  # limma
)

# Install missing Bioconductor packages
for (p in bioc_packages) {
  if (!p %in% rownames(installed.packages())) {
    # Package missing
    print(paste("Installing missing Bioconductor package...", p))
    BiocManager::install(p, ask = FALSE, update = FALSE, Ncpus = use_ncpus)
  }
}

# not on CRAN
remotes::install_github("satijalab/seurat-wrappers")
