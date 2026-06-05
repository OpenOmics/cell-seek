#!/usr/bin/env Rscript
# Load a Visium HD dataset into Seurat, run QC, normalization, clustering, UMAP,
# and save the resulting Seurat object as an RDS file.

library(Seurat)
library(hdf5r)
library(ggplot2)
library(future)
library(argparse)

# give good tracebacks on errors
options(error = function() traceback(2))
err <- function(...) {
  cat(sprintf(...), sep = "\n", file = stderr())
}
fatal <- function(...) {
  err(...)
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
parser <- ArgumentParser(
  description = "Ingest a 10x Visium HD dataset into Seurat and produce an RDS file."
)
parser$add_argument(
  "-d", "--data-dir",
  dest    = "data_dir",
  help    = "Path to the Visium HD output directory (contains binned_outputs/ and spatial/).",
  required = TRUE
)
parser$add_argument(
  "-o", "--outdir",
  dest    = "outdir",
  help    = "Output directory for the RDS file and PDF plots (created if absent).",
  default = "."
)
parser$add_argument(
  "-b", "--bin-size",
  dest    = "bin_size",
  type    = "integer",
  help    = "Bin size in µm to load (2, 8, or 16). Default: 16.",
  default = 16
)
parser$add_argument(
  "-s", "--slice",
  dest    = "slice",
  help    = "Slice name passed to Load10X_Spatial(). Default: basename of --data-dir.",
  default = NULL
)
parser$add_argument(
  "--assay",
  dest    = "assay",
  help    = "Assay name passed to Load10X_Spatial(). Default: Spatial.",
  default = "Spatial"
)
parser$add_argument(
  "--min-features",
  dest    = "min_features",
  type    = "integer",
  help    = "Minimum number of features (genes) per bin to retain. Default: 200.",
  default = 200
)
parser$add_argument(
  "--max-features",
  dest    = "max_features",
  type    = "integer",
  help    = "Maximum number of features (genes) per bin to retain. Default: 8000.",
  default = 8000
)
parser$add_argument(
  "--max-pct-mt",
  dest    = "max_pct_mt",
  type    = "double",
  help    = "Maximum mitochondrial read percentage per bin to retain. Default: 20.",
  default = 20
)
parser$add_argument(
  "--nfeatures-hvg",
  dest    = "nfeatures_hvg",
  type    = "integer",
  help    = "Number of highly variable genes for FindVariableFeatures(). Default: 3000.",
  default = 3000
)
parser$add_argument(
  "--npcs",
  dest    = "npcs",
  type    = "integer",
  help    = "Number of principal components for PCA / UMAP / clustering. Default: 30.",
  default = 30
)
parser$add_argument(
  "--resolution",
  dest    = "resolution",
  type    = "double",
  help    = "Clustering resolution for FindClusters(). Default: 0.5.",
  default = 0.5
)
parser$add_argument(
  "--output-rds",
  dest    = "output_rds",
  help    = paste0("Filename (not path) of the output RDS file. ",
                   "Default: visium_hd_<slice>_<bin_size>um.rds"),
  default = NULL
)
parser$add_argument(
  "--memory-gb",
  dest    = "memory_gb",
  type    = "double",
  help    = "future globals memory cap in GB. Default: 32.",
  default = 32
)

args <- parser$parse_args()

# ---------------------------------------------------------------------------
# Resolve derived values
# ---------------------------------------------------------------------------
options(future.globals.maxSize = args$memory_gb * 1024^3)

data_dir    <- normalizePath(args$data_dir, mustWork = FALSE)
outdir      <- normalizePath(args$outdir,   mustWork = FALSE)
bin_size    <- args$bin_size
slice_name  <- if (!is.null(args$slice)) args$slice else basename(data_dir)

if (is.null(args$output_rds)) {
  output_rds <- file.path(outdir, paste0("visium_hd_", slice_name, "_", bin_size, "um.rds"))
} else {
  output_rds <- file.path(outdir, args$output_rds)
}

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if (!dir.exists(data_dir)) {
  fatal(paste0("Data directory not found: ", data_dir))
}
binned_dir <- file.path(data_dir, "binned_outputs")
spatial_dir <- file.path(data_dir, "spatial")
if (!dir.exists(binned_dir)) {
  fatal(paste0("binned_outputs/ not found inside data directory: ", binned_dir))
}
if (!dir.exists(spatial_dir)) {
  fatal(paste0("spatial/ not found inside data directory: ", spatial_dir))
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat(paste0("INFO - Created output directory: ", outdir, "\n"))
}

# List available bin sizes
available_bins <- list.dirs(binned_dir, full.names = FALSE, recursive = FALSE)
cat("INFO - Available bin sizes:", paste(available_bins, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat(paste0("INFO - Loading Visium HD data: ", data_dir,
           " | bin size: ", bin_size, "µm | slice: ", slice_name, "\n"))

visium_hd <- Load10X_Spatial(
  data.dir      = data_dir,
  bin.size      = bin_size,
  assay         = args$assay,
  slice         = slice_name,
  filter.matrix = TRUE
)

cat(paste0("INFO - Seurat object loaded: ", ncol(visium_hd),
           " bins, ", nrow(visium_hd), " features\n"))
cat(paste0("INFO - Metadata columns: ",
           paste(colnames(visium_hd@meta.data), collapse = ", "), "\n"))

default_assay <- DefaultAssay(visium_hd)
cat(paste0("INFO - Default assay: ", default_assay, "\n"))

ncount_col  <- paste0("nCount_",   default_assay)
nfeature_col <- paste0("nFeature_", default_assay)

# ---------------------------------------------------------------------------
# QC metrics
# ---------------------------------------------------------------------------
cat("INFO - Calculating QC metrics...\n")
visium_hd[["percent.mt"]]   <- PercentageFeatureSet(visium_hd, pattern = "^MT-")
visium_hd[["percent.ribo"]] <- PercentageFeatureSet(visium_hd, pattern = "^RP[SL]")

pdf(file.path(outdir, "qc_metrics.pdf"), width = 12, height = 4)
VlnPlot(visium_hd,
        features = c(ncount_col, nfeature_col, "percent.mt"),
        pt.size  = 0.1,
        ncol     = 3)
dev.off()

pdf(file.path(outdir, "spatial_qc.pdf"), width = 12, height = 8)
print(SpatialFeaturePlot(visium_hd, features = ncount_col)   + ggtitle("UMI Counts per Bin"))
print(SpatialFeaturePlot(visium_hd, features = nfeature_col) + ggtitle("Genes per Bin"))
print(SpatialFeaturePlot(visium_hd, features = "percent.mt") + ggtitle("Mitochondrial %"))
dev.off()

# ---------------------------------------------------------------------------
# Filter low-quality bins
# ---------------------------------------------------------------------------
cat(paste0("INFO - Filtering bins: features [", args$min_features, ", ",
           args$max_features, "], percent.mt < ", args$max_pct_mt, "\n"))

metadata  <- visium_hd@meta.data
keep_bins <- metadata[[nfeature_col]] > args$min_features &
             metadata[[nfeature_col]] < args$max_features &
             metadata$percent.mt      < args$max_pct_mt
keep_bins[is.na(keep_bins)] <- FALSE

cat(paste0("INFO - Bins passing QC: ", sum(keep_bins),
           " of ", ncol(visium_hd), "\n"))

if (sum(keep_bins) == 0) {
  fatal("Pre-flight ERROR: no bins passed QC filters. Relax --min-features, --max-features, or --max-pct-mt.")
}

visium_hd <- visium_hd[, keep_bins]
cat(paste0("INFO - Bins retained after filtering: ", ncol(visium_hd), "\n"))

# ---------------------------------------------------------------------------
# Normalization, variable features, scaling
# ---------------------------------------------------------------------------
cat("INFO - Running LogNormalize...\n")
visium_hd <- NormalizeData(visium_hd,
                           normalization.method = "LogNormalize",
                           scale.factor         = 10000)
visium_hd <- FindVariableFeatures(visium_hd,
                                  selection.method = "vst",
                                  nfeatures        = args$nfeatures_hvg)
visium_hd <- ScaleData(visium_hd, features = VariableFeatures(visium_hd))

# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------
cat(paste0("INFO - Running PCA (", args$npcs, " PCs)...\n"))
visium_hd <- RunPCA(visium_hd,
                    features = VariableFeatures(visium_hd),
                    verbose  = FALSE)

pdf(file.path(outdir, "elbow_plot.pdf"), width = 6, height = 4)
ElbowPlot(visium_hd, ndims = 50)
dev.off()

# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------
cat(paste0("INFO - Finding neighbors and clusters (resolution = ",
           args$resolution, ")...\n"))
visium_hd <- FindNeighbors(visium_hd, reduction = "pca", dims = 1:args$npcs)
visium_hd <- FindClusters(visium_hd,  resolution = args$resolution, verbose = FALSE)

# ---------------------------------------------------------------------------
# UMAP
# ---------------------------------------------------------------------------
cat("INFO - Running UMAP...\n")
visium_hd <- RunUMAP(visium_hd, reduction = "pca", dims = 1:args$npcs, verbose = FALSE)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
pdf(file.path(outdir, "spatial_clusters.pdf"), width = 10, height = 8)
print(SpatialDimPlot(visium_hd, label = TRUE, label.size = 3) +
      ggtitle(paste0("Spatial Clusters - ", slice_name)))
dev.off()

pdf(file.path(outdir, "umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(visium_hd, reduction = "umap", label = TRUE))
dev.off()

# ---------------------------------------------------------------------------
# Save RDS
# ---------------------------------------------------------------------------
cat(paste0("INFO - Saving Seurat object to: ", output_rds, "\n"))
saveRDS(visium_hd, file = output_rds)

cat("\nINFO - Analysis complete!\n")
cat(paste0("INFO - RDS saved to:  ", output_rds, "\n"))
cat("INFO - Plots written:\n")
for (f in c("qc_metrics.pdf", "spatial_qc.pdf", "elbow_plot.pdf",
            "spatial_clusters.pdf", "umap_clusters.pdf")) {
  cat(paste0("  - ", file.path(outdir, f), "\n"))
}

cat("\nSession Info:\n")
sessionInfo()

