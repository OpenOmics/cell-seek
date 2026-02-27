#!/usr/bin/env Rscript

# Enable better error tracebacks
options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    cat("Traceback:\n")
    traceback(2)
  }
  if (!interactive()) {
    quit(status = 1)
  }
})

library(Signac)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(optparse)
library(dplyr)
library(stringr)
library(cluster)

option_list <- list(
  make_option(
    "--matrix",
    type = "character",
    default = NA,
    help = "Comma-separated paths to the matrix file(s) (one per sample)"
  ),
  make_option(
    "--fragments",
    type = "character",
    action = "store",
    default = NA,
    help = "Comma-separated paths to the fragments file(s) (one per sample)"
  ),
  make_option(
    "--barcodes",
    type = "character",
    action = "store",
    default = NA,
    help = "Comma-separated paths to the barcodes file(s) (one per sample)"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    action = "store",
    default = NA,
    help = "Path to the output directory"
  ),
  make_option(
    c("-p", "--project"),
    type = "character",
    action = "store",
    default = "project",
    help = "Project name to be stored within Seurat metadata"
  ),
  make_option(
    c("-f", "--filterfile"),
    type = "character",
    action = "store",
    default = NA,
    help = "CSV file containing filters to be applied (one row per sample)"
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    action = "store",
    default = NA,
    help = "Metadata file with information to add to samples"
  ),
  make_option(
    "--genes",
    type = "character",
    action = "store",
    default = NA,
    help = "Path to GTF annotation file (e.g., gencode.v32.primary_assembly.annotation.gtf.gz)"
  ),
  make_option(
    c("-g", "--genome"),
    type = "character",
    action = "store",
    default = "hg38",
    help = "Genome name (e.g., hg38, mm10) - used for seqlevels style"
  ),
  make_option(
    c("-n", "--minfragments"),
    type = "integer",
    action = "store",
    default = 500,
    help = "Minimum number of fragments for initial cell filtering (default: 500)"
  ),
  make_option(
    c("-s", "--samples"),
    type = "character",
    action = "store",
    help = "Comma-separated list of sample names"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Expand output path to full absolute path
if (!is.na(opt$output)) {
  opt$output <- normalizePath(opt$output, mustWork = FALSE)
}

## ----Parse Sample Information----
if (grepl(",", opt$samples)) {
  samples <- strsplit(opt$samples, ",")[[1]] %>% trimws()
} else {
  samples <- trimws(opt$samples)
}

## ----Parse Comma-Separated Input Paths----
parse_csv_arg <- function(arg) {
  if (grepl(",", arg)) {
    strsplit(arg, ",")[[1]] %>% trimws()
  } else {
    trimws(arg)
  }
}

matrix_paths   <- parse_csv_arg(opt$matrix)
fragment_paths <- parse_csv_arg(opt$fragments)
barcode_paths  <- parse_csv_arg(opt$barcodes)

# Validate that the number of paths matches the number of samples
if (length(matrix_paths) != length(samples)) {
  stop(paste0(
    "Number of --matrix paths (", length(matrix_paths),
    ") does not match the number of --samples (", length(samples), ")"
  ))
}
if (length(fragment_paths) != length(samples)) {
  stop(paste0(
    "Number of --fragments paths (", length(fragment_paths),
    ") does not match the number of --samples (", length(samples), ")"
  ))
}
if (length(barcode_paths) != length(samples)) {
  stop(paste0(
    "Number of --barcodes paths (", length(barcode_paths),
    ") does not match the number of --samples (", length(samples), ")"
  ))
}

## ----Load Reference Annotations----
if (is.na(opt$genes)) {
  stop("GTF annotation file path is required (use -r/--reference)")
}

cat(paste0("Loading annotations from ", opt$genes, "...\n"))
annotations <- import(opt$genes)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- opt$genome

# Rename gene_type to gene_biotype if needed (GENCODE uses gene_type)
if (
  "gene_type" %in%
    colnames(mcols(annotations)) &
    !"gene_biotype" %in% colnames(mcols(annotations))
) {
  annotations$gene_biotype <- annotations$gene_type
}

if (dir.exists(opt$output)) {
  unlink(opt$output, recursive = TRUE)
} else {
  cat(paste0("Creating output directory ", opt$output, "...\n"))
}
dir.create(opt$output, recursive = TRUE)

figures <- list()

## ----Load Individual Samples and Create Union Peak Set----
cat("Loading samples and creating union peak set...\n")

sample_objects <- list()
all_peaks <- GRangesList()

for (i in 1:length(samples)) {
  cat(paste0("Loading sample ", samples[i], "...\n"))

  counts <- Read10X_h5(matrix_paths[i])
  cellranger_metrics <- read.csv(
    file = barcode_paths[i],
    header = TRUE,
    row.names = 1
  )

  # Filter low-depth cells before creating ChromatinAssay
  cells_before <- nrow(cellranger_metrics)
  cellranger_metrics <- cellranger_metrics[
    cellranger_metrics$passed_filters > opt$minfragments,
  ]
  cells_after <- nrow(cellranger_metrics)
  cat(paste0(
    "  Filtered ",
    cells_before - cells_after,
    " cells with <",
    opt$minfragments,
    " fragments\n"
  ))

  # Subset counts matrix to cells that passed filter AND exist in counts matrix
  cells_to_keep <- intersect(colnames(counts), rownames(cellranger_metrics))
  counts <- counts[, cells_to_keep]
  cellranger_metrics <- cellranger_metrics[cells_to_keep, ]
  cat(paste0("  Retained ", length(cells_to_keep), " cells for analysis\n"))

  # Get peaks for this sample
  peaks <- StringToGRanges(rownames(counts), sep = c(":", "-"))
  all_peaks[[samples[i]]] <- peaks

  # Create ChromatinAssay with sample-specific peaks for now
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = genome(annotations)[1],
    fragments = fragment_paths[i],
    min.cells = 1,
    min.features = 1
  )

  seur <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = opt$project,
    meta.data = cellranger_metrics
  )

  seur$Sample <- samples[i]

  # Add custom metadata
  if (!is.na(opt$metadata)) {
    sample_metadata <- read.csv(opt$metadata)
    index <- which(
      sample_metadata[, grep(
        "sample",
        colnames(sample_metadata),
        ignore.case = T
      )] ==
        samples[i]
    )
    if (length(index) > 0) {
      for (header in grep(
        "sample",
        colnames(sample_metadata),
        ignore.case = T,
        invert = T,
        value = T
      )) {
        seur[[header]] <- sample_metadata[index, header]
      }
    }
  }

  sample_objects[[samples[i]]] <- seur
}

# Create union peak set
cat("Creating union peak set...\n")
union_peaks <- reduce(unlist(all_peaks))

# Filter out peaks that are too wide or too narrow
peakwidths <- width(union_peaks)
union_peaks <- union_peaks[peakwidths < 10000 & peakwidths > 20]
cat(paste0(
  "Union peak set contains ",
  length(union_peaks),
  " peaks after filtering\n"
))

# Requantify each sample with union peaks
cat("Requantifying samples with union peak set...\n")
for (i in 1:length(samples)) {
  cat(paste0("Requantifying sample ", samples[i], "...\n"))

  # Get fragment path
  fragment_path <- fragment_paths[i]

  # Create Fragment object
  frags <- CreateFragmentObject(
    path = fragment_path,
    cells = colnames(sample_objects[[samples[i]]])
  )

  # Count fragments in union peaks
  counts <- FeatureMatrix(
    fragments = frags,
    features = union_peaks,
    cells = colnames(sample_objects[[samples[i]]])
  )

  # Create new ChromatinAssay with union peaks
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    fragments = fragment_path,
    genome = genome(annotations)[1],
    min.cells = 1,
    min.features = 1
  )

  # Preserve metadata
  metadata <- sample_objects[[samples[i]]]@meta.data

  # Create new Seurat object
  sample_objects[[samples[i]]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = opt$project,
    meta.data = metadata
  )
}

## ----Merge Samples----
cat("Merging samples...\n")
if (length(sample_objects) > 1) {
  seur <- merge(
    x = sample_objects[[1]],
    y = sample_objects[2:length(sample_objects)],
    add.cell.ids = samples
  )
} else {
  seur <- sample_objects[[1]]
}

## ----Add Annotations----
Annotation(seur) <- annotations

## ----Create Genome-Wide Region for FragmentHistogram----
# GTF files don't contain chromosome lengths, so we'll get them from the Seurat object
# The ChromatinAssay should have genome information from when it was created
chrom_assay <- seur[["peaks"]]
genome_info <- seqinfo(chrom_assay)
cat("Genome info from ChromatinAssay:\n")
print(genome_info)

# Extract chromosome names and lengths
chr_names <- seqlevels(genome_info)
chr_lengths <- seqlengths(genome_info)

# Check if we have valid seqlengths
if (all(is.na(chr_lengths))) {
  stop("ERROR: No chromosome lengths available in the genome information. ",
       "Please ensure the genome parameter is set correctly when creating the ChromatinAssay, ",
       "or provide a genome annotation with seqlength information.")
}

# Remove any chromosomes with NA lengths
valid_chrs <- !is.na(chr_lengths)
chr_names <- chr_names[valid_chrs]
chr_lengths <- chr_lengths[valid_chrs]

if (length(chr_names) == 0) {
  stop("ERROR: No valid chromosomes with length information found.")
}

# Create region strings in the format "chr-start-end" for each chromosome
genome_region <- paste0(chr_names, "-1-", chr_lengths)
cat("Sample genome regions:\n")
cat(paste(head(genome_region, 3), collapse = "\n"), "\n")

## ----Calculate QC Metrics----
cat("Calculating QC metrics...\n")
seur <- NucleosomeSignal(object = seur)
seur <- TSSEnrichment(object = seur, fast = FALSE)

# Calculate blacklist ratio using Signac's FractionCountsInRegion
# This computes the fraction directly from fragments overlapping blacklist regions,
# rather than relying on Cell Ranger singlecell.csv columns (blacklist_region_fragments,
# peak_region_fragments) which may not exist after requantification with union peaks.
if (opt$genome %in% c("hg38", "hg2024")) {
  seur$blacklist_ratio <- FractionCountsInRegion(
    object = seur,
    assay = "peaks",
    regions = blacklist_hg38_unified
  )
} else if (opt$genome %in% c("mm10", "mm2024")) {
  seur$blacklist_ratio <- FractionCountsInRegion(
    object = seur,
    assay = "peaks",
    regions = blacklist_mm10
  )
} else {
  warning("Genome not recognized for blacklist ratio calculation. Setting blacklist_ratio to NA.")
  seur$blacklist_ratio <- NA
}

# Calculate pct_reads_in_peaks from Cell Ranger singlecell.csv metadata
# These columns are preserved in metadata from the initial load
if (
  "peak_region_fragments" %in% colnames(seur@meta.data) &&
    "passed_filters" %in% colnames(seur@meta.data)
) {
  seur$pct_reads_in_peaks <- seur$peak_region_fragments /
    seur$passed_filters *
    100
  cat("Calculated pct_reads_in_peaks from Cell Ranger metadata.\n")
} else {
  # Fallback: compute fraction of counts in peaks from the assay
  warning(
    "peak_region_fragments and/or passed_filters columns not found in metadata. ",
    "Computing pct_reads_in_peaks as fraction of counts in peak regions."
  )
  total_counts <- Matrix::colSums(GetAssayData(seur, assay = "peaks", slot = "counts"))
  seur$pct_reads_in_peaks <- (total_counts / seur$nCount_peaks) * 100
}


## ----Pre-Filter Feature Plot----
for (sample in samples) {
  seur_subset <- subset(seur, subset = Sample == sample)
  # Only use group.by if there are multiple samples
  if (length(samples) > 1) {
    plot1 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nCount_peaks",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot2 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nCount_peaks",
      feature2 = "nFeature_peaks"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot3 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nucleosome_signal",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
  } else {
    plot1 <- FeatureScatter(
      seur_subset,
      feature1 = "nCount_peaks",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot2 <- FeatureScatter(
      seur_subset,
      feature1 = "nCount_peaks",
      feature2 = "nFeature_peaks"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot3 <- FeatureScatter(
      seur_subset,
      feature1 = "nucleosome_signal",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
  }

  png(
    file.path(opt$output, paste0("PreFilter_Feature_Plot_", sample, ".png")),
    height = 5,
    width = 15,
    units = "in",
    res = 300
  )
  print(plot1 + plot2 + plot3)
  dev.off()

  figures[[paste0("PreFilter_Feature_Plot_", sample)]] <- plot1 + plot2
}

## ----Cell Quality Thresholds - Default----
defaultThreshold <- function(seur, sample_name) {
  thresh <- list()
  sample_cells <- colnames(seur)[seur$Sample == sample_name]

  thresh["nFeature_peaks_low"] <- expm1(
    median(log1p(seur$nFeature_peaks[sample_cells])) -
      3 * mad(log1p(seur$nFeature_peaks[sample_cells]))
  ) %>%
    round
  thresh["nFeature_peaks_high"] <- expm1(
    median(log1p(seur$nFeature_peaks[sample_cells])) +
      3 * mad(log1p(seur$nFeature_peaks[sample_cells]))
  ) %>%
    round
  thresh["nCount_peaks_low"] <- expm1(
    median(log1p(seur$nCount_peaks[sample_cells])) -
      3 * mad(log1p(seur$nCount_peaks[sample_cells]))
  ) %>%
    round
  thresh["nCount_peaks_high"] <- expm1(
    median(log1p(seur$nCount_peaks[sample_cells])) +
      3 * mad(log1p(seur$nCount_peaks[sample_cells]))
  ) %>%
    round
  thresh["TSS.enrichment_low"] <- max(
    expm1(
      median(log1p(seur$TSS.enrichment[sample_cells])) -
        3 * mad(log1p(seur$TSS.enrichment[sample_cells]))
    ),
    1
  )
  thresh["nucleosome_signal_high"] <- expm1(
    median(log1p(seur$nucleosome_signal[sample_cells])) +
      3 * mad(log1p(seur$nucleosome_signal[sample_cells]))
  )
  thresh["pct_reads_in_peaks_low"] <- expm1(
    median(log1p(seur$pct_reads_in_peaks[sample_cells])) -
      3 * mad(log1p(seur$pct_reads_in_peaks[sample_cells]))
  )

  cellsToRemove <- sample_cells[which(
    seur$nFeature_peaks[sample_cells] < thresh["nFeature_peaks_low"] |
      seur$nFeature_peaks[sample_cells] > thresh["nFeature_peaks_high"]
  )]
  cellsToRemove <- union(
    cellsToRemove,
    sample_cells[which(
      seur$nCount_peaks[sample_cells] < thresh["nCount_peaks_low"] |
        seur$nCount_peaks[sample_cells] > thresh["nCount_peaks_high"]
    )]
  )
  cellsToRemove <- union(
    cellsToRemove,
    sample_cells[which(
      seur$TSS.enrichment[sample_cells] < thresh["TSS.enrichment_low"]
    )]
  )
  cellsToRemove <- union(
    cellsToRemove,
    sample_cells[which(
      seur$nucleosome_signal[sample_cells] > thresh["nucleosome_signal_high"]
    )]
  )
  cellsToRemove <- union(
    cellsToRemove,
    sample_cells[which(
      seur$pct_reads_in_peaks[sample_cells] < thresh["pct_reads_in_peaks_low"]
    )]
  )

  thresh["numCellsRemove"] <- length(cellsToRemove)
  thresh["pctCellsRemove"] <- length(cellsToRemove) / length(sample_cells) * 100
  return(list(threshold = thresh, filter = cellsToRemove))
}

# Calculate thresholds for each sample
all_thresholds <- list()
all_cellsToRemove <- character()

for (sample in samples) {
  if (!is.na(opt$filterfile)) {
    thresholds <- read.csv(opt$filterfile)
    index <- grep("Sample", colnames(thresholds), ignore.case = T)
    if (sum(thresholds[, index] == sample) == 1) {
      thresh_orig <- thresholds[which(thresholds[, index] == sample), ]
      thresh_orig[index] <- NULL

      thresh <- list()
      cellsToRemove <- character()
      sample_cells <- colnames(seur)[seur$Sample == sample]

      for (i in colnames(thresh_orig)) {
        try({
          if (
            length(grep("_high$", i, ignore.case = T, value = T)) > 0 |
              length(grep("_low$", i, ignore.case = T, value = T)) > 0
          ) {
            colname <- (i %>% strsplit(split = "_"))[[1]] %>%
              head(n = -1) %>%
              paste(collapse = "_")
            direction <- (i %>% strsplit(split = "_"))[[1]] %>%
              tail(n = 1) %>%
              str_to_lower
            if (direction == "low") {
              cellsToRemove <- union(
                cellsToRemove,
                sample_cells[which(
                  seur[[colname]][sample_cells, ] < thresh_orig[[i]]
                )]
              )
              thresh[i] <- thresh_orig[i]
            }
            if (direction == "high") {
              cellsToRemove <- union(
                cellsToRemove,
                sample_cells[which(
                  seur[[colname]][sample_cells, ] > thresh_orig[[i]]
                )]
              )
              thresh[i] <- thresh_orig[i]
            }
          }
        })
      }
      thresh["numCellsRemove"] <- length(cellsToRemove)
      thresh["pctCellsRemove"] <- length(cellsToRemove) /
        length(sample_cells) *
        100
    } else {
      result <- defaultThreshold(seur, sample)
      thresh <- result$threshold
      cellsToRemove <- result$filter
    }
  } else {
    result <- defaultThreshold(seur, sample)
    thresh <- result$threshold
    cellsToRemove <- result$filter
  }

  thresh$Sample <- sample
  all_thresholds[[sample]] <- thresh
  all_cellsToRemove <- union(all_cellsToRemove, cellsToRemove)
}

# Combine thresholds into a data frame
thresh_df <- do.call(
  rbind,
  lapply(all_thresholds, function(x) as.data.frame(x))
)
write.csv(thresh_df, file.path(opt$output, "cell_filter_info.csv"), row.names = FALSE)

## ----Pre-Filter ATAC Violin Plot-------
doVlnPlot <- function(aspect, seur, thresh, sample_name) {
  seur_subset <- subset(seur, subset = Sample == sample_name)
  # Only use group.by if there are multiple samples
  if (length(unique(seur$Sample)) > 1) {
    temp_plot <- VlnPlot(seur_subset, group.by = "Sample", features = aspect) +
      NoLegend() +
      ggtitle(sample_name)
  } else {
    temp_plot <- VlnPlot(seur_subset, features = aspect) +
      NoLegend() +
      ggtitle(sample_name)
  }
  if (
    length(grep(paste0("^", aspect, "_low$"), names(thresh), ignore.case = T)) >
      0
  ) {
    try(
      temp_plot <- temp_plot +
        geom_hline(
          yintercept = thresh[grep(
            paste0("^", aspect, "_low$"),
            names(thresh),
            ignore.case = T
          )][[1]],
          linetype = "dashed"
        )
    )
  }
  if (
    length(grep(
      paste0("^", aspect, "_high$"),
      names(thresh),
      ignore.case = T
    )) >
      0
  ) {
    try(
      temp_plot <- temp_plot +
        geom_hline(
          yintercept = thresh[grep(
            paste0("^", aspect, "_high$"),
            names(thresh),
            ignore.case = T
          )][[1]],
          linetype = "dashed"
        )
    )
  }
  return(temp_plot)
}

for (sample in samples) {
  plots <- sapply(
    c("nFeature_peaks", "nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
    function(x) {
      doVlnPlot(
        aspect = x,
        seur = seur,
        thresh = all_thresholds[[sample]],
        sample_name = sample
      )
    }
  )

  png(
    file.path(opt$output, paste0("PreFilter_VlnPlot_ATAC_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  do.call("grid.arrange", c(plots, nrow = 1))
  dev.off()

  figures[[paste0("PreFilter_VlnPlot_ATAC_", sample)]] <- do.call(
    "grid.arrange",
    c(plots, nrow = 1)
  )

  plots <- sapply(c("pct_reads_in_peaks", "blacklist_ratio"), function(x) {
    doVlnPlot(
      aspect = x,
      seur = seur,
      thresh = all_thresholds[[sample]],
      sample_name = sample
    )
  })

  png(
    file.path(opt$output, paste0("PreFilter_VlnPlot_QC_", sample, ".png")),
    height = 7,
    width = 5,
    units = "in",
    res = 300
  )
  do.call("grid.arrange", c(plots, nrow = 1))
  dev.off()
}

## ----TSS Enrichment and Fragment Length Plots----
for (sample in samples) {
  seur_subset <- subset(seur, subset = Sample == sample)

  png(
    file.path(opt$output, paste0("PreFilter_TSS_Enrichment_Plot_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  # Only use group.by if there are multiple samples
  if (length(samples) > 1) {
    print(
      TSSPlot(seur_subset, group.by = "Sample") +
        ggtitle(paste0("TSS Enrichment - ", sample)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
  } else {
    print(
      TSSPlot(seur_subset) +
        ggtitle(paste0("TSS Enrichment - ", sample)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
  }
  dev.off()

  png(
    file.path(opt$output, paste0("PreFilter_Fragment_Length_Plot_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  
  # Use genome-wide region for FragmentHistogram
  # For single sample, omit group.by to avoid faceting issues
  if (length(samples) > 1) {
    print(FragmentHistogram(object = seur_subset, region = genome_region, group.by = "Sample"))
  } else {
    print(FragmentHistogram(object = seur_subset, region = genome_region))
  }
  dev.off()
}

## ----Pre-Filter UMAP Plot-------
cat("Generating pre-filter UMAP...\n")
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff = "q0")
seur <- RunSVD(seur)
seur <- RunUMAP(object = seur, reduction = "lsi", dims = 2:30)
seur <- FindNeighbors(object = seur, reduction = "lsi", dims = 2:30)
seur <- FindClusters(
  object = seur,
  verbose = FALSE,
  algorithm = 3,
  resolution = 0.8
)

png(file.path(opt$output, "PreFilter_UMAP_ATAC.png"), width = 1800, height = 1600, res = 300)
DimPlot(seur, reduction = "umap", label = TRUE, group.by = "Sample") +
  ggtitle("Pre-Filter UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(file.path(opt$output, "PreFilter_UMAP_ATAC_Clusters.png"), width = 1800, height = 1600, res = 300)
DimPlot(seur, reduction = "umap", label = TRUE) +
  ggtitle("Pre-Filter UMAP - Clusters") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

for (sample in samples) {
  sample_cells <- colnames(seur)[seur$Sample == sample]
  filtered_cells <- intersect(sample_cells, all_cellsToRemove)
  kept_cells <- setdiff(sample_cells, filtered_cells)
  other_cells <- setdiff(colnames(seur), sample_cells)

  png(
    file.path(opt$output, paste0("PreFilter_UMAP_ATAC_Filter_", sample, ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(
    DimPlot(
      seur,
      reduction = "umap",
      label = TRUE,
      cells.highlight = list(
        "Filtered Cells" = filtered_cells,
        "Kept Cells" = kept_cells,
        "Other Samples" = other_cells
      )
    ) +
      ggtitle(paste0("Pre-Filter UMAP - Filtered Cells - ", sample)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(
        labels = c("Other Samples", "Kept Cells", "Filtered Cells"),
        values = c("grey90", "grey30", "#DE2D26")
      )
  )
  dev.off()

  figures[[paste0("PreFilter_UMAP_ATAC_Filter_", sample)]] <- DimPlot(
    seur,
    reduction = "umap",
    label = TRUE,
    cells.highlight = list(
      "Filtered Cells" = filtered_cells,
      "Kept Cells" = kept_cells,
      "Other Samples" = other_cells
    )
  ) +
    ggtitle(paste0("Pre-Filter UMAP - Filtered Cells - ", sample)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(
      labels = c("Other Samples", "Kept Cells", "Filtered Cells"),
      values = c("grey90", "grey30", "#DE2D26")
    )
}

figures$PreFilter_UMAP_ATAC <- DimPlot(
  seur,
  reduction = "umap",
  label = TRUE,
  group.by = "Sample"
) +
  ggtitle("Pre-Filter UMAP") +
  theme(plot.title = element_text(hjust = 0.5))

## ----Pre-Filter QC Feature Plots----
for (sample in samples) {
  seur_subset <- subset(seur, subset = Sample == sample)

  png(
    file.path(opt$output, paste0("PreFilter_FeaturePlot_Counts_", sample, ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c("nCount_peaks", "nFeature_peaks"),
    ncol = 2
  ))
  dev.off()

  png(
    file.path(opt$output, paste0("PreFilter_FeaturePlot_QC_", sample, ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c(
      "TSS.enrichment",
      "nucleosome_signal",
      "pct_reads_in_peaks",
      "blacklist_ratio"
    ),
    ncol = 2
  ))
  dev.off()

  figures[[paste0("PreFilter_FeaturePlot_Counts_", sample)]] <- FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c("nCount_peaks", "nFeature_peaks"),
    ncol = 2
  )
  figures[[paste0("PreFilter_FeaturePlot_QC_", sample)]] <- FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c(
      "TSS.enrichment",
      "nucleosome_signal",
      "pct_reads_in_peaks",
      "blacklist_ratio"
    ),
    ncol = 2
  )
}

## ----Depth vs TSS enrichment plot----
png(
  file.path(opt$output, "PreFilter_DepthCor_Plot.png"),
  height = 7,
  width = 10,
  units = "in",
  res = 300
)
DepthCor(seur)
dev.off()

cat(paste0("Filtering ", length(all_cellsToRemove), " cells...\n"))
seur <- subset(seur, cells = all_cellsToRemove, invert = T)

## ----Post-Filter Feature Plot----
for (sample in samples) {
  sample_cells <- colnames(seur)[seur$Sample == sample]
  if (length(sample_cells) == 0) {
    next
  }

  seur_subset <- subset(seur, subset = Sample == sample)
  # Only use group.by if there are multiple samples
  if (length(samples) > 1) {
    plot1 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nCount_peaks",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot2 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nCount_peaks",
      feature2 = "nFeature_peaks"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot3 <- FeatureScatter(
      seur_subset,
      group.by = "Sample",
      feature1 = "nucleosome_signal",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
  } else {
    plot1 <- FeatureScatter(
      seur_subset,
      feature1 = "nCount_peaks",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot2 <- FeatureScatter(
      seur_subset,
      feature1 = "nCount_peaks",
      feature2 = "nFeature_peaks"
    ) +
      ggtitle(sample) +
      NoLegend()
    plot3 <- FeatureScatter(
      seur_subset,
      feature1 = "nucleosome_signal",
      feature2 = "TSS.enrichment"
    ) +
      ggtitle(sample) +
      NoLegend()
  }

  png(
    file.path(opt$output, paste0("PostFilter_Feature_Plot_", sample, ".png")),
    height = 5,
    width = 15,
    units = "in",
    res = 300
  )
  print(plot1 + plot2 + plot3)
  dev.off()

  figures[[paste0("PostFilter_Feature_Plot_", sample)]] <- plot1 + plot2
}

## ----Post-Filter ATAC Violin Plot-------
for (sample in samples) {
  sample_cells <- colnames(seur)[seur$Sample == sample]
  if (length(sample_cells) == 0) {
    next
  }

  seur_subset <- subset(seur, subset = Sample == sample)

  png(
    file.path(opt$output, paste0("PostFilter_VlnPlot_ATAC_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  # Only use group.by if there are multiple samples
  if (length(samples) > 1) {
    print(VlnPlot(
      seur_subset,
      group.by = "Sample",
      features = c(
        "nFeature_peaks",
        "nCount_peaks",
        "TSS.enrichment",
        "nucleosome_signal"
      ),
      ncol = 4
    ))
  } else {
    print(VlnPlot(
      seur_subset,
      features = c(
        "nFeature_peaks",
        "nCount_peaks",
        "TSS.enrichment",
        "nucleosome_signal"
      ),
      ncol = 4
    ))
  }
  dev.off()

  figures[[paste0("PostFilter_VlnPlot_ATAC_", sample)]] <- VlnPlot(
    seur_subset,
    features = c(
      "nFeature_peaks",
      "nCount_peaks",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    ncol = 4
  )

  plots <- sapply(c("pct_reads_in_peaks", "blacklist_ratio"), function(x) {
    doVlnPlot(
      aspect = x,
      seur = seur,
      thresh = all_thresholds[[sample]],
      sample_name = sample
    )
  })

  png(
    file.path(opt$output, paste0("PostFilter_VlnPlot_QC_", sample, ".png")),
    height = 7,
    width = 5,
    units = "in",
    res = 300
  )
  do.call("grid.arrange", c(plots, nrow = 1))
  dev.off()
}

## ----Post-Filter TSS and Fragment Plots----
for (sample in samples) {
  sample_cells <- colnames(seur)[seur$Sample == sample]
  if (length(sample_cells) == 0) {
    next
  }

  seur_subset <- subset(seur, subset = Sample == sample)

  png(
    file.path(opt$output, paste0("PostFilter_TSS_Enrichment_Plot_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  # Only use group.by if there are multiple samples
  if (length(samples) > 1) {
    print(
      TSSPlot(seur_subset, group.by = "Sample") +
        ggtitle(paste0("TSS Enrichment - ", sample)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
  } else {
    print(
      TSSPlot(seur_subset) +
        ggtitle(paste0("TSS Enrichment - ", sample)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
  }
  dev.off()

  png(
    file.path(opt$output, paste0("PostFilter_Fragment_Length_Plot_", sample, ".png")),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
  
  # Use genome-wide region for FragmentHistogram
  # For single sample, omit group.by to avoid faceting issues
  if (length(samples) > 1) {
    print(FragmentHistogram(object = seur_subset, region = genome_region, group.by = "Sample"))
  } else {
    print(FragmentHistogram(object = seur_subset, region = genome_region))
  }
  dev.off()
}

## ----ATAC Normalizing and Clustering----
cat("Running normalization and clustering on filtered data...\n")
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff = "q0")
seur <- RunSVD(seur)
seur <- RunUMAP(object = seur, reduction = "lsi", dims = 2:30)
seur <- FindNeighbors(object = seur, reduction = "lsi", dims = 2:30)

## ----Post-Filter QC Feature Plots----
for (sample in samples) {
  sample_cells <- colnames(seur)[seur$Sample == sample]
  if (length(sample_cells) == 0) {
    next
  }

  seur_subset <- subset(seur, subset = Sample == sample)

  png(
    file.path(opt$output, paste0("PostFilter_FeaturePlot_Counts_", sample, ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c("nCount_peaks", "nFeature_peaks"),
    ncol = 2
  ))
  dev.off()

  png(
    file.path(opt$output, paste0("PostFilter_FeaturePlot_QC_", sample, ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c(
      "TSS.enrichment",
      "nucleosome_signal",
      "pct_reads_in_peaks",
      "blacklist_ratio"
    ),
    ncol = 2
  ))
  dev.off()

  figures[[paste0("PostFilter_FeaturePlot_Counts_", sample)]] <- FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c("nCount_peaks", "nFeature_peaks"),
    ncol = 2
  )
  figures[[paste0("PostFilter_FeaturePlot_QC_", sample)]] <- FeaturePlot(
    seur_subset,
    reduction = "umap",
    features = c(
      "TSS.enrichment",
      "nucleosome_signal",
      "pct_reads_in_peaks",
      "blacklist_ratio"
    ),
    ncol = 2
  )
}

coord <- Embeddings(seur, reduction = "lsi")[, 2:30]
d <- dist(coord, method = "euclidean")
for (resolution in c(0.1, seq(0.2, 1.0, 0.2), 1.5, 2.0)) {
  seur <- FindClusters(seur, resolution = resolution)

  try({
    clusters <- Idents(seur)
    sil <- silhouette(as.numeric(clusters), dist = d)
    write.csv(
      sil,
      file.path(opt$output, paste0("SilhouetteResult_res.", resolution, ".csv")),
      row.names = F,
      quote = F
    )
    pdf(file.path(opt$output, paste0("SilhouettePlot_res.", resolution, ".pdf")))
    print(plot(
      sil,
      col = as.factor(clusters[order(clusters, decreasing = FALSE)]),
      main = paste(
        "Silhouette plot of Signac clustering - resolution ",
        resolution,
        sep = ""
      ),
      lty = 2
    ))
    print(abline(v = mean(sil[, 3]), col = "red4", lty = 2))
    dev.off()
  })
}

## ----Elbow Plot (LSI components)----
png(file.path(opt$output, "ElbowPlot_LSI.png"), height = 7, width = 7, units = "in", res = 300)
ElbowPlot(seur, reduction = "lsi", ndims = 30)
dev.off()

figures$ElbowPlot_LSI <- ElbowPlot(seur, reduction = "lsi")

## ----UMAP Plots----
for (resolution in grep("_res.", colnames(seur@meta.data), value = T)) {
  png(
    file.path(opt$output, paste0("UMAP_", gsub("_snn", "", resolution), ".png")),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(
    DimPlot(seur, reduction = "umap", group.by = resolution, label = TRUE) +
      ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
  figures[[paste0("UMAP_", gsub("_snn", "", resolution))]] <- DimPlot(
    seur,
    reduction = "umap",
    group.by = resolution,
    label = TRUE
  ) +
    ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Also create UMAP colored by Sample
png(file.path(opt$output, "UMAP_by_Sample.png"), width = 1800, height = 1600, res = 300)
DimPlot(seur, reduction = "umap", group.by = "Sample", label = FALSE) +
  ggtitle("UMAP - Colored by Sample") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

figures$UMAP_by_Sample <- DimPlot(
  seur,
  reduction = "umap",
  group.by = "Sample",
  label = FALSE
) +
  ggtitle("UMAP - Colored by Sample") +
  theme(plot.title = element_text(hjust = 0.5))

saveRDS(seur, file.path(opt$output, "seur_cluster.rds"))

writeLines(capture.output(devtools::session_info()), file.path(opt$output, "sessionInfo.txt"))
