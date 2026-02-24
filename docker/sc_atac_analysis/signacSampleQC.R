#!/usr/bin/env Rscript
library(Signac)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(cluster)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(optparse)
library(dplyr)
library(stringr)

options(error = function() traceback(2))

option_list <- list(
  ###### Sample ID ######
  make_option(
    c("-s", "--sample"),
    default = NA,
    type = 'character',
    help = "Sample name to be stored within Seurat metadata"
  ),

  ###### Project Name ######
  make_option(
    c("--project"),
    type = 'character',
    help = "Project name to be stored within Seurat metadata"
  ),

  ###### SAMPLE INFORMATION ######
  # Example **filter file** format:
  #   > Sample,nFeature_RNA_low,nFeature_RNA_high,percent.mito_high
  #   > sample1,500,6000,15
  #   > sample2,500,6000,5
  #   > sample4,500,6000,5
  make_option(
    c("--filterfile"),
    type = 'character',
    default = NA,
    help = "CSV file containing filters to be applied"
  ),
  # Example **metadata** file format:
  #   > Sample,Type,batch
  #   > sample1,tumor,1
  #   > sample2,normal,1
  #   > sample4,tumor,2
  make_option(
    c("--metadata"),
    type = 'character',
    default = NA,
    help = "Metadata file with information to add to samples"
  ),
  # Example **barcodes** file format:
  #   > NO_BARCODE,18170307,843531,2466,3590401,1250975,81911,5063,12395960,0,0,0,0,0,0,0,0,0,0
  #   > AAACGAAAGAAACGCC-1,36,6,0,3,10,3,0,14,0,0,1,0,0,0,1,0,4,8
  #   > AAACGAAAGAAAGCAG-1,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1,1
  #   > AAACGAAAGAAAGGGT-1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
  #   > AAACGAAAGAAATGGG-1,2,0,0,0,0,0,0,2,0,3,1,0,0,0,1,0,1,2
  #   > AAACGAAAGAAATTCG-1,1,0,0,0,0,0,0,1,0,2,1,0,0,0,1,0,0,0
  #   > AAACGAAAGAACAGGA-1,1,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,0
  make_option(
    c("--barcodes"),
    type = "character",
    help = "Barcode file path from cellranger outputs (singlecell.csv)"
  ),
  # Example **fragments** file format: BED5
  #   > chr1    10081   10267   GCTCCTAAGGGTTCCC-1      1
  #   > chr1    10083   10266   TGATTTCCAGATTAAG-1      1
  #   > chr1    10085   10198   AACCAACCAGCTGATT-1      1
  #   > chr1    10085   10267   GGTACCGTCTTAATCC-1      1
  #   > chr1    10090   10204   TAGGAGGAGCTTACCA-1      1
  #   > chr1    10091   10198   CCTCCCTCACCAAGGA-1      1
  #   > chr1    10095   10204   GTTACGAAGATCGATA-1      1
  make_option(
    c("--fragments"),
    type = "character",
    help = "Fragments file path (BED5/6) from cellranger outputs (fragments.tsv.gz)"
  ),

  ###### GENOMIC REFERENCES ######
  # See config/genome.json
  make_option(
    c("--genes"),
    type = "character",
    default = NA,
    help = "Path to GTF annotation file (e.g., gencode.v32.primary_assembly.annotation.gtf.gz)"
  ),
  # See config/genome.json
  make_option(
    c("-g", "--genome"),
    type = "character",
    default = "hg38",
    help = "Genome name (e.g., hg38, mm10) - used for seqlevels style"
  ),
  ###### filtered_peak_bc_matrix.h5 cellranger output ######
  make_option(
    c("-f", "--featurematrix"),
    type = 'character',
    default = NA,
    help = "Path to the 10X feature matrix (from cellranger)"
  ),

  ###### OUTPUT DESIGNATION ######
  make_option(
    c("-o", "--output"),
    type = 'character',
    default = file.path(getwd()),
    help = "Output directory path (default: current working directory)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (file.exists(opt$output) == FALSE) {
  dir.create(opt$output, recursive = TRUE)
}

if (is.na(opt$sample)) {
  stop("Sample name must be provided with -s/--sample")
}

## ----Load Reference Annotations----

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

## ----Load 10X Data----
counts <- Read10X_h5(file.path(opt$featurematrix))
cellranger_metrics <- read.csv(
  file = file.path(opt$barcodes),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = genome(annotations)[1],
  fragments = file.path(opt$fragments),
  min.cells = 10,
  min.features = 200
)

seur <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = opt$project,
  meta.data = cellranger_metrics
)

seur$Sample <- opt$sample

## ----Add Sample Metadata----
if (!is.na(opt$metadata)) {
  sample_metadata <- read.csv(opt$metadata)
  index <- which(
    sample_metadata[, grep(
      'sample',
      colnames(sample_metadata),
      ignore.case = T
    )] ==
      opt$sample
  )
  if (length(index) > 0) {
    for (header in grep(
      'sample',
      colnames(sample_metadata),
      ignore.case = T,
      invert = T,
      value = T
    )) {
      seur[[header]] <- sample_metadata[index, header]
    }
  }
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

figures <- list()

## ----Calculate QC Metrics----
# Compute nucleosome signal score per cell
seur <- NucleosomeSignal(object = seur)

# Compute TSS enrichment score per cell
seur <- TSSEnrichment(object = seur, fast = FALSE)

# Add blacklist ratio
if (opt$genome == "hg38" || opt$genome == "hg2024") {
  seur$blacklist_ratio <- (seur$blacklist_region_fragments / seur$peak_region_fragments) * 100
} else if (opt$genome == "mm10" || opt$genome == "mm2024") {
  seur$blacklist_ratio <- (seur$blacklist_region_fragments / seur$peak_region_fragments) * 100
} else {
  warning("Genome not recognized for blacklist ratio calculation. Setting blacklist_ratio to NA.")
  seur$blacklist_ratio <- NA
}

# Calculate fraction of reads in peaks (FRiP)
seur$pct_reads_in_peaks <- seur$peak_region_fragments /
  seur$passed_filters *
  100

## ----Pre-Filter Gene Plot----
plot1 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nCount_peaks",
  feature2 = "TSS.enrichment"
) +
  NoLegend()
plot2 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nCount_peaks",
  feature2 = "nFeature_peaks"
) +
  NoLegend()
plot3 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nucleosome_signal",
  feature2 = "TSS.enrichment"
) +
  NoLegend()

png(
  file.path(opt$output, "PreFilter_Feature_Plot.png"),
  height = 5,
  width = 15,
  units = 'in',
  res = 300
)
print(plot1 + plot2 + plot3)
dev.off()

figures$PreFilter_Feature_Plot <- plot1 + plot2

## ----Cell Quality Thresholds - Default----
thresh <- list()
defaultThreshold <- function(seur) {
  thresh <- list()
  thresh['nFeature_peaks_low'] <- expm1(
    median(log1p(seur$nFeature_peaks)) - 3 * mad(log1p(seur$nFeature_peaks))
  ) %>%
    round
  thresh['nFeature_peaks_high'] <- expm1(
    median(log1p(seur$nFeature_peaks)) + 3 * mad(log1p(seur$nFeature_peaks))
  ) %>%
    round
  thresh['nCount_peaks_low'] <- expm1(
    median(log1p(seur$nCount_peaks)) - 3 * mad(log1p(seur$nCount_peaks))
  ) %>%
    round
  thresh['nCount_peaks_high'] <- expm1(
    median(log1p(seur$nCount_peaks)) + 3 * mad(log1p(seur$nCount_peaks))
  ) %>%
    round
  thresh['TSS.enrichment_low'] <- max(
    expm1(
      median(log1p(seur$TSS.enrichment)) - 3 * mad(log1p(seur$TSS.enrichment))
    ),
    1
  )
  thresh['nucleosome_signal_high'] <- expm1(
    median(log1p(seur$nucleosome_signal)) +
      3 * mad(log1p(seur$nucleosome_signal))
  )
  thresh['pct_reads_in_peaks_low'] <- expm1(
    median(log1p(seur$pct_reads_in_peaks)) -
      3 * mad(log1p(seur$pct_reads_in_peaks))
  )

  cellsToRemove <- colnames(seur)[which(
    seur$nFeature_peaks < thresh['nFeature_peaks_low'] |
      seur$nFeature_peaks > thresh['nFeature_peaks_high']
  )]
  cellsToRemove <- union(
    cellsToRemove,
    colnames(seur)[which(
      seur$nCount_peaks < thresh['nCount_peaks_low'] |
        seur$nCount_peaks > thresh['nCount_peaks_high']
    )]
  )
  cellsToRemove <- union(
    cellsToRemove,
    colnames(seur)[which(seur$TSS.enrichment < thresh['TSS.enrichment_low'])]
  )
  cellsToRemove <- union(
    cellsToRemove,
    colnames(seur)[which(
      seur$nucleosome_signal > thresh['nucleosome_signal_high']
    )]
  )
  cellsToRemove <- union(
    cellsToRemove,
    colnames(seur)[which(
      seur$pct_reads_in_peaks < thresh['pct_reads_in_peaks_low']
    )]
  )

  thresh['numCellsRemove'] <- length(cellsToRemove)
  thresh['pctCellsRemove'] <- length(cellsToRemove) / dim(seur)[2] * 100
  return(list(threshold = thresh, filter = cellsToRemove))
}

# Calculate thresholds for the sample provided via command line
if (!is.na(opt$filterfile)){
  thresholds <- read.csv(opt$filterfile)
  index <- grep("Sample", colnames(thresholds), ignore.case=T)
  if (sum(thresholds[,index] == opt$sample) == 1) {
    thresh_orig <- thresholds[which(thresholds[,index] == opt$sample),]
    thresh_orig[index] <- NULL
    
    thresh <- list()
    cellsToRemove <- character()
    
    for (i in colnames(thresh_orig)) {
      try({
        if (length(grep("_high$", i, ignore.case=T, value=T)) > 0 | length(grep("_low$", i, ignore.case=T, value=T)) > 0) {
          colname <- (i %>% strsplit(split="_"))[[1]] %>% head(n=-1) %>% paste(collapse="_")
          direction <- (i %>% strsplit(split="_"))[[1]] %>% tail(n=1) %>% str_to_lower
          if(direction == "low") {
            cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur[[colname]] < thresh_orig[[i]])])
            thresh[i] <- thresh_orig[i]
          }
          if (direction == "high") {
            cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur[[colname]] > thresh_orig[[i]])])
            thresh[i] <- thresh_orig[i]
          }
        }
      })
    }
    thresh["numCellsRemove"] <- length(cellsToRemove)
    thresh["pctCellsRemove"] <- length(cellsToRemove) / ncol(seur) * 100
    thresh["Sample"] <- opt$sample
  } else {
    result <- defaultThreshold(seur)
    thresh <- result$threshold
    cellsToRemove <- result$filter
    thresh["Sample"] <- opt$sample
  }
} else {
  result <- defaultThreshold(seur)
  thresh <- result$threshold
  cellsToRemove <- result$filter
  thresh["Sample"] <- opt$sample
}

# Write threshold information to CSV
if(!is.null(thresh) && length(thresh) > 0) {
  write.csv(as.data.frame(thresh), file.path(opt$output, "cell_filter_info.csv"), row.names=FALSE)
}

## ----Pre-Filter ATAC Violin Plot-------
doVlnPlot <- function(aspect, seur, thresh) {
  temp_plot <- VlnPlot(seur, group.by = 'Sample', features = aspect) +
    NoLegend()
  if (
    length(grep(paste0('^', aspect, '_low$'), names(thresh), ignore.case = T)) >
      0
  ) {
    try(
      temp_plot <- temp_plot +
        geom_hline(
          yintercept = thresh[grep(
            paste0('^', aspect, '_low$'),
            names(thresh),
            ignore.case = T
          )][[1]],
          linetype = "dashed"
        )
    )
  }
  if (
    length(grep(
      paste0('^', aspect, '_high$'),
      names(thresh),
      ignore.case = T
    )) >
      0
  ) {
    try(
      temp_plot <- temp_plot +
        geom_hline(
          yintercept = thresh[grep(
            paste0('^', aspect, '_high$'),
            names(thresh),
            ignore.case = T
          )][[1]],
          linetype = "dashed"
        )
    )
  }
  return(temp_plot)
}

plots <- sapply(
  c("nFeature_peaks", "nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
  function(x) doVlnPlot(aspect = x, seur = seur, thresh = thresh)
)

png(
  file.path(opt$output, "PreFilter_VlnPlot_ATAC.png"),
  height = 7,
  width = 10,
  units = 'in',
  res = 300
)
do.call("grid.arrange", c(plots, nrow = 1))
dev.off()

figures$PreFilter_VlnPlot_ATAC <- do.call("grid.arrange", c(plots, nrow = 1))
plots <- sapply(c("pct_reads_in_peaks", "blacklist_ratio"), function(x) {
  doVlnPlot(aspect = x, seur = seur, thresh = thresh)
})

png(
  file.path(opt$output, "PreFilter_VlnPlot_QC.png"),
  height = 7,
  width = 5,
  units = 'in',
  res = 300
)
do.call("grid.arrange", c(plots, nrow = 1))
dev.off()


## ----TSS Enrichment and Fragment Length Plots----
png(
  file.path(opt$output, "PreFilter_TSS_Enrichment_Plot.png"),
  height = 7,
  width = 10,
  units = 'in',
  res = 300
)
print(TSSPlot(seur, group.by = 'Sample') +
  ggtitle("TSS Enrichment") +
  theme(plot.title = element_text(hjust = 0.5)))
dev.off()

if (ncol(seur) > 0) {
  tryCatch({
    png(
      file.path(opt$output, "PreFilter_Fragment_Length_Plot.png"),
      height = 7,
      width = 10,
      units = 'in',
      res = 300
    )
    # Only create fragment object if it doesn't exist
    if (length(Fragments(seur)) == 0) {
      Fragments(seur) <- CreateFragmentObject(
        path = opt$fragments,
        cells = colnames(seur),
        validate.fragments = FALSE
      )
    }
    # Use genome-wide region for FragmentHistogram
    print(FragmentHistogram(seur, region = genome_region))
    dev.off()
  }, error = function(e) {
    warning("Failed to create pre-filter fragment histogram: ", e$message)
  })
} else {
  warning("No cells available. Skipping pre-filter fragment histogram.")
}

## ----Pre-Filter UMAP Plot-------
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff = 'q0')
seur <- RunSVD(seur)
seur <- RunUMAP(object = seur, reduction = 'lsi', dims = 2:30)
seur <- FindNeighbors(object = seur, reduction = 'lsi', dims = 2:30)
seur <- FindClusters(
  object = seur,
  verbose = FALSE,
  algorithm = 3,
  resolution = 0.8
)

png(
  file.path(opt$output, "PreFilter_UMAP_ATAC.png"),
  width = 1800,
  height = 1600,
  res = 300
)
print(DimPlot(seur, reduction = 'umap', label = TRUE) +
  ggtitle("Pre-Filter UMAP") +
  theme(plot.title = element_text(hjust = 0.5)))
dev.off()

figures$PreFilter_UMAP_ATAC <- DimPlot(seur, reduction = 'umap', label = TRUE) +
  ggtitle("Pre-Filter UMAP") +
  theme(plot.title = element_text(hjust = 0.5))

png(
  file.path(opt$output, "PreFilter_UMAP_ATAC_Filter.png"),
  width = 1800,
  height = 1600,
  res = 300
)
print(DimPlot(
  seur,
  reduction = 'umap',
  label = TRUE,
  cells.highlight = list("Filtered Cells" = cellsToRemove)
) +
  ggtitle("Pre-Filter UMAP - Filtered Cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    labels = c("Kept Cells", "Filtered Cells"),
    values = c("grey", "#DE2D26")
  ))
dev.off()

figures$PreFilter_UMAP_ATAC_Filter <- DimPlot(
  seur,
  reduction = 'umap',
  label = TRUE,
  cells.highlight = list("Filtered Cells" = cellsToRemove)
) +
  ggtitle("Pre-Filter UMAP - Filtered Cells") +
  theme(plot.title = element_text(hjust = 0.5))

## ----Pre-Filter QC Feature Plots----
png(
  file.path(opt$output, "PreFilter_FeaturePlot_Counts.png"),
  width = 1800,
  height = 1600,
  res = 300
)
FeaturePlot(
  seur,
  reduction = 'umap',
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2
)
dev.off()

png(
  file.path(opt$output, "PreFilter_FeaturePlot_QC.png"),
  width = 1800,
  height = 1600,
  res = 300
)
FeaturePlot(
  seur,
  reduction = 'umap',
  features = c(
    "TSS.enrichment",
    "nucleosome_signal",
    "pct_reads_in_peaks",
    "blacklist_ratio"
  ),
  ncol = 2
)
dev.off()

figures$PreFilter_FeaturePlot_Counts <- FeaturePlot(
  seur,
  reduction = 'umap',
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2
)
figures$PreFilter_FeaturePlot_QC <- FeaturePlot(
  seur,
  reduction = 'umap',
  features = c(
    "TSS.enrichment",
    "nucleosome_signal",
    "pct_reads_in_peaks",
    "blacklist_ratio"
  ),
  ncol = 2
)

## ----Depth vs TSS enrichment plot----
png(
  file.path(opt$output, "PreFilter_DepthCor_Plot.png"),
  height = 7,
  width = 10,
  units = 'in',
  res = 300
)
DepthCor(seur)
dev.off()

seur <- subset(seur, cells = cellsToRemove, invert = T)

## ----Post-Filter Feature Plot----
plot1 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nCount_peaks",
  feature2 = "TSS.enrichment"
) +
  NoLegend()
plot2 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nCount_peaks",
  feature2 = "nFeature_peaks"
) +
  NoLegend()
plot3 <- FeatureScatter(
  seur,
  group.by = 'Sample',
  feature1 = "nucleosome_signal",
  feature2 = "TSS.enrichment"
) +
  NoLegend()

png(
  file.path(opt$output, "PostFilter_Feature_Plot.png"),
  height = 5,
  width = 15,
  units = 'in',
  res = 300
)
plot1 + plot2 + plot3
dev.off()

figures$PostFilter_Feature_Plot <- plot1 + plot2

## ----Post-Filter ATAC Violin Plot-------
try({
  png(
    file.path(opt$output, "PostFilter_VlnPlot_ATAC.png"),
    height = 7,
    width = 10,
    units = 'in',
    res = 300
  )
  print(VlnPlot(
    seur,
    group.by = 'Sample',
    features = c(
      "nFeature_peaks",
      "nCount_peaks",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    ncol = 4
  ))
  dev.off()

  figures$PostFilter_VlnPlot_ATAC <- VlnPlot(
    seur,
    features = c(
      "nFeature_peaks",
      "nCount_peaks",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    ncol = 4
  )
})

try({
  plots <- sapply(c("pct_reads_in_peaks", "blacklist_ratio"), function(x) {
    doVlnPlot(aspect = x, seur = seur, thresh = thresh)
  })

  png(
    file.path(opt$output, "PostFilter_VlnPlot_QC.png"),
    height = 7,
    width = 5,
    units = 'in',
    res = 300
  )
  do.call("grid.arrange", c(plots, nrow = 1))
  dev.off()
})

## ----Post-Filter TSS and Fragment Plots----
try({
  png(
    file.path(opt$output, "PostFilter_TSS_Enrichment_Plot.png"),
    height = 7,
    width = 10,
    units = 'in',
    res = 300
  )
  print(TSSPlot(seur, group.by = 'Sample') +
    ggtitle("TSS Enrichment") +
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
})

## ----Post-Filter Fragment Length Plot----
if (ncol(seur) > 0) {
  tryCatch({
    png(
      file.path(opt$output, "PostFilter_Fragment_Length_Plot.png"),
      height = 7,
      width = 10,
      units = 'in',
      res = 300
    )
    # Use genome-wide region for FragmentHistogram
    print(FragmentHistogram(seur, region = genome_region))
    dev.off()
  }, error = function(e) {
    warning("Failed to create post-filter fragment histogram: ", e$message)
  })
} else {
  warning("No cells remaining after filtering. Skipping post-filter fragment histogram.")
}

## ----ATAC Normalizing and Clustering----
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff = 'q0')
seur <- RunSVD(seur)
seur <- RunUMAP(object = seur, reduction = 'lsi', dims = 2:30)
seur <- FindNeighbors(object = seur, reduction = 'lsi', dims = 2:30)

## ----Post-Filter QC Feature Plots----
png(
  file.path(opt$output, "PostFilter_FeaturePlot_Counts.png"),
  width = 1800,
  height = 1600,
  res = 300
)
FeaturePlot(
  seur,
  reduction = 'umap',
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2
)
dev.off()

png(
  file.path(opt$output, "PostFilter_FeaturePlot_QC.png"),
  width = 1800,
  height = 1600,
  res = 300
)
FeaturePlot(
  seur,
  reduction = 'umap',
  features = c(
    "TSS.enrichment",
    "nucleosome_signal",
    "pct_reads_in_peaks",
    "blacklist_ratio"
  ),
  ncol = 2
)
dev.off()

figures$PostFilter_FeaturePlot_Counts <- FeaturePlot(
  seur,
  reduction = 'umap',
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2
)
figures$PostFilter_FeaturePlot_QC <- FeaturePlot(
  seur,
  reduction = 'umap',
  features = c(
    "TSS.enrichment",
    "nucleosome_signal",
    "pct_reads_in_peaks",
    "blacklist_ratio"
  ),
  ncol = 2
)

coord <- Embeddings(seur, reduction = 'lsi')[, 2:30]
d <- dist(coord, method = "euclidean")
for (resolution in c(0.1, seq(0.2, 1.0, 0.2), 1.5, 2.0)) {
  seur <- FindClusters(seur, resolution = resolution)
  # Calculate silhouette scores and generate plots
  try({
    clusters <- Idents(seur)
    sil <- silhouette(as.numeric(clusters), dist = d)
    write.csv(
      sil,
      file.path(
        opt$output,
        paste0('SilhouetteResult_res.', resolution, '.csv')
      ),
      row.names = F,
      quote = F
    )
    pdf(file.path(
      opt$output,
      paste0("SilhouettePlot_res.", resolution, ".pdf")
    ))
    plot(
      sil,
      col = as.factor(clusters[order(clusters, decreasing = FALSE)]),
      main = paste(
        "Silhouette plot of Signac clustering - resolution ",
        resolution,
        sep = ""
      ),
      lty = 2
    )
    abline(v = mean(sil[, 3]), col = "red4", lty = 2)
    dev.off()
  })
}

## ----Elbow Plot (LSI components)----
png(
  file.path(opt$output, "ElbowPlot_LSI.png"),
  height = 7,
  width = 7,
  units = 'in',
  res = 300
)
ElbowPlot(seur, reduction = "lsi", ndims = 30)
dev.off()

figures$ElbowPlot_LSI <- ElbowPlot(seur, reduction = "lsi")

## ----UMAP Plots----
for (resolution in grep('_res.', colnames(seur@meta.data), value = T)) {
  png(
    file.path(
      opt$output,
      paste0("UMAP_", gsub("_snn", "", resolution), ".png")
    ),
    width = 1800,
    height = 1600,
    res = 300
  )
  print(
    DimPlot(seur, reduction = 'umap', group.by = resolution, label = TRUE) +
      ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
  figures[[paste0("UMAP_", gsub("_snn", "", resolution))]] <- DimPlot(
    seur,
    reduction = 'umap',
    group.by = resolution,
    label = TRUE
  ) +
    ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) +
    theme(plot.title = element_text(hjust = 0.5))
}

saveRDS(seur, file.path(opt$output, 'seur_cluster.rds'))
# saveRDS(figures, file.path(opt$output, 'seur_figures.rds'))

writeLines(
  capture.output(devtools::session_info()),
  file.path(opt$output, 'sessionInfo.txt')
)
