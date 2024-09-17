library(Seurat)
library(ggplot2)
library(gridExtra)
library(SingleR)
library(scRNAseq)
library(scater)
library(cluster)
library(optparse)
library(dplyr)
library(stringr)

option_list <- list(
  make_option(c("-w", "--workdir"), type='character', action='store', default=NA,
              help="Path to the working directory"),
  make_option(c("-d", "--datapath"), type='character', action='store', default=NA,
              help="Path to the 10X output filtered features directory"),
  make_option(c("-s", "--sample"), type='character', action='store', default="sample",
              help="Sample name to be stored within Seurat metadata"),
  make_option(c("-p", "--project"), type='character', action='store', default="project",
              help="Project name to be stored within Seurat metadata"),
  make_option(c("-f", "--filterfile"), type='character', action='store', default=NA,
              help="CSV file containing filters to be applied"),
  make_option(c("-m", "--metadata"), type='character', action='store', default=NA,
              help="Metadata file with information to add to samples")
)

opt <- parse_args(OptionParser(option_list=option_list))

data <- Read10X(opt$datapath)

seur <- CreateSeuratObject(counts=data, project = opt$project)

seur$Sample <- opt$sample

if (!is.na(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
  index <- which(metadata[,grep('sample', colnames(metadata), ignore.case=T)] == opt$sample)
  if (length(index) > 0) {
    for (header in grep('sample', colnames(metadata), ignore.case=T, invert=T, value=T)) {
      seur[[header]] <- metadata[index,header]
    }
  }
}

dir.create(opt$workdir, recursive=TRUE)
setwd(opt$workdir)

figures <- list()

## ----Pre-Filter Gene Plot----
seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^[Mm][Tt]-")

plot1 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot2 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot3 <- FeatureScatter(seur, group.by='Sample', feature1 = "nFeature_RNA", feature2 = "percent.mito") + NoLegend()

png("PreFilter_Gene_Plot.png", height=5, width=15, units='in', res=300)
plot1+plot3+plot2
dev.off()

figures$PreFilter_Gene_Plot <- plot1+plot2


## ----Cell Quality Thresholds - Default----
thresh <- list()
defaultThreshold <- function(seur) {
  thresh <- list()
  thresh['nFeature_RNA_low'] <- expm1(median(log1p(seur$nFeature_RNA)) - 3*mad(log1p(seur$nFeature_RNA))) %>% round
  thresh['nFeature_RNA_low'] <- expm1(median(log1p(seur$nFeature_RNA)) - 3*mad(log1p(seur$nFeature_RNA))) %>% round
  thresh['nFeature_RNA_high'] <- expm1(median(log1p(seur$nFeature_RNA)) + 3*mad(log1p(seur$nFeature_RNA))) %>% round
  thresh['nCount_RNA_low'] <- expm1(median(log1p(seur$nCount_RNA)) - 3*mad(log1p(seur$nCount_RNA))) %>% round
  thresh['nCount_RNA_high'] <- expm1(median(log1p(seur$nCount_RNA)) + 3*mad(log1p(seur$nCount_RNA))) %>% round
  thresh['percent.mito_high'] = min(expm1(median(log1p(seur$percent.mito)) + 3*mad(log1p(seur$percent.mito))) %>% round, 100)
  
  cellsToRemove <- colnames(seur)[which(seur$nFeature_RNA < thresh['nFeature_RNA_low'] | seur$nFeature_RNA > thresh['nFeature_RNA_high'])]
  cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur$nCount_RNA < thresh['nCount_RNA_low'] | seur$nCount_RNA > thresh['nCount_RNA_high'])])
  cellsToRemove <- union(cellsToRemove,  colnames(seur)[which(seur$percent.mito > thresh['percent.mito_high'])])
  
  
  thresh['numCellsRemove'] <- length(cellsToRemove)
  thresh['pctCellsRemove'] <- length(cellsToRemove) / dim(seur)[2] * 100
  return(list(threshold=thresh, filter=cellsToRemove))
}

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
        if (length(grep('_high$', i, ignore.case=T, value=T)) > 0 | length(grep('_low$', i, ignore.case=T, value=T)) > 0) {
          colname <- (i %>% strsplit(split='_'))[[1]] %>% head(n=-1) %>% paste(collapse='_')
          direction <- (i %>% strsplit(split='_'))[[1]] %>% tail(n=1) %>% str_to_lower
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
    thresh['numCellsRemove'] <- length(cellsToRemove)
    thresh['pctCellsRemove'] <- length(cellsToRemove) / dim(seur)[2] * 100
  } else {
    result <- defaultThreshold(seur)
    thresh <- result$threshold
    cellsToRemove <- result$filter
  }
} else {
  result <- defaultThreshold(seur)
  thresh <- result$threshold
  cellsToRemove <- result$filter
}

write.csv(thresh, 'cell_filter_info.csv', row.names=FALSE)

## ----Pre-Filter RNA Violin Plot-------
doVlnPlot <- function(aspect, seur, thresh) {
  temp_plot <- VlnPlot(seur, group.by='Sample', features=aspect) + NoLegend()
  if (length(grep(paste0('^', aspect, '_low$'), names(thresh), ignore.case=T)) > 0) {
    try(
      temp_plot <- temp_plot + geom_hline(yintercept=thresh[grep(paste0('^', aspect, '_low$'), names(thresh), ignore.case=T)][[1]], linetype="dashed")
    )
  }
  if (length(grep(paste0('^', aspect, '_high$'), names(thresh), ignore.case=T)) > 0) {
    try(
      temp_plot <- temp_plot + geom_hline(yintercept=thresh[grep(paste0('^', aspect, '_high$'), names(thresh), ignore.case=T)][[1]], linetype="dashed")
    )
  }
  return(temp_plot)
}

plots <- sapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(x) doVlnPlot(aspect=x, seur=seur, thresh=thresh))

png("PreFilter_VlnPlot_RNA.png", height=7, width=7, units='in', res=300)
do.call("grid.arrange", c(plots, nrow=1))
dev.off()

figures$PreFilter_VlnPlot_RNA <- do.call("grid.arrange", c(plots, nrow=1))

## ----Pre-Filter UMAP Plot-------
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)
seur <- RunPCA(seur, features = VariableFeatures(object = seur))
seur <- FindNeighbors(seur, dims = 1:30)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA')
seur <- FindClusters(seur, resolution = 0.8, algorithm=3, verbose = FALSE)

png("PreFilter_UMAP_RNA.png", width=1800, height=1600, res = 300)
DimPlot(seur, reduction='umap', label = TRUE) + ggtitle("Pre-Filter UMAP") + theme(plot.title = element_text(hjust=0.5))
dev.off()

png("PreFilter_UMAP_RNA_Filter.png", width=1800, height=1600, res = 300)
DimPlot(seur, reduction='umap', label = TRUE, cells.highlight=list("Filtered Cells"=cellsToRemove)) + ggtitle("Pre-Filter UMAP - Filtered Cells") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(labels = c("Kept Cells", "Filtered Cells"), values = c("grey", "#DE2D26"))
dev.off()

figures$PreFilter_UMAP_RNA <- DimPlot(seur, reduction='umap', label = TRUE) + ggtitle("Pre-Filter UMAP") + theme(plot.title = element_text(hjust=0.5))
figures$PreFilter_UMAP_RNA_Filter <- DimPlot(seur, reduction='umap', label = TRUE, cells.highlight=list("Filtered Cells"=cellsToRemove)) + ggtitle("Pre-Filter UMAP - Filtered Cells") + theme(plot.title = element_text(hjust = 0.5))


seur <- subset(seur, cells = cellsToRemove, invert=T)

## ----Post-Filter Gene Plot----
#Post-Filter Plots
plot1 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot2 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot3 <- FeatureScatter(seur, group.by='Sample', feature1 = "nFeature_RNA", feature2 = "percent.mito") + NoLegend()
png("PostFilter_Gene_Plot.png", height=5, width=10, units='in', res=300)
plot1+plot3+plot2
dev.off()

figures$PostFilter_Gene_Plot <- plot1+plot2


## ----Post-Filter RNA Violin Plot-------
png("PostFilter_VlnPlot_RNA.png", height=7, width=7, units='in', res=300)
VlnPlot(seur, group.by='Sample', features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

figures$PostFilter_VlnPlot_RNA <- VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

## ----RNA Normalizing and Clustering----
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)
seur <- RunPCA(seur, npcs=50, features = VariableFeatures(object = seur))
seur <- FindNeighbors(seur, dims = 1:30)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA')


coord <- Embeddings(seur, reduction='pca')[,1:30]
d <- dist(coord, method="euclidean")
for(resolution in c(seq(0.2,1.0,0.2), 1.5, 2.0)){
  seur <- FindClusters(seur, resolution = resolution)
  
  #Calculate silhouette scores and generate plots
  try({
    clusters <- Idents(seur)
    sil<-silhouette(as.numeric(clusters), dist=d)  
    pdf(paste0("SilhouettePlot_res.",resolution,".pdf"))
    print(plot(sil, col=as.factor(clusters[order(clusters, decreasing=FALSE)]), main=paste("Silhouette plot of Seurat clustering - resolution ", resolution, sep=""), lty=2))
    print(abline(v=mean(sil[,3]), col="red4", lty=2))
    dev.off()
    write.csv(sil, paste0('SilhouetteResult_res.', resolution, '.csv'), row.names=F, quote=F)
  })
}

## ----Elbow Plot----
png("ElbowPlot.png", height=7, width=7, units='in', res=300)
ElbowPlot(seur, ndims=50)
dev.off()

figures$ElbowPlot <- ElbowPlot(seur)

## ----UMAP Plots----
for (resolution in grep('_res.', colnames(seur@meta.data), value=T)) {
  png(paste0("UMAP_", gsub("_snn", "", resolution), ".png"), width=1800, height=1600, res = 300)
  print(DimPlot(seur, reduction='umap', group.by=resolution, label=TRUE) + ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  figures[[paste0("UMAP_", gsub("_snn", "", resolution))]] <- DimPlot(seur, reduction='umap', group.by=resolution, label=TRUE) + ggtitle(paste0("UMAP ", gsub("_snn_", " - ", resolution))) + theme(plot.title = element_text(hjust = 0.5))
}


saveRDS(seur, 'seur_cluster.rds')
#saveRDS(figures, 'seur_figures.rds')

writeLines(capture.output(devtools::session_info()), 'sessionInfo.txt')

