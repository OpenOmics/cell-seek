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
              help="Metadata file with information to add to samples"),
  make_option(c("-a", "--adtthresh"), type='double', action='store', default=10,
              help="Lower threshold for filtering ADT features")
)

opt <- parse_args(OptionParser(option_list=option_list))


dir.create(opt$workdir, recursive=TRUE)
setwd(opt$workdir)


rdata <- Read10X(opt$datapath)

seur <- CreateSeuratObject(counts=rdata$`Gene Expression`, project = opt$project)

# Add in ADT assay (trying to filter out potential hashtags by grepping for HTO)
filtered_cite <- list()
adt_thresh <- opt$adtthresh

adt_assay <- CreateAssayObject(counts=rdata$`Antibody Capture`[grep('^HTO[-_]', grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE, invert=TRUE), value=TRUE, ignore.case=TRUE, invert=TRUE),])
filtered_cite[['ADT']] <- names(which(apply(GetAssayData(adt_assay, slot='counts'), 1, max) <= adt_thresh))
adt_names <- names(which(apply(GetAssayData(adt_assay, slot='counts'), 1, max) > adt_thresh))
seur[['ADT']] <- CreateAssayObject(counts=GetAssayData(adt_assay, slot='counts')[adt_names,])

# Add in HTO assay if features with HTO was found
hashtag = FALSE
if (length(as.character(c(grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE), grep('^HTO[-_]', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE)))) > 0) {
  hto_assay <- CreateAssayObject(counts=rdata$`Antibody Capture`[unique(as.character(c(grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE), grep('^HTO[-_]', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE)))),])
  filtered_cite[['HTO']] <- names(which(apply(GetAssayData(hto_assay, slot='counts'), 1, max) <= adt_thresh))
  hto_names <- names(which(apply(GetAssayData(hto_assay, slot='counts'), 1, max) > adt_thresh))
  seur[['HTO']] <- CreateAssayObject(counts=GetAssayData(hto_assay, slot='counts')[hto_names,])
  hashtag = TRUE
}

write.table(adt_thresh, 'CITE_threshold.txt', col.names = FALSE, row.names=FALSE)

ddd <- data.frame(a=I(unlist(lapply(filtered_cite,paste,collapse=","))))
write.table(ddd,file="CITE_excluded.csv", sep=',', quote=FALSE, col.names=FALSE)

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


figures <- list()

## ----Pre-Filter Gene Plot----
seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^[Mm][Tt]-")

plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

plot1 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot2 <- FeatureScatter(seur, group.by='Sample', feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot3 <- FeatureScatter(seur, group.by='Sample', feature1 = "nFeature_RNA", feature2 = "percent.mito") + NoLegend()
png("PreFilter_Gene_Plot.png", height=5, width=10, units='in', res=300)
plot1+plot3+plot2
dev.off()

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
  
#  thresh['nFeature_ADT_low'] <- expm1(median(log1p(seur$nFeature_ADT)) - 3*mad(log1p(seur$nFeature_ADT))) %>% round
#  thresh['nFeature_ADT_high'] <- expm1(median(log1p(seur$nFeature_ADT)) + 3*mad(log1p(seur$nFeature_ADT))) %>% round
  thresh['nCount_ADT_low'] <- expm1(median(log1p(seur$nCount_ADT)) - 3*mad(log1p(seur$nCount_ADT))) %>% round
  thresh['nCount_ADT_high'] <- expm1(median(log1p(seur$nCount_ADT)) + 3*mad(log1p(seur$nCount_ADT))) %>% round
  
  if (hashtag) {
    thresh['nCount_HTO_high'] <- expm1(median(log1p(seur$nCount_HTO)) + 3*mad(log1p(seur$nCount_HTO))) %>% round
  } 
  
  cellsToRemove <- colnames(seur)[which(seur$nFeature_RNA < thresh['nFeature_RNA_low'] | seur$nFeature_RNA > thresh['nFeature_RNA_high'])]
  cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur$nCount_RNA < thresh['nCount_RNA_low'] | seur$nCount_RNA > thresh['nCount_RNA_high'])])
  cellsToRemove <- union(cellsToRemove,  colnames(seur)[which(seur$percent.mito > thresh['percent.mito_high'])])
  
#  cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur$nFeature_ADT < thresh['nFeature_ADT_low'] | seur$nFeature_ADT > thresh['nFeature_ADT_high'])])
  cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur$nCount_ADT < thresh['nCount_ADT_low'] | seur$nCount_ADT > thresh['nCount_ADT_high'])])
  
  if (hashtag) {
    cellsToRemove <- union(cellsToRemove, colnames(seur)[which(seur$nCount_HTO > thresh['nCount_HTO_high'])])
  }
  
  
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

## ----Pre-Filter GEX Violin Plot
plots <- sapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(x) doVlnPlot(aspect=x, seur=seur, thresh=thresh))

png("PreFilter_VlnPlot_RNA.png", height=7, width=7, units='in', res=300)
do.call("grid.arrange", c(plots, nrow=1))
dev.off()

figures$PreFilter_VlnPlot_RNA <- do.call("grid.arrange", c(plots, nrow=1))


## ----Pre-Filter ADT Violin Plot
plots <- sapply(c("nFeature_ADT", "nCount_ADT"), function(x) doVlnPlot(aspect=x, seur=seur, thresh=thresh))

png("PreFilter_VlnPlot_ADT.png", height=7, width=5, units='in', res=300)
do.call("grid.arrange", c(plots, nrow=1))
dev.off()

figures$PreFilter_VlnPlot_ADT <- do.call("grid.arrange", c(plots, nrow=1))


## ----Pre-Filter HTO Violin Plot
if (hashtag) {
  plots <- sapply(c("nFeature_HTO", "nCount_HTO"), function(x) doVlnPlot(aspect=x, seur=seur, thresh=thresh))
  
  png("PreFilter_VlnPlot_HTO.png", height=7, width=5, units='in', res=300)
  do.call("grid.arrange", c(plots, nrow=1))
  dev.off()
  
  figures$PreFilter_VlnPlot_HTO
}

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


## ----Post-Filter ADT Violin Plot

png("PostFilter_VlnPlot_ADT.png", height=7, width=5, units='in', res=300)
VlnPlot(seur, group.by='Sample', features=c("nFeature_ADT", "nCount_ADT"), ncol=2)
dev.off()

figures$PostFilter_VlnPlot_ADT <- VlnPlot(seur, group.by='Sample', features=c("nFeature_ADT", "nCount_ADT"), ncol=2)


## ----Post-Filter HTO Violin Plot
if (hashtag) {
  
  png("PostFilter_VlnPlot_HTO.png", height=7, width=5, units='in', res=300)
  print(VlnPlot(seur, group.by='Sample', features=c("nFeature_HTO", "nCount_HTO"), ncol=2))
  dev.off()
  
  figures$PostFilter_VlnPlot_HTO <- VlnPlot(seur, group.by='Sample', features=c("nFeature_HTO", "nCount_HTO"), ncol=2)
}

## ----Normalize HTO Data----
if (hashtag) {
  hto_quantile <- 0.99
  seur <- NormalizeData(seur, assay = "HTO", normalization.method = "CLR")
  seur <- ScaleData(seur, assay = "HTO", model.use = "linear")
  result <- tryCatch({
    seur <- HTODemux(seur, assay = "HTO", positive.quantile = hto_quantile)
    hashIndex <- 'hash.ID'
    write.table(c("HTODemux", hto_quantile), 'hto_threshold.csv', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    c(seur, hashIndex)
  }, error = function(err){
    #seur <- MULTIseqDemux(seur, assay = "HTO", quantile = hto_quantile, autoThresh=TRUE)
    con <- file('MULTIseqDemux.log')
    sink(con, append=T)
    sink(con, append=T, type="message")
    seur <- MULTIseqDemux(seur, assay = "HTO", autoThresh=TRUE)
    sink()
    sink(type="message")
    hashIndex <- 'MULTI_ID'
    hto_quantile <- tail(str_split(grep("quantile", readLines('MULTIseqDemux.log'), value=T) %>% tail(n=1), pattern=' ')[[1]], n=1)
    write.table(c("MULTIseqDemux", hto_quantile), 'hto_threshold.csv', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    return(c(seur, hashIndex))
  })
  
  seur <- result[[1]]
  hashIndex <- result[[2]]


  entries <- entries <- levels(seur[[hashIndex]][,1])
  seur[[hashIndex]] <- factor(seur[[hashIndex]][,1], levels=c(str_sort(entries[!entries %in% c("Doublet", "Negative")], numeric=T), c("Doublet", "Negative")))
  #seur[[hashIndex]] <- factor(seur[[hashIndex]], levels=levels(seur[[hashIndex]])[order(levels(seur[[hashIndex]]))])

  write.csv(table(seur[[hashIndex]]), paste0('HTO_', hashIndex, '_count.csv'), row.names=F)
}


## ----ADT Normalizing----
DefaultAssay(seur) <- 'ADT'
VariableFeatures(seur) <- rownames(seur[["ADT"]])
seur <- NormalizeData(seur, assay = "ADT", normalization.method = "CLR")
seur <- ScaleData(seur, assay = "ADT", model.use = "linear")
seur <- RunPCA(seur, assay="ADT", reduction.name = 'apca')
seur <- FindNeighbors(seur, dims = 1:min(length(rownames(seur[['ADT']]))-1, 20), reduction = "apca")
seur <- FindClusters(seur, graph.name = "ADT_snn", algorithm = 3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = 'apca', dims = 1:min(length(rownames(seur[['ADT']]))-1, 20), assay = 'ADT', 
                reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

png("UMAP_ADT.png", width=1800, height=1600, res = 300)
DimPlot(seur, reduction='adt.umap', label = TRUE) + ggtitle("ADT")
dev.off()


## ----RNA Normalizing and Clustering----
DefaultAssay(seur) <- 'RNA'
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)
seur <- RunPCA(seur, npcs=50, features = VariableFeatures(object = seur))
seur <- FindNeighbors(seur, dims = 1:30)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA')

coord <- Embeddings(seur, reduction='pca')[,1:30]
d <- dist(coord, method="euclidean")
for(resolution in c(0.1, seq(0.2,1.0,0.2), 1.5, 2.0)){
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

if (hashtag) {
  png(paste0("UMAP_HTO_", hashIndex, ".png"), width=1800, height=1600, res=300)
  print(DimPlot(seur, group.by=hashIndex) + ggtitle(paste0("UMAP ", hashIndex)) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  ncol <- seur[[hashIndex]][,1] %>% levels %>% length() %>% sqrt() %>% ceiling
  nrow <- (seur[[hashIndex]][,1] %>% levels %>% length())/ncol %>% ceiling
  png(paste0("UMAP_HTO_", hashIndex, "-split.png"), width=200+(800*ncol), height=800*nrow, res=300)
  print(DimPlot(seur, split.by=hashIndex, group.by=hashIndex, ncol=ncol) + ggtitle(paste0("UMAP ", hashIndex, " Split")) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

## ----ADT Ridge Plots----

suppressWarnings(dir.create('ADT'))
#pdf('./Ridgeplots.pdf', width=21, height=21)
count <- 1
sil_files <- Sys.glob("SilhouetteResult_res*.csv")
resolutions <- sapply(sil_files, function(x) (x %>% basename %>% tools::file_path_sans_ext() %>% strsplit(split='res.'))[[1]][2]) %>% str_sort(numeric=T)
sil_mean <- sapply(resolutions, function(resolution) {filename <- grep(paste0('res.', resolution, '.csv'), sil_files, value=T); if(length(filename) > 0) { data <- read.csv(filename); mean(data$sil_width)} else{-Inf}})

if (hashtag) {
  group <- hashIndex
} else if(max(sil_mean) > -Inf){
  group <- paste0("RNA_snn_res.", sil_mean %>% which.max %>% names)
}else {
  group <- "RNA_snn_res.0.8"
}
for (i in seq(1,length(names(which(rowSums(seur[['ADT']]) > adt_thresh))), by=25)) {
  png(paste0("ADT/ADT_Ridgeplots_", count, ".png"), width=21, height=4.2*length(i:min(i+24,length(rownames(seur[['ADT']]))))/5, res=300, units='in')
  print(RidgePlot(seur, sort(names(which(rowSums(seur[['ADT']]) > adt_thresh)))[i:min(i+24,length(rownames(seur[['ADT']])))], assay="ADT", ncol=5, group.by=group))
  dev.off()
  count <- count + 1
}

## ----HTO Ridge Plot---
if (hashtag) {
#png('HTO_Ridge_Plot.png', units='in', width=21, height=4.2*length(i:min(i+24,length(rownames(seur[['HTO']]))))/5, res=300)
#png('HTO_Ridge_Plot.png', units='in', width=21, height=4.2*length(1:min(1+24,length(rownames(seur[['HTO']]))))/5, res=300)
png('HTO_Ridge_Plot.png', units='in', width=12, height=9, res=300)
#  for (i in seq(1,length(rownames(seur[["HTO"]])), by=25)) {
print(RidgePlot(seur, sort(rownames(seur[['HTO']]))[1:min(1+24,length(rownames(seur[['HTO']])))], assay="HTO", ncol=min(5, ceiling(sqrt(length(rownames(seur[['HTO']]))))), group.by=hashIndex))
dev.off()
}


saveRDS(seur, 'seur_cluster.rds')
#saveRDS(figures, 'seur_figures.rds')
                      
## ----Matrix export----
hto_mat <- GetAssayData(object = seur, assay = "HTO", slot = "data")
adt_mat <- GetAssayData(object = seur, assay = "ADT", slot = "data")
dir.create(file.path(opt$workdir, "cite-seq-matrix"), showWarnings = FALSE)
write.csv(adt_mat, file = file.path(opt$workdir, "cite-seq-matrix", paste0(opts$project, "_ADT_matrix.csv")))
write.csv(hto_mat, file = file.path(opt$workdir, "cite-seq-matrix", paste0(opts$project, "_HTO_matrix.csv")))

writeLines(capture.output(devtools::session_info()), 'sessionInfo.txt')

