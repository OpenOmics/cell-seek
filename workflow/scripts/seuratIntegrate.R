library(Seurat)
library(ggplot2)
library(optparse)
library(dplyr)
library(stringr)

option_list <- list(
  make_option(c("-w", "--workdir"), type='character', action='store', default=NA,
              help="Path to the working directory"),
  make_option(c("-r", "--rdsfiles"), type='character', action='store', default=NA,
              help="List of rds files for the samples to integrate")
)

opt <- parse_args(OptionParser(option_list=option_list))


seur_list <- sapply(strsplit(opt$rdsfiles, ',')[[1]], readRDS)
names(seur_list) <- sapply(seur_list, function(x) x$Sample[1])


setwd(opt$workdir)

obj <- merge(seur_list[[1]], seur_list[2:length(seur_list)], add.cell.ids=names(seur_list))

temp <- obj

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "merged_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.merged")

options(future.globals.maxSize = 8000 * 1024^2)

for (method in c('CCAIntegration', 'RPCAIntegration', 'HarmonyIntegration', 'JointPCAIntegration')) {
  try({
  reduction <- paste0("integrated.", gsub('Integration', '', method) %>% tolower())
  clusters <- paste0(gsub('Integration', '', method) %>% tolower(), '_clusters')
  umap <- gsub('integrated', 'umap', reduction)
  obj <- IntegrateLayers(
    object = obj, method = method,
    orig.reduction = "pca", new.reduction = reduction,
    verbose = FALSE
    )
    obj <- FindClusters(obj, resolution = 0.2, cluster.name = clusters)
    obj <- RunUMAP(obj, reduction = reduction, dims = 1:30, reduction.name = umap)
  })
}

saveRDS(obj, 'integrated_norm.rds')

groupwidth <- 2/3
plotwidth <- 19/3
columns <- obj$Sample %>% unique %>% length %>% sqrt %>% ceiling
rows <- ((obj$Sample %>% unique %>% length)/columns) %>% ceiling
for (reduction in obj@reductions %>% names %>% grep(pattern='umap', value=T)) {
  name <- gsub(pattern='umap.', replacement='', x=reduction)
  groupcols <- (obj[[paste0(name, '_clusters')]][,1] %>% unique %>% length/20) %>% ceiling
  ggsave(filename=paste0('UMAP_Norm_', reduction, '_Sample.png'), plot=DimPlot(obj, reduction=reduction, group.by='Sample', label=T), units='in', width=7, height=7, dpi=300)
  ggsave(filename=paste0('UMAP_Norm_', reduction, '_Sample-Split.png'), plot=DimPlot(obj, reduction=reduction, group.by=paste0(name, '_clusters'), split.by='Sample', label=T, ncol=columns), units='in', width=(columns*plotwidth)+groupwidth*groupcols, height=rows*7, dpi=300)
  ggsave(filename=paste0('UMAP_Norm_', reduction, '_Cluster.png'), plot=DimPlot(obj, reduction=reduction, group.by=paste0(name, '_clusters'), label=T), units='in', width=plotwidth+groupwidth*groupcols, height=7, dpi=300)
}

obj <- SCTransform(temp)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "merged_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.merged")

for (method in c('CCAIntegration', 'RPCAIntegration', 'HarmonyIntegration', 'JointPCAIntegration')) {
  try({
  reduction <- paste0("integrated.", gsub('Integration', '', method) %>% tolower())
  clusters <- paste0(gsub('Integration', '', method) %>% tolower(), '_clusters')
  umap <- gsub('integrated', 'umap', reduction)
  obj <- IntegrateLayers(
    object = obj, method = method,
    orig.reduction = "pca", new.reduction = reduction,
    verbose = FALSE, normalization.method='SCT'
    )
    obj <- FindNeighbors(obj, reduction = reduction, dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.2, cluster.name = clusters)
    obj <- RunUMAP(obj, reduction = reduction, dims = 1:30, reduction.name = umap)
  })
}

saveRDS(obj, 'integrated_sct.rds')

columns <- obj$Sample %>% unique %>% length %>% sqrt %>% ceiling
rows <- ((obj$Sample %>% unique %>% length)/columns) %>% ceiling
for (reduction in obj@reductions %>% names %>% grep(pattern='umap', value=T)) {
  name <- gsub(pattern='umap.', replacement='', x=reduction)
  groupcols <- (obj[[paste0(name, '_clusters')]][,1] %>% unique %>% length/20) %>% ceiling
  ggsave(filename=paste0('UMAP_SCT_', reduction, '_Sample.png'), plot=DimPlot(obj, reduction=reduction, group.by='Sample', label=T), units='in', width=7, height=7, dpi=300)
  ggsave(filename=paste0('UMAP_SCT_', reduction, '_Sample-Split.png'), plot=DimPlot(obj, reduction=reduction, group.by=paste0(name, '_clusters'), split.by='Sample', label=T, ncol=columns), units='in', width=(columns*plotwidth)+groupwidth*groupcols, height=rows*7, dpi=300)
  ggsave(filename=paste0('UMAP_SCT_', reduction, '_Cluster.png'), plot=DimPlot(obj, reduction=reduction, group.by=paste0(name, '_clusters'), label=T), units='in', width=plotwidth+groupwidth*groupcols, height=7, dpi=300)
}

writeLines(capture.output(devtools::session_info()), 'seuratIntegrate_sessionInfo.txt')

