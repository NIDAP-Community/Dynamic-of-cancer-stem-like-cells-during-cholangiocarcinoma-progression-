recluster <- function(SO,
                      prepend_text = "old",
                      old_columns_to_save = c(),
                      number_of_pcs = 20,
                      cluster_resolution_low_range = 0.2,
                      cluster_resolution_high_range = 1.2,
                      cluster_resolution_range_bins = 0.2,
                      reduction = "tsne"){
  
library(Seurat)
library(dplyr)
library(cowplot)

## --------------- ##
## Main Code Block ##
## --------------- ##

meta <- SO@meta.data

## Check columns.
if(length(old_columns_to_save) > 0 & !all(old_columns_to_save %in% colnames(SO[[]]))){
  colnames(SO[[]]) <- gsub("\\.","_",colnames(SO[[]]))
  if(!all(old_columns_to_save %in% colnames(SO[[]]))){
    stop("Could not find requested metadata columns!")
  }
}

## Get columns.
for(i in seq_along(old_columns_to_save)){
  old_column_name <- old_columns_to_save[i]
  new_colume_name <- paste(prepend_text,old_column_name, sep="_")
  SO@meta.data[[new_colume_name]] <- SO@meta.data[[old_column_name]]
}

## Remove original cluster columns because they are inaccurate
columns_to_remove <- c("seurat_clusters",grep("^SCT_snn_res",colnames(SO[[]]), value=TRUE))
SO@meta.data %>% select(-one_of(columns_to_remove)) -> SO@meta.data

## Find new clusters.
SO <- FindVariableFeatures(object = SO, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)

# Method for automatically picking NPCS
if(number_of_pcs == 0){
  SO <- RunPCA(object = SO, npcs = 30, verbose = FALSE,seed.use = 42) # initial run
  sumpcsd = sum(SO@reductions$pca@stdev)
  pcvar = (SO@reductions$pca@stdev/sumpcsd)*100
  cumu <- cumsum(pcvar)
  co1 <- which(cumu > 80 & pcvar < 5)[1]
  co2 <- sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 0.1), decreasing = T)[1] + 1
  number_of_pcs = min(co1,co2)
  print(number_of_pcs)
}

## Dimensionality reduction
SO <- RunPCA(object = SO, npcs = number_of_pcs, verbose = FALSE,seed.use = 42)
SO <- RunUMAP(object = SO, reduction = "pca", dims = 1:number_of_pcs, seed.use=42)
SO <- RunTSNE(object = SO, reduction = "pca", dim.embed = 2, dims = 1:number_of_pcs, seed.use = 1)
SO <- FindNeighbors(SO, dims = 1:number_of_pcs)

## Find Clusters.
resolutions <- seq(cluster_resolution_low_range, cluster_resolution_high_range, cluster_resolution_range_bins)
for (r in resolutions) {
  SO <- FindClusters(SO, resolution = r, algorithm = 1)
}
print("Clustering successful!")

## Fix orig_ident back to orig.ident.
colnames(SO@meta.data)[colnames(SO@meta.data) == "orig_ident"] <- "orig.ident"

plot.list <- lapply(paste0("SCT_snn_res.",resolutions),function(z) DimPlot(SO, reduction=reduction, group.by=z) +
                      labs(title=z)+
                      theme(plot.title = element_text(hjust = 0.5)))

g <- plot_grid(plotlist = plot.list)
print(g)

return(SO)
}
