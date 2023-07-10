library(TSCAN)
library(Seurat)
library(tidyr)
library(mclust)
library(igraph)
library(ggplot2)
library(tidyverse)

# Modify TSCAN exprmclust function to display Seurat cell embeddings
exprmclust_mod <- function (data, clusternum = 2:9, modelNames = "VVV", 
                            reduce = T, custom_embedding = NULL) 
{
  
  set.seed(12345)
  if (reduce) {
    sdev <- prcomp(t(data), scale = T)$sdev[1:20]
    x <- 1:20
    optpoint <- which.min(sapply(2:10, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(sdev ~ x + x2)$residuals^2)
    }))
    pcadim = optpoint + 1
    tmpdata <- t(apply(data, 1, scale))
    colnames(tmpdata) <- colnames(data)
    tmppc <- prcomp(t(tmpdata), scale = T)
    if (!is.null(custom_embedding)){
      pcareduceres <- custom_embedding
    } else {
      pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
    }
  }
  else {
    pcareduceres <- t(data)
  }
  clusternum <- clusternum[clusternum > 1]
  res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, 
                                 modelNames = "VVV"))
  clusterid <- apply(res$z, 1, which.max)
  clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = res$G)
  for (cid in 1:res$G) {
    clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == 
                                                                cid]), , drop = F])
  }
  dp <- as.matrix(dist(clucenter))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, 
       clucenter = clucenter)
}

# Try displaying branching paths based on start and end clusters 
plotmclust_mod <- function (mclustobj, x = 1, y = 2, MSTorder = NULL, 
                            show_tree = T, show_branches = F,
                            TSCANorder_results = NULL, show_cell_names = T, 
                            cell_name_size = 3, markerexpr = NULL, 
                            use_custom_path = F, start_cluster = 1, 
                            end_cluster = 2)
{
  color_by = "State"
  lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid, 
                                     sample_name = names(mclustobj$clusterid))
  # Color by cluster number
  lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
  
  # Cell Embedding
  S_matrix <- mclustobj$pcareduceres
  pca_space_df <- data.frame(S_matrix[, c(x, y)])
  colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
  pca_space_df$sample_name <- row.names(pca_space_df)
  
  # edge_df is a dataframe with sample name(barcoded cells) pca and state data
  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "sample_name")
  edge_df$markerexpr <- markerexpr[edge_df$sample_name]
  if (!is.null(markerexpr)) {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2, 
                                    size = markerexpr))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
  } else {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, 
                        size = 3)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_tree) {
    # clucenter are the centers of each cluster
    clucenter <- mclustobj$clucenter[, c(x, y)]
    clulines <- NULL
    if (is.null(MSTorder)) {
      # Find all the shortest paths from one cluster to another
      allsp <- shortest.paths(mclustobj$MSTtree)
      
      # Of all the "shortest" paths, which one has the longest distance
      ## clusters will be considered the endpoints for get.shortest.paths
      longestsp <- which(allsp == max(allsp), arr.ind = T)
      
      if(use_custom_path == FALSE){
        # Find shortest paths - default option using longest paths method
        MSTorder <- get.shortest.paths(mclustobj$MSTtree, 
                                       longestsp[1, 1], 
                                       longestsp[1, 2])$vpath[[1]]
        
        # Custom trajectories - user specified
      } else {
        start_cluster <- start_cluster
        MSTorder_mod <- get.all.shortest.paths(mclustobj$MSTtree, start_cluster, 
                                               1:nrow(clucenter))
        
        # Select which custom paths to display - res[m], where m is the endpoint
        end_cluster <- end_cluster
        MSTorder <- MSTorder_mod$res[[end_cluster]]
      }
    }
    if (show_branches == TRUE){
      order_list <- names(TSCANorder_results)
      order_list_processed <- lapply(regmatches(
        order_list, gregexpr("[[:digit:]]+", order_list)), as.numeric)
      for (j in 1:length(order_list_processed)){
        MSTorder <- order_list_processed[[j]]
        clulines = NULL
        # retrieve center coordinates for each cluster
        for (i in 1:(length(MSTorder) - 1)) {
          clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
          ], clucenter[MSTorder[i + 1], ]))
        }
        # annotate columns of clulines
        clulines <- data.frame(x = clulines[, 1], xend = clulines[, 3], 
                               y = clulines[, 2], yend = clulines[, 4])
        
        g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                         y = "y", yend = "yend", size = NULL), 
                              data = clulines, size = 1) 
        
      }} else {
        for (i in 1:(length(MSTorder) - 1)) {
          clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
          ], clucenter[MSTorder[i + 1], ]))
        }
        clulines <- data.frame(x = clulines[, 1], xend = clulines[, 3],
                               y = clulines[, 2], yend = clulines[, 4])
        
        g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                         y = "y", yend = "yend", size = NULL), 
                              data = clulines, size = 1) 
      }
    clucenter <- data.frame(x = clucenter[, 1], 
                            y = clucenter[,2], id = 1:nrow(clucenter))
    
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
    xlab(paste0("PCA_dimension_", x)) +
    ylab(paste0("PCA_dimension_",y)) + 
    theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_blank()) + 
    theme(legend.position = "top", legend.key.size = unit(0.3,"in"),
          legend.text = element_text(size = 20), 
          legend.title = element_text(size = 20)) + 
    theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(axis.text.x = element_text(size = 17, color = "darkred"), 
          #axis.text.y = element_text(size = 17, color = "black"), 
          axis.title.x = element_text(size = 20, vjust = -1), 
          axis.title.y = element_text(size = 20, vjust = 1), 
          plot.margin = unit(c(1, 1, 1, 1), "cm"))
  g
}

# Bring in object
so_malign <- readRDS("~/Manuscript/Data/ccbr1119_malign.rds")

# Subset by variable features
variable_genes <- VariableFeatures(so_malign)

# Preprocess data
expr_mtx <- as.matrix(so_malign@assays$SCT@data[variable_genes,])
procdata <- expr_mtx

# Remove columns with zero variance
expr_mtx <- expr_mtx[,apply(expr_mtx, 2, var, na.rm=TRUE) != 0]

# Use embedding and cluster information from Seurat
custom_embed <- so_malign@reductions[["pca"]]@cell.embeddings

# Set number of clusters detected by TSCAN (vector of IDs to be used for unsupervised clustering)
#   will be used to construct minimum spanning tree
clus_num <- 1:100

custom_tscan <- exprmclust_mod(data = expr_mtx, clusternum = clus_num, 
                               custom_embedding = NULL)
plot_prelim <- plotmclust_mod(custom_tscan, show_cell_names = F)
print(plot_prelim)

#Selected trajectory differentially expressed genes under qvlue cutoff of 0.05
lpsorder <- TSCANorder(custom_tscan)

diffval <- difftest(procdata,lpsorder)

TSCAN_diff_genes <- diffval[diffval$qval < 0.05,]

# cluster_id
idents <- so_malign@meta.data$custom_cluster[match(names(custom_tscan$clusterid)
                                             ,so_malign@meta.data$Barcode)] 
names(idents) <- so_malign@meta.data$Barcode[match(names(custom_tscan$clusterid),
                                                   so_malign@meta.data$Barcode)]

# use custom annotations from seurat instead of TSCAN defined states
orig_pseudo <- TSCANorder(custom_tscan, flip=TRUE, orderonly=FALSE)

custom_tscan$clusterid <- idents
plot_custom <- plotmclust_mod(custom_tscan, show_cell_names = F)
print(plot_custom)

# Sox9 - Pseudotime scatter plot
gene_expr_df <- as.data.frame.matrix(as.matrix(
  so_malign@assays$SCT@data["Sox9",]))
gene_expr_df$barcode <- rownames(gene_expr_df)

sox9_expr <- gene_expr_df$V1
names(sox9_expr) <- gene_expr_df$barcode

singlegeneplot(sox9_expr, TSCANorder(custom_tscan, flip=TRUE, orderonly=FALSE))

