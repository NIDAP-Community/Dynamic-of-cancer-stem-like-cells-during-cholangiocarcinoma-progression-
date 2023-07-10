#' @title Compute ModScore
#' @description Returns Seurat-class object with metadata containing 
#'              ModuleScores and Likely_CellType calls
#' @details Analyzed features are binned based on averaged expression; 
#'          control features are randomly selected from each bin.
#' 
#' @param object Seurat-class object
#' @param samples.subset List of samples to subset the data by
#' @param sample.to.display List of samples to depict on dimension plot, 
#'                          samples not in the list would be colored gray in 
#'                          the background
#' @param marker.table Table of marker genes for each celltype 
#'                          (column names of the table), append "_prot" or 
#'                          "_neg" for proteins or negative markers
#' @param cite.seq Set to TRUE if there are CITE-seq markers in 
#'                          marker.table (Default: FALSE)
#' @param celltypes Vector of celltypes from marker.table to 
#'                             screen for
#' @param threshold Specify bimodal thresholds for cell classification, 
#'                         should be of the same length as celltypes 
#'                         vector
#' @param general.class Base population of cells to classify
#' @param multi.lvl Toggle to TRUE if there are subpopulations of cells 
#'                          you want to screen for (Default: FALSE)
#' @param lvl.df Table of subpopulation levels and parent-child 
#'                         information (e.g. Tcells-CD4, Tcells-CD8)
#' @param reduction Choose among tsne, umap, and pca (Default: tsne)
#' @param nbins Number of bins for storing control features and analyzing 
#'              average expression (Default: 10)
#' @param gradient.ft.size Set size of axis labels on gradient 
#'                                   density plot of ModuleScore distribution
#'                                   (Default: 6)
#' @param violin.ft.size Set size of axis labels on violin plot of 
#'                             ModuleScore distribution (Default: 6)
#' @param step.size Set step size of distribution plots (Default: 0.1)
#' 
#' @import Seurat
#' @import tidyverse
#' @import gridExtra
#' @import quantmod
#' @import grid
#' @import data.table
#' @import utils
#' @importFrom dplyr select
#'   
#' @export
#' @example Do not run: moduleScore(object = seurat,
#'                                  samples.subset = c("mouse1","mouse2"),
#'                                  sample.to.display = c("mouse1"),
#'                                  marker.table = immuneCellMarkers,
#'                                  celltypes = c("CD4_T","Treg",Monocytes"),
#'                                  threshold = c(0.1,0.4, 0.3),
#'                                  multi.lvl = FALSE
#'                                  )
#'                                  
#' @example Do not run: moduleScore(object = seurat,
#'                                  samples.subset = c("mouse1","mouse2"),
#'                                  sample.to.display = c("mouse1"),
#'                                  marker.table = immuneCellMarkers,
#'                                  celltypes = c("CD4_T","Treg",Monocytes"),
#'                                  threshold = c(0.1,0.4, 0.3),
#'                                  general.class = c("CD_T","Monocytes"),
#'                                  multi.lvl = TRUE,
#'                                  lvl.df = parentChildTable
#'                                  )

#' @return List containing annotated dimension plot with ModuleScore 
#'         distribution of cell marker gene, Seurat Object with cell 
#'         classification metadata

modScore <- function(object, 
                     samples.subset, 
                     sample.to.display, 
                     marker.table, 
                     celltypes,
                     threshold = c(0), 
                     general.class, 
                     multi.lvl = FALSE, 
                     lvl.df, 
                     nbins = 10){

  library(Seurat)
  library(dplyr)
  library(gridExtra)
  library(quantmod)
  library(grid)
  library(data.table)
  library(utils)
  library(dplyr)
  
  # Helper Functions
  # Give each cell an identity based on modulescores and bimodal threshold
  .modScoreCall <- function(ms.meta,threshold,reject){
    
    thres.ls <- list()
    for (i in 1:ncol(ms.meta)){
      thres.ls[[i]]<- rep(threshold[i],nrow(ms.meta))
    }
    thres.df <- data.frame(matrix(unlist(thres.ls),nrow = nrow(ms.meta)))
    
    thres.filter <- ms.meta > thres.df
    ms.meta.filt <- ms.meta * thres.filter
    
    # Find column number with highest modscore
    max.col.vec <- max.col(ms.meta.filt)
    
    # If a row contains all zeroes, they will be labeled with unknown
    zero.filt <- as.integer(!apply(ms.meta.filt, 1,
                                         function(find_zero_rows) 
                                           all(find_zero_rows == 0)))
    
    # Final filtering: 
    final.filt <- (max.col.vec * zero.filt) + 1
    
    # Adjust rejected cells based on upper level
    append.name <- c(reject, names(ms.meta))
    
    # Added the names into a Likely_CellType Column
    dupl.data <- ms.meta
    dupl.data[,"Likely_CellType"] <- append.name[final.filt]
    return(dupl.data)
  }
  
  # Main Code Block
  marker.tab <- unlist(marker.table)
  
  # Create a Barcode column if none is detected
  if (!"Barcode" %in% colnames(object@meta.data)){
    object@meta.data$Barcode <- rownames(object@meta.data)
  }
  
  if (length(samples.subset) == 0) {
    samples.subset = unique(object@meta.data$sample.name)
  }
  
  colnames(object@meta.data) <- gsub("orig_ident","orig.ident",
                                     colnames(object@meta.data))
  
  if("active.ident" %in% slotNames(object)){
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name)=names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  } else {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name)=names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  } 
  
  # Remove original unprocessed object
  rm(object)
  
  # Retrive markers from list
  marker = select(marker.table, celltypes)
  marker.list = as.list(marker)
  
  # Error checking and messages 
  if (sum(unlist(marker.list) %in% rownames(object.sub@assays$SCT@data)) == 0){
    stop("No genes from list was found in data")
  }
  
  if (length(threshold) != length(celltypes)){
    if (sum(threshold) == 0){
      threshold <- rep(0, length(celltypes))
      print("Manual threshold set to zero - outputing preliminary data")
    } else {
      stop("Threshold length does not match # celltypes to analyze")
    }}
  
  names(threshold) <- celltypes
  
  # Check if all markers for particular celltype is present
  figures <- list()
  exclude_cells <- c()
  
  h = 0
  j = 1
  
  for (h in seq_along(marker.list)) {
    print(names(marker.list[h]))
    
    present=lapply(marker.list[[h]], function(x) x %in% rownames(object.sub))
    
    absentgenes = unlist(marker.list[[h]])[present==FALSE];absentgenes=
      absentgenes[is.na(absentgenes)==F]
    
    presentgenes = unlist(marker.list[[h]])[present==TRUE];presentgenes=
      presentgenes[is.na(presentgenes)==F]
    
    print(paste0("Genes not present: ",paste0(absentgenes,collapse=",")))
    print(paste0("Genes present: ",paste0(presentgenes,collapse=",")))
    
    if(length(presentgenes) == 0){
      print(paste0(names(marker.list[h]), 
                   " genes were not found in object and will not be analyzed"))
      exclude_cells[j] <- h
      j = j + 1
    }}  
  
  if (length(exclude_cells) > 0){
    marker.list <- marker.list[-exclude_cells]} else {
      marker.list <- marker.list
    }  
  
  # Calculate modulescores for each celltype in marker table
  for (i in seq_along(marker.list)) { 
    object.sub=AddModuleScore(object.sub,marker.list[i],
                              name = names(marker.list[i]),
                              nbin = nbins)
    
    m = paste0(names(marker.list[i]),"1")
    object.sub@meta.data[[m]] <- scales::rescale(object.sub@meta.data[[m]], 
                                                 to=c(0,1))
    }
  
  # Get rid of "1" at the end of MS columns
  colnames(object.sub@meta.data)[colnames(object.sub@meta.data) %in% paste0(
    names(marker.list),1)] <- names(marker.list)
  
  # Analyze and plot only the class of cells found in the metadata
  general.class <- general.class[general.class %in% 
                                   colnames(object.sub@meta.data)]
  
  # Subset the columns of the metadata containing module scores only
  trunc.meta.gen <- object.sub@meta.data[general.class]
  
  gen.thrs.vec <- threshold[general.class]
  
  # Set elements below threshold to zero. Keep elements above threshold
  call.res <- .modScoreCall(trunc.meta.gen,gen.thrs.vec,reject = "unknown")
  call.res$Barcode <- rownames(call.res)
  
  # Hierarchical Classification
  if (multi.lvl){   
    
    for (k in 1:ncol(lvl.df)){ 
      
      # Temporarily keep results from subpopulation calls
      sub.class.call <- list()
      
      ## Subclass Identification
      # Remove any NAs from comparisons
      store.sub.class <- lvl.df[[k]][!is.na(lvl.df[[k]])]
      
      parent.class <- unique(gsub("(.*)-(.*)","\\1",store.sub.class))
      parent.class <- parent.class[parent.class != ""]
      
      for (parent in parent.class){
        sub.class <- store.sub.class[grepl(parent,store.sub.class)]
        children_class <- gsub("(.*)-(.*)","\\2",sub.class)
        
        # Subset out cells predicted to be parents
        parents <- call.res$Barcode[call.res$Likely_CellType == parent]
        trunc.meta.parent <- object.sub@meta.data[parents,] %>% 
          select(children_class)
        
        # Create output table without information of parent population
        trunc.meta.no.parent <- call.res[!call.res$Likely_CellType == parent,]
        non.parent <- rownames(trunc.meta.no.parent)
        
        # Repeat Module Score Comparison and Cell Prediction with Child Subset:
        child.thres.vec <- threshold[children_class]
        
        sub.class.call[[match(parent,parent.class)]] <- .modScoreCall(
          trunc.meta.parent,child.thres.vec,reject = parent) %>% 
          select(Likely_CellType)
      }
      
      # Reappend subclassification results back to general output 
      sub.class.call <- do.call(rbind,sub.class.call)
      sub.class.call$Barcode <- rownames(sub.class.call)
      
      # Update Likely_CellType column in Calls Output 
      call.res$temp.call <- sub.class.call$Likely_CellType[
        match(call.res$Barcode,sub.class.call$Barcode)]
      
      call.res <- call.res %>% mutate(Likely_CellType = case_when(
        is.na(temp.call) ~ Likely_CellType,
        TRUE ~ temp.call
      ))
      
      # Remove temp.call
      call.res$temp.call <- NULL
    }}
  
  # Updating CellType(s) in metadata with subclass calls
  object.sub@meta.data$Likely_CellType <- call.res$Likely_CellType[
    match(object.sub@meta.data$Barcode,call.res$Barcode)]
  
  return(object.sub)
 }
 