
MAST <- function(SO.sub,
                 metadata_table, 
                 parameter_to_test = "",
                 contrasts = c(""),
                 test_to_use = "MAST",
                 log_fc_threshold = 0,
                 assay_to_use = "SCT",
                 use_log_2 = TRUE,
                 latent_vars = c()
                 ) {
  
  #define contrasts
  newcont <- list()
  for (i in 1:length(contrasts)){
    newcont[[i]] <- c(paste(unlist(strsplit(contrasts[i],"-"))))
  }
  contrasts <- newcont
  
  #ERROR CATCHING
  #collect valid names of valid columns
  validColumns <- character()
  for (i in colnames(metadata_table)) {
    if (!any(is.na(metadata_table[[i]]))) {
      validColumns <-c(validColumns,i)
    }
  }
  
  param2test <- parameter_to_test
  
  if (param2test =="") {
    mcols = colnames(SO.sub@meta.data)
    param2test <-mcols[grepl("RNA_snn",mcols)][[1]]
    print(paste("No parameter selected, defaulting to",param2test))
  }
  
  contrastTarget <- SO.sub@meta.data[[param2test]]
  contrastType <- param2test
  contrastCounts = as.data.frame(table(contrastTarget))
  validContrasts = subset(contrastCounts, Freq>2)[[1]]
  
  #catch malformed contrasts
  for (i in contrasts) {
    if (!(i[[1]] %in% contrastTarget)) {
      print(paste(i[[1]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[2]] %in% contrastTarget) & (i[[2]] != "all")) {
      print(paste(i[[2]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (length(i)>2) {
      print("Contrasts are as follows..")
      print(i)
      stop("The console says there are too many inputs in your contrasts. A contrast should only contain Group1-Group2, but the console thinks you have inputed Group1-Group2-Group3")
    } else if (!(i[[2]] %in% validContrasts) & (i[[2]] != "all")) {
      print(paste(i[[2]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[1]] %in% validContrasts)) {
      print(paste(i[[1]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    }
  }
  
  #print out contrast cell contrastCounts
  for (i in seq_along(contrasts)) {
    firstGroup <- contrasts[[i]][[1]]
    firstGroupCount <- subset(contrastCounts, contrastTarget == firstGroup)$Freq
    if  (contrasts[[i]][[2]]!= "all") {
      secondGroup <-contrasts[[i]][[2]]
      secondGroupCount <-subset(contrastCounts, contrastTarget == secondGroup)$Freq      
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. cluster",secondGroup,"with",secondGroupCount,"cells."))
    } else {
      secondGroupCount <-ncol(SO.sub)-firstGroupCount
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. all other clusters, totalling",secondGroupCount,"cells."))
    } 
  }
  
  #define and call function for running DEG
  get_deg_table <- function(n) {
    library(Seurat)
    
    firstCluster <-n[1]
    secondCluster <- n[2]
    
    if (n[2]=='all') {
      secondCluster <- NULL
    }
    
    Idents(SO.sub) <- param2test
    
    #workaround for log2/ln changes:
    if (use_log_2) { log_fc_threshold <- log_fc_threshold/log2(exp(1)) }
    
    markers = FindMarkers(SO.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = test_to_use, logfc.threshold = log_fc_threshold, verbose=FALSE, assay = assay_to_use, latent.vars = eval(parse(text = "latent_vars")))
    colnames(markers) <- chartr(old=" ",new="_",paste(colnames(markers), n[1],"vs",n[2],sep = "_"))
    
    if (use_log_2) { markers[grep("avg_logFC_", colnames(markers))] <- markers[grep("avg_logFC_", colnames(markers))]*log2(exp(1)) }
    
    return(markers)
  }
  
  
  deg_tables <- lapply(contrasts, get_deg_table) 
  
  for(i in seq_along(deg_tables)){
    degtab <- deg_tables[[i]]
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos 
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
    print(paste0("The number of upregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(pos[1])
    print(paste0("The number of downregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(neg[1]) 
  }
  
  #Merge the deg tables together
  out_df <- NULL
  for (i in deg_tables) {
    if (is.null(out_df)) {
      out_df <- deg_tables[1]
      out_df <- as.data.frame(out_df)
    } else {
      out_df <- merge(out_df, i, by="row.names", all=TRUE)
      rownames(out_df) <- out_df$Row.names #set the rownames
      out_df$Row.names <- NULL #drop the row.names columns which we no longer need
    }
  }
  
  out_df$Gene <- rownames(out_df)
  out_df$Row.names <- NULL
  out_df <- out_df %>% dplyr::select(Gene, everything())
  return(out_df)
  
}
