
  require(ggplot2)
  require(gridExtra)
  library(Seurat)
  library(dplyr)
  
  so = readRDS("Data/Supplementary/Zhang_supp.rds")
  #print(rownames(so[['Protein']]))
  
  samples = c("X_ICC_18_Tumor","X_ICC_20_Tumor","X_ICC_23_Tumor","X_ICC_24_Tumor1","X_ICC_24_Tumor2")
  
  ## Goal is to have column 1 of the new metadata be named "orig.ident" for downstream compatibility.
  ## Check new metadata for "orig.ident" column, else fix the "orig_ident" column name, else print an error message.
  if ("orig.ident" %in% colnames(so@meta.data)) { ## If orig.ident already is the first column ...
    print("Found orig.ident in column 1 of SO metadata.")
  } else if ("orig_ident" %in% colnames(so@meta.data)) { ## Else if "orig_ident" is the first column ...
    colnames(so@meta.data)[colnames(so@meta.data) == "orig_ident"] <- "orig.ident"
    print("Found orig_ident in column 1 of new metadata table. Changed to orig.ident for downstream compatibility.")
  } else { ## Else print an error message explaining we expect one of the two above as the first column in the new metadata.
    print("ERROR: Found neither orig.ident nor orig_ident in column 1 of new metadata table. Please try again with a new metadata table with one of these as the column name of the first column in the dataframe.")
  }
  
  if("active.ident" %in% slotNames(so)){
    sample_name = as.factor(so@meta.data$orig.ident)
    names(sample_name)=names(so@active.ident)
    so@active.ident <- as.factor(vector())
    so@active.ident <- sample_name
    so.sub = subset(so, ident = samples)
  } else {
    sample_name = as.factor(so@meta.data$orig.ident)
    names(sample_name)=names(so@active.ident)
    so@active.ident <- as.factor(vector())
    so@active.ident <- sample_name
    so.sub = subset(so, ident = samples)
  }
  
  gg.overlay <- function(so.sub,df,marker1,marker2){
    
    df %>% dplyr::arrange(mark1.scale) -> df    
    xmin = min(df$dr1) - 0.1*min(df$dr1)
    xmax = max(df$dr1) + 0.1*min(df$dr1)
    
    gg.z1 <- ggplot(df, aes(dr1,dr2))+
      geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape=16,size=0.5, alpha=0.5)+
      theme_classic() +
      xlab("tsne-1") + 
      ylab("tsne-2") +
      ggtitle(marker1) +
      coord_fixed()
    
    df %>% dplyr::arrange(mark2.scale) -> df
    
    gg.z2 <- ggplot(df, aes(dr1,dr2))+
      geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape=16,size=0.5, alpha=0.5)+
      theme_classic() +
      xlab("tsne-1") + 
      ylab("tsne-2") +
      ggtitle(marker2) +
      coord_fixed()
    
    df %>% dplyr::mutate(avg = mark2.scale+mark1.scale) %>% dplyr::arrange(avg) -> df
    
    gg <- ggplot(df, aes(dr1,dr2))+
      geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape=16,size=0.5, alpha=0.5)+
      theme_classic() +
      xlab("tsne-1") + 
      ylab("tsne-2") +
      ggtitle("Combined") +
      coord_fixed()
    
    return(list(gg.z1,gg.z2,gg))
  }
  
  t1 = 0.4
  t2 = 0.60
  addlines = TRUE
  
  gg.overlay2 <- function(so.sub,df,marker1,marker2){
    
    df %>% dplyr::arrange(mark1.scale) -> df    
    
    # Create unscaled axis labels
    display_unscaled <- TRUE
    
    if(display_unscaled == TRUE){
      label1_min <- paste("unscaled min:", round(min(mark1),digits = 2))
      label1_max <- paste("unscaled max:", round(max(mark1),digits = 2))
      label1 <- paste(as.character(marker1), label1_min, label1_max, sep = "\n")
      
      label2_min <- paste("unscaled min:", round(min(mark2),digits = 2))
      label2_max <- paste("unscaled max:", round(max(mark2),digits = 2))
      label2 <- paste(as.character(marker2), label2_min, label2_max, sep = "\n")} else {
        
        label1 <- as.character(marker1)
        label2 <- as.character(marker2)
      }
    
    gg.z1 <- ggplot(df, aes(mark1.scale,mark2.scale))+
      geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape = 20,size=0.5)+
      theme_classic() +
      xlab(label1) + 
      ylab(label2) +
      coord_fixed()
    
    df %>% dplyr::arrange(mark2.scale) -> df
    
    gg.z2 <- ggplot(df, aes(mark1.scale,mark2.scale))+
      geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape = 20,size=0.5)+
      theme_classic() +
      xlab(label1) + 
      ylab(label2) +
      coord_fixed()
    
    df %>% dplyr::mutate(avg = mark2.scale+mark1.scale) %>% dplyr::arrange(avg) -> df
    
    gg <- ggplot(df, aes(mark1.scale,mark2.scale))+
      geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape = 20,size=0.5)+
      theme_classic() +
      xlab(label1) + 
      ylab(label2) +
      coord_fixed()
    
    if(addlines==TRUE){
      gg.z1 <- gg.z1 + 
        geom_vline(xintercept=t1,linetype="dashed") +
        geom_hline(yintercept=t2,linetype="dashed") 
      gg.z2 <- gg.z2 + 
        geom_vline(xintercept=t1,linetype="dashed") +
        geom_hline(yintercept=t2,linetype="dashed") 
      gg <- gg +
        geom_vline(xintercept=t1,linetype="dashed") +
        geom_hline(yintercept=t2,linetype="dashed") 
    }
    
    return(list(gg.z1,gg.z2,gg))
  }
  
  marker1 <- "TM4SF1"
  marker2 <- "cytotrace"
  scale1 <- TRUE
  scale2 <- TRUE
  
  mark1 = so.sub@assays$SCT_cytotrace@data[marker1,]
  if(scale1){
    q1 = quantile(mark1,0.99)
    q0 = quantile(mark1,1-0.99)
    mark1[mark1<q0]=q0
    mark1[mark1>q1]=q1
  }
  mark1.scale <- scales::rescale(mark1, to=c(0,1))
  
  mark2 = so.sub@assays$SCT_cytotrace@data[marker2,]
  if(scale2){
    q1 = quantile(mark2,0.99)
    q0 = quantile(mark2,1-0.99)
    mark2[mark2<q0]=q0
    mark2[mark2>q1]=q1
  }
  mark2.scale <- scales::rescale(mark2, to=c(0,1))
  
  df <- data.frame(cbind(dr1=so.sub@reductions$tsne@cell.embeddings[,1],
                         dr2=so.sub@reductions$tsne@cell.embeddings[,2],
                         mark1.scale,mark2.scale))
  gg.list <- gg.overlay(so.sub,df,marker1,marker2)
  gg.list2 <- gg.overlay2(so.sub,df,marker1,marker2)
  
  x <- df$mark1.scale
  y <- df$mark2.scale
  
  df_heatmap <- data.frame(x = x, y = y,
                           d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
  
  p <- ggplot(df_heatmap, aes(x, y, col = d)) +
    geom_point(size = 1) +
    scale_color_identity() + xlab(marker1) + ylab(marker2) +
    theme_bw() + geom_smooth(method="auto", se=TRUE, colour="black", size=0.5, fullrange=FALSE, level=0.95)
  
  grid.arrange(gg.list[[1]], gg.list[[2]], gg.list[[3]], gg.list2[[1]],gg.list2[[2]],gg.list2[[3]],p,ncol=3)
  