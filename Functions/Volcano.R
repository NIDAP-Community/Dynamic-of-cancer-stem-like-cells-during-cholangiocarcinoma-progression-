
Volcano <- function(df,
                    label.col = "",
                    sig.col = "",
                    pCutoff  = 0.001,
                    lfc.col = "",
                    FCcutoff = 0.5){
  
  library(stringr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  
  df$delabel <- NA
  df$diffexpressed <- "NO"
  
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  df$diffexpressed[df[,lfc.col] > FCcutoff & df[,sig.col] < pCutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df[,lfc.col] < -FCcutoff & df[,sig.col] < pCutoff] <- "DOWN"
  
  df$delabel[df$diffexpressed != "NO"] <- df[,label.col][df$diffexpressed != "NO"]
  
  # plot adding up all layers we have seen so far
  
  volc_plot <- ggplot(data=df, aes(x=eval(parse(text = lfc.col)),
                      y=-log10(eval(parse(text = sig.col))), col=diffexpressed,
                      label=delabel)) +
    geom_point() + 
    ggtitle(gsub(".*_(?=Patient)", "", 
                 colnames(df)[grepl("avg_log2FC", colnames(df))], 
                 perl = TRUE)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    xlab(colnames(df)[grepl( "avg_log2FC", colnames(df))]) +
    ylab("log10_pval") +
    xlab("logFC_PTCy") +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")

  return(volc_plot)
}
