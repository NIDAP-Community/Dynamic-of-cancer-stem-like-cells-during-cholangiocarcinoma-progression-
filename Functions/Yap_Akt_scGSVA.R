scGSVA_wrapper <- function(dir){
  
  library(dplyr)
  library(GSVA)
  library(matrixStats)
  library(limma)
  library(stringr)
  library(edgeR)
  library(tibble)
  library(tidyr)
  
  malign <- readRDS(dir)
  
  #Select Sammples to run Differential Expression
  malign_count <- as.data.frame.matrix(as.matrix(malign$SCT@counts))
  malign_count$Gene <- rownames(malign_count)
  
  # Move last column to the start
  malign_count <- malign_count %>%
    select(Gene, everything())
  
  es <- malign_count
  samples <- colnames(es)[-1]
  es %>% dplyr::select(Gene,samples) %>% dplyr::group_by(Gene) %>% dplyr::summarize_all(sum) -> es
  es %>% dplyr::select(Gene) -> gene
  es %>% dplyr::select(samples) %>% as.matrix() -> es.mat
  rownames(es.mat) <- gene$Gene
  
  x <- DGEList(counts=es.mat, genes=gene$Gene)     
  
  #Run Differential Expression Analysis on group:
  df <- malign@meta.data
  df <- df[match(colnames(es.mat),df$Barcode),]
  dm.formula <- as.formula(paste("~0 +", paste(c("custom_cluster")), sep="+", collapse="+"))
  design=model.matrix(dm.formula, df)
  
  colnames(design) <- str_replace_all(colnames(design), "custom_cluster", "")
  
  v <- voom(x,design=design,normalize="quantile",plot=TRUE)
  
  rownames(v$E) <- v$genes$genes
  as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
  var.df = apply(v$E,1,var)
  v = v[var.df >0,]
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  DEG <- topTable(fit, number = length(v$E), adjust="BH") 
  print(DEG)
  genes = rownames(DEG)
  
  constraints = "Mouse"
  v = length(constraints) == 1
  df.msig <- read.csv("Data/Supplementary/ccbr1119_msigDB.csv")
  
  geneset_list.list <- split(df.msig$gene_symbol, df.msig$gene_set_name) 
  
  df.gsva <- GSVA::gsva(es.mat[genes,],geneset_list.list,method="gsva", min.sz = 15, max.sz = 1200,
                        parallel.sz=4L)
  
  as.data.frame(df.gsva) %>% rownames_to_column("Geneset") -> gsva.df
  fit <- lmFit(df.gsva, design)
  
  cm <- makeContrasts(contrasts = c("Tum1-Tum2",
                                    "Tum2-Tum3",
                                    "Tum1-Tum3"),
                      levels=design)
  print(colnames(design))
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  logFC = fit2$coefficients
  colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
  tstat = fit2$t
  colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
  FC = 2^fit2$coefficients
  FC = ifelse(FC<1,-1/FC,FC)
  colnames(FC)=paste(colnames(FC),"FC",sep="_")
  pvalall=fit2$p.value
  colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
  pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
  colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
  
  finalres=as.data.frame(cbind(df.gsva,FC, logFC, tstat, pvalall, pvaladjall))
  finalres %>% rownames_to_column("Geneset") -> finalres
  
  ### Barplots of Stemness Enrichment Scores ###
  stem_es <- finalres[finalres$Geneset == "custom_stemness_spike_in",]
  
  gsva_annot <- malign@meta.data
  
  # Grab just the GSVA enrichment scores for stemness spike in - barplot
  stem_path <- stem_es[1,1:(ncol(stem_es)-15)]
  
  stem_df <- data.frame(Barcode = colnames(stem_path)[-1], gsva_stm_sc = as.numeric(stem_path[,-1]))
  stem_df$custom_cluster <- gsva_annot$custom_cluster[match(stem_df$Barcode,gsva_annot$Barcode)]
  
  Barplt_df <- stem_df %>% group_by(custom_cluster) %>% summarise(Score_GSVA = mean(gsva_stm_sc))
  
  p <- ggplot(data=Barplt_df, aes(x=custom_cluster, y=Score_GSVA, fill = custom_cluster)) +
    geom_bar(stat="identity") + theme_classic() + theme(text = element_text(size = 40))  
  
  plot(p)
}
  