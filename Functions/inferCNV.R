# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("infercnv")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpmisc)
library(Rfast)
library(quantmod)
library(grid)
library(data.table)
library(dplyr)
library(infercnv)
library(biomaRt)

so_for_cnv <- readRDS("~/Manuscript/Data/Seurat_rds/ccbr1119_merged_processed.rds")

Idents(so_for_cnv) <- so_for_cnv@meta.data$Likely_CellType

Matrix <- as.data.frame(GetAssayData(so_for_cnv,slot='counts'))

# Create annotation table of genes and chromosomal start / end position
biolist <- rownames(Matrix)
ensembl <- useMart("ensembl")

ensembl <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl", 
                      mirror = "useast")
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

t2g <- getBM(attributes=c('mgi_symbol','chromosome_name','start_position',
                          'end_position'), mart = ensembl)
geneInfor <- t2g[match(biolist, t2g$mgi_symbol),]
geneInfor <- drop_na(geneInfor[!duplicated(geneInfor[,1]),])

rownames(geneInfor) <- NULL
colnames(geneInfor) <- NULL

Matrix_final <- Matrix[rownames(Matrix) %in% geneInfor[,1],]

Annotation <- data.frame(v1 = colnames(Matrix_final),
                              v2 = so_for_cnv@meta.data$infercnv_Likely_CellType[
                                so_for_cnv@meta.data$Barcode %in% 
                                  colnames(Matrix_final)])

infercnv_obj=CreateInfercnvObject(raw_counts_matrix=Matrix_final,
                                  annotations_file=Annotation,
                                  gene_order_file=geneInfor,
                                  ref_group_names = c("CD4","CD8","B_cells")
) 

rm(list=setdiff(ls(), c("infercnv_obj")))

infercnv_obj = infercnv::run(infercnv_obj,
                             # use 1 for smart-seq, 0.1 for 10x-genomics
                             cutoff=0.1,  
                             out_dir='CNVresult_CCBR_1119_all',  
                             cluster_by_groups=T,
                             num_threads=6,
                             sd_amplifier = 1.5,
                             denoise=T,
                             HMM=F,no_plot=F)

#########################
# CNV Score Calculation #
#########################
cnv_values <- infercnv_obj@expr.data
cnv_values_scaled <- scales::rescale(as.matrix(cnv_values), to = c(-1,1))

cnv_val_epi <- cnv_values_scaled[,Annotation_file$v1[Annotation_file$v2 == 
                                                       "Epithelial"]]
cnv_score_epi <- apply(cnv_val_epi, MARGIN = 2, function(x) sum((x - mean(x))^2))

cnv_meta_epi <- epithelial_so@meta.data
cnv_meta_epi$cnv_scores <- cnv_score_epi[match(cnv_meta_epi$Barcode, 
                                               names(cnv_score_epi))]

cnv_meta_epi <- cnv_meta_epi %>% mutate(Likely_CellType = case_when(
  pooled_ident_coarse == "Ctrl" ~ "Cholangiocyte",
  cnv_scores > quantile(cnv_scores)[3]  ~ "Malignant",
  TRUE ~ "Cholangiocyte"
))

# # For updating epithelial cell calls in ccbr1119_merged_processed.rds
# write.csv(cnv_meta_epi,"ccbr1119_inferCNV.csv")