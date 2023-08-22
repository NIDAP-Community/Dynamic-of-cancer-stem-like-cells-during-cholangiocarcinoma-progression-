### Load general functions for making figures ###
source("~/Manuscript/Functions/plotMeta.R")
source("~/Manuscript/Functions/Recluster.R")
source("~/Manuscript/Functions/ModuleScore.R")
source("~/Manuscript/Functions/Color_by_Gene.R")
source("~/Manuscript/Functions/Violin_Plots_by_Metadata.R")

## Run code for preprocessing h5 files (Filter, PCA Normalization, Combine Renormalize)
SO <- source("~/Manuscript/CCBR_1119/Functions/Preprocessing.R")

## Cell Classification / update metadata of already preprocessed seurat ##
SO <- modScore(object = SO, 
        samples.subset = unique(SO$orig.ident), 
        sample.to.display = unique(SO$orig.ident), 
        marker.table = read.csv("~/Manuscript/Data/Supplementary/ccbr1119_markers.csv"), 
        celltypes = c("Epithelial","Hepatocytes","B_cells","T_cells","NK_cells",
                      "TAM","Dendritic","Endothelial","TAM_M1","TAM_M2","CD8",
                      "CD4","Treg","Fibroblasts"),
        threshold = c(0.06,0.055,0.2,0.2,0.2,0.2,0.25,0.25,0.122,0.15,0.21,0.25,
                      0.1,0.04), 
        general.class = c("Epithelial","Hepatocytes","B_cells","T_cells",
                          "NK_cells","TAM","Dendritic","Endothelial",
                          "Fibroblasts"), 
        multi.lvl = TRUE, 
        lvl.df = read.csv("~/Manuscript/Data/Supplementary/ccbr1119_levels.csv"), 
        nbins = 24)

# From running inferCNV.R
cnv_res <- read.csv("~/Manuscript/Data/Supplementary/ccbr1119_inferCNV.csv")

so_metadata <- SO@meta.data

so_metadata$cnv_calls <- cnv_res$Likely_CellType[match(
  so_metadata$Barcode,cnv_res$Barcode)]

# Metadata preprocessing - Append cnv calls for epithelial cells, set timepoints
so_metadata <- so_metadata %>% mutate(Likely_CellType2 = case_when(
  Likely_CellType == "Epithelial" ~ cnv_calls,
  TRUE ~ Likely_CellType
)) %>% mutate(pooled_ident_pre = gsub(".{1}$","",orig.ident)) %>% 
  mutate(pooled_ident =case_when(
  pooled_ident_pre == "Eight_week" ~ "Week_8",
  pooled_ident_pre == "Five_week" ~ "Week_5",
  TRUE ~ "Ctrl"
))

SO$Likely_CellType <- so_metadata$Likely_CellType2
SO$pooled_ident_pre <- so_metadata$pooled_ident_pre
SO$pooled_ident <- so_metadata$pooled_ident

cell_cnt <- table(SO@meta.data$Likely_CellType)
SO@meta.data$celltype_w_cnt <- paste(SO@meta.data$Likely_CellType, 
                                     cell_cnt[match(SO@meta.data$Likely_CellType, 
                                                    names(cell_cnt))], sep = " ")

keep_cells <- setdiff(SO$Barcode,SO$Barcode[is.na(SO$Likely_CellType)])
SO <- subset(SO, cells = keep_cells)

### Supplementary Figure 2 ###
Idents(SO) <- SO$pooled_ident_pre

ctrl_so <- subset(SO, idents = c("Ctrl","Ctrl_Eight_week","Ctrl_Five_week"))

ctrl_so <- recluster(ctrl_so,
                     prepend_text = "old",
                     old_columns_to_save = c(),
                     number_of_pcs = 20,
                     reduction = "pca")

ctrl_so$ctrl_pca <- case_when(ctrl_so$pooled_ident_pre == "Ctrl_Five_week" ~ "Week_5_Ctrl",
                              ctrl_so$pooled_ident_pre == "Ctrl_Eight_week" ~ "Week_8_Ctrl",
                              TRUE ~ "Ctrl")

plotMetadata(ctrl_so, metadata.to.plot = "ctrl_pca",
             reduction.type = "pca")

# Filter out Ctrl 5wk and Ctrl 8wk, since they will not be used for downstream analysis 
Idents(SO) <- SO$orig.ident
SO <- subset(SO, idents = c("Ctrl1","Ctrl2","Ctrl3","Eight_week1",
                            "Eight_week2","Eight_week3","Five_week1",
                            "Five_week2","Five_week3"))

### Figure 1B ###
plotMetadata(object = SO,
             metadata.to.plot = c("celltype_w_cnt","pooled_ident"),
             columns.to.summarize = "c()",
             summarization.cut.off = 5,
             reduction.type = "umap")

### Table 1C ###
table(SO$Likely_CellType)

### Figure 1D ###
colorByGene(object = SO,
            gene = c("Krt19","Fabp1","Cd79a","Cd3d","Foxp3",
                     "Cd7","Cd14","Clec9a","Eng","Acta2"),
            reduction.type = "umap",
            number.of.rows = 2,
            assay = "SCT")

## Filter down to Malignant Cells ##
Idents(SO) <- SO$Likely_CellType
malign <- subset(SO, idents = "Malignant")

### Recluster Malignant Cluster with SCT_snn_res.0.1 ###
malign <- recluster(malign,
          prepend_text = "old",
          old_columns_to_save = c(),
          number_of_pcs = 20,
          cluster_resolution_low_range = 0.1,
          cluster_resolution_high_range = 0.4,
          cluster_resolution_range_bins = 0.1,
          reduction = "pca")

# Merged clusters 0 and 1 due to high stemness GSVA scores 
malign_meta <- malign@meta.data

malign_meta$SCT_snn_res.0.1_new <- as.character(malign_meta$SCT_snn_res.0.1)

malign_meta <- malign_meta %>% mutate(custom_cluster = case_when(
  SCT_snn_res.0.1_new == "0" | SCT_snn_res.0.1_new == "3" ~ "Tum1",
  SCT_snn_res.0.1_new == "1" ~ "Tum2",
  SCT_snn_res.0.1_new == "2" ~ "Tum3"))

malign$custom_cluster <- malign_meta$custom_cluster[match(
  malign@meta.data$Barcode, malign_meta$Barcode)]

Tm4sf1_sct <- malign$SCT@scale.data["Tm4sf1",]
malign$Tm4sf1_expr <- Tm4sf1_sct[match(malign_meta$Barcode,names(Tm4sf1_sct))]

### Figure 2A ###
plotMetadata(object = malign,
             metadata.to.plot = c("custom_cluster","pooled_ident"),
             reduction.type = "pca")

### Figure 2B ###
table(malign$custom_cluster, malign$pooled_ident)

### Supplementary Figure 4 ###
colorByGene(object = malign,
            gene = c("Cat","Atp1b1","Ccnd1","Cxcl12"),
            reduction.type = "pca",
            number.of.rows = 0,
            return.seurat.object = FALSE,
            color = "red",
            point.size = 1,
            point.shape = 16,
            point.transparency = 0.5,
            use.cite.seq.data = FALSE,
            assay = "SCT")

### Supplementary Figure 5A ###
Idents(malign) <- malign$custom_cluster
malign_Tum1 <- subset(malign, ident = "Tum1")

### SCT_snn_res.0.4 ###
malign_Tum1 <- recluster(malign_Tum1,
          prepend_text = "old",
          old_columns_to_save = c(),
          number_of_pcs = 20,
          cluster_resolution_low_range = 0.1,
          cluster_resolution_high_range = 0.4,
          cluster_resolution_range_bins = 0.1,
          reduction = "pca")

### Figure 2C, 3E, Supplementary Figure 6A and 6B ~ GSVA ###
source("~/Manuscript/Functions/scGSVA_wrapper.R")

scGSVA_wrapper(seurat_dir = "~/Manuscript/Data/Seurat_rds/ccbr1119_malign.rds",
               msigDB_dir = "~Manuscript/Data/Supplementary/ccbr1119_msigDB.csv",
               clust_ident = "custom_cluster")
scGSVA_wrapper(seurat_dir = "~/Manuscript/Data/Supplementary/ccbr1119_kras_supplement.rds",
               msigDB_dir = "~Manuscript/Data/Supplementary/ccbr1119_msigDB.csv",
               clust_ident = "Tum_clus")
scGSVA_wrapper(seurat_dir = "~/Manuscript/Data/Supplementary/ccbr1119_yap_supplement.rds",
               msigDB_dir = "~Manuscript/Data/Supplementary/ccbr1119_msigDB.csv",
               clust_ident = "Tum_clus")

### Figure 2D Weighted Coexpression ###
source("~/Manuscript/Functions/scWCGNA.R")
wcgna_res <- scWCGNA_wrapper(so)

### Figure 2D Kaplan-Meier Plots, use corresponding genes from wcgna_res ###
source("~/Manuscript/Functions/Survival_Analysis.R")
survival_wrapper(wcgna_mod = "mod1")
survival_wrapper(wcgna_mod = "mod3")

### Figure 3A and 3B ###
source("~/Manuscript/Functions/TSCAN_Trajectory.R")

### Figure 3C and Supplementary Figure 5B ~ Regulon Expression ###
source("~/Manuscript/Functions/Malign_SCENIC.R")
source("~/Manuscript/Functions/Malign_Tum1_SCENIC.R")

### Violin Plots Figures 3D, 4A, Supp. 6A, 6B
violinPlot(object = so, 
          assay = "SCT", 
          slot = "scale.data", 
          group.by = "custom_cluster", 
          group.subset = c("Tum1","Tum2","Tum3"),
          genes.of.interest = c("Mki67","Tm4sf1"),
          filter.outliers = TRUE)

### Figure 4D ###
source("~/Manuscript/Functions/Scatter_Density.R")

### Figure 5A ###
Idents(SO) <- SO$Likely_CellType
immune_so <- subset(SO, idents = c("B_cells","CD4","CD8","Dendritic","NK_cells",
                            "T_cells","TAM","TAM_M1","TAM_M2",
                            "Treg"))

plotMetadata(object = immune_so,
             metadata.to.plot = c("celltype_w_cnt","pooled_ident"),
             reduction.type = "umap")

### Figure 5B ###
source("~/Manuscript/Functions/Sankey_Plot.R")

### Figure 5C ###
source("~/Manuscript/Functions/CellChat_Indv_Wrapper.R")
cellchat_indv_wrapper(SO,
                      timepoint = "Ctrl",
                      celltypes = c('Endothelial','TAM_M1','B_cells',
                                    'NK_cells','CD4','CD8','Fibroblasts',
                                    'TAM_M2','Dendritic','Hepatocytes',
                                    'Treg','Cholangiocyte'))

cellchat_indv_wrapper(SO,
                      timepoint = "Week_5",
                      celltypes = c('Endothelial','TAM_M1','B_cells',
                                    'NK_cells','CD4','CD8','Fibroblasts',
                                    'TAM_M2','Dendritic','Hepatocytes',
                                    'Treg','Malignant','Cholangiocyte'))

cellchat_indv_wrapper(SO,
                      timepoint = "Week_8",
                      celltypes = c('Endothelial','TAM_M1','B_cells',
                                    'NK_cells','CD4','CD8','Fibroblasts',
                                    'TAM_M2','Dendritic','Hepatocytes',
                                    'Treg','Malignant','Cholangiocyte'))

### Figure 5D ###
source("~/Manuscript/Functions/CellChat_Diff_Wrapper.R")
cellchat_diff_wrapper(SO,
                      time1 = "Ctrl",
                      time2 = "Five")

cellchat_diff_wrapper(SO,
                      time1 = "Ctrl",
                      time2 = "Eight")

cellchat_diff_wrapper(SO,
                      time1 = "Five",
                      time2 = "Eight")

### Figure 5E ###
# see cellphoneDB file (cpdb_script.txt) and supplementary data 
#   (cpdb_LR / cpdb_cell_interactions)

# Generate counts and metadata for cpdb
Idents(SO) <- SO@meta.data$Likely_CellType
SO <- subset(SO, idents = c("Malignant","Endothelial","Fibroblasts",
                            "B_cells","NK_cells","CD4","CD8",
                            "TAM_M1","TAM_M2","Dendritic","Treg"))

# Designate malignant cells as early (high stem) or late (low stem)
#   based on GSVA stemness and pseudotime
meta <- SO@meta.data

meta <- meta %>% mutate(stemness = case_when(
  Likely_CellType == "Malignant" & custom_cluster == "Tum1" ~ "Malign_Early",
  Likely_CellType == "Malignant" & custom_cluster == "Tum3" ~ "Malign_Late"
  TRUE ~ Likely_CellType
))

celltype_annot <- data.frame(Cell = meta$Barcode, cell_type = meta$Likely_CellType)
#write.table(celltype_annot, "CCBR_1119_cpdb_meta.txt", 
#   sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Capitalize genes for CelphoneDB v2.1.7
Counts <- as.data.frame(SO$SCT@scale.data)
rownames(Counts) <- toupper(rownames(Counts))
#write.table(Counts, "CCBR_1119_cpdb_counts.txt", sep = "\t", 
#   col.names = NA, row.names = TRUE, quote = FALSE)

### Supplementary Figure 3 Tum1 vs all other clusters###
source("~/Manuscript/Functions/MAST.R")
source("~/Manuscript/Functions/Volcano.R")

malign <- readRDS("~/Manuscript/Data/Seurat_rds/ccbr1119_malign.rds")
Tum1_v_all <- MAST(SO.sub = malign,
     metadata_table = malign@meta.data, 
     parameter_to_test = "custom_cluster",
     contrasts = c("Tum1-all"),
     test_to_use = "MAST")

Volcano(df = yap_Tum1_v_all,
        label.col = "Gene",
        sig.col = "p_val_Tum1_vs_all",
        pCutoff  = 0.00001,
        lfc.col = "avg_log2FC_Tum1_vs_all",
        FCcutoff = 0.5)

### Supplementary Analysis ###
# Kruskal-Wallis, Dunn PostHoc test for expression differences across Tum subclusters
source("~/Manuscript/Functions/DunnPostHoc.R")

# Map WGCNA module genes to pathways
source("~/Manuscript/Functions/l2p.R")

# Generates table of pathways (along with descriptions) that can be mapped to WGCNA modules
l2p_wrapper(wgcna_module = "1")
l2p_wrapper(wgcna_module = "2")
l2p_wrapper(wgcna_module = "3")



