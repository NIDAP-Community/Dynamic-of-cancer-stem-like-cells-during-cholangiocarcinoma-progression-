# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::version()
# # If your bioconductor version is previous to 4.0, see the section bellow
# 
# ## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
# 
# ## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# 
# devtools::install_github("aertslab/SCENIC") 
# packageVersion("SCENIC")

library(SCENIC)
library(Seurat)

so = readRDS("~/Manuscript/Data/Seurat_rds/ccbr1119_malign.rds")

Idents(so) <- so@meta.data$custom_cluster

## Get data from sce object:
exprMat <- as.matrix(so@assays$SCT@data)
cellInfo <- data.frame(cellType = so@meta.data$custom_cluster)
rownames(cellInfo) <- so@meta.data$Barcode

org <- "mgi" # or hgnc, or dmel
dbDir="~/" # RcisTarget databases location
myDatasetTitle <- "CCBR_1119_Malign_SCENIC" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]

# scenicOptions contains information on the run parameters for subsequent steps of SCENIC
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

### Co-expression network
# Generate expression matrix
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

# Calculate correlation (whether targets are positively/negatively regulated by Transcription Factors)
runCorrelation(exprMat_filtered, scenicOptions)

# Run only if data has not been log-normalized
#exprMat_filtered_log <- log2(exprMat_filtered+1) 

# procedure that aims at recovering a gene regulatory network from multifactorial expression data
# Use GRNBoost to save time to retrieve list of potential Transcription Factor Targets
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 5
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=NULL)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

library(grid)
# Returns heatmap of co-factor expression across clusters
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                        show_row_names = T, row_names_gp = gpar(fontsize = 6))

# Filters out low-confidence extended_regulons
filtered_regulons <- regulonActivity_byCellType_Scaled[
  !grepl("extended",rownames(regulonActivity_byCellType_Scaled)),]

# Trim results further
filtered_genes <- list()
for (cluster in colnames(filtered_regulons)){
  filtered_genes[[cluster]] <- names(filtered_regulons[,cluster][
    filtered_regulons[,cluster] > quantile(filtered_regulons[,cluster], 0.9)])
}

filtered_regulons2 <- filtered_regulons[unique(unlist(filtered_genes)),]
colnames(filtered_regulons2) <- c("Tum1","Tum2","Tum3")

ComplexHeatmap::Heatmap(filtered_regulons2, name="Regulon activity",
                        show_row_names = T, row_names_gp = gpar(fontsize = 10), 
                        column_order = c("Tum1","Tum2","Tum3"))