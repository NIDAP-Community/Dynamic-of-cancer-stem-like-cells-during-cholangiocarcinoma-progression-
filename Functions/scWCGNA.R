
scWCGNA_wrapper <- function(seurat_obj){
  
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(scWGCNA)
  library(corrplot)
  library(harmony)
  
  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())
  seurat_obj = so
  
  # Prepare data for scWGCNA
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    #features = VariableFeatures(seurat_obj),
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "ccbr1119" # the name of the scWGCNA experiment
  )
  
  # construct metacells  in each group
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("custom_cluster", "pooled_ident"), # specify the columns in seurat_obj@meta.data to group by
    k = 5, # nearest-neighbors parameter
    ident.group = 'custom_cluster' # set the Idents of the metacell seurat object
  )
  
  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)
  
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = "Tum1", # the name of the group of interest in the group.by column
    group.by='custom_cluster' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  )
  
  # Test different soft powers:
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    setDatExpr = FALSE, # set this to FALSE since we did this above
  )
  
  # plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)
  
  # assemble with patchwork
  wrap_plots(plot_list, ncol=2)
  
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)
  
  # construct co-expression network:
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=18,
    setDatExpr=FALSE
  )
  
  # expression matrix for all the WGCNA genes:
  seurat_obj <- Seurat::ScaleData(
    seurat_obj,
    features = GetWGCNAGenes(seurat_obj)
    #vars.to.regress = c('nFeature_SCT', 'percent_mt')
  )
  
  # compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(
    seurat_obj
    #group.by.vars="pooled_ident"
  )
  
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  
  # compute intramodular connectivity:
  seurat_obj <- ModuleConnectivity(seurat_obj)
  
  # rename the modules
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = "Module"
  )
  
  # get the module assignment table:
  modules <- GetModules(seurat_obj)
  
  # grab metadata with genes / module identity
  module_meta <- seurat_obj@misc$ccbr1119$wgcna_modules
  
  # remove genes that are in the grey ~ presumably non-coexpressed group
  module_meta <- module_meta[module_meta$color != "grey",]
  
  unique(module_meta$module)
  
  #setwd("CCBR_1119_WGCNA")
  genelist <- list()
  for (module in unique(module_meta$module)){
    mod_genes <- module_meta$gene_name[module_meta$module == module]
    
    mod_heat <- DoHeatmap(object = seurat_obj, group.by = "custom_cluster", features = mod_genes, slot = "scale.data", size = 2) + ylab(module) +
      theme(axis.text.y = element_blank(), axis.title.y = element_text(size=10, face="bold"))
    
    genelist[[module]] <- mod_genes
    
    print(mod_heat)
  }
  
  return(genelist)
}
