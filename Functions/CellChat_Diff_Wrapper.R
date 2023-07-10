cellchat_diff_wrapper <- function(so,
                                  time1,
                                  time2){
  
  library(CellChat)
  library(patchwork)
  
  # Subset by relevant celltype
  Idents(so) <- so@meta.data$Likely_CellType
  
  # Change and append name for cholangiocyte and malignant cells so CellChat can run differential expression
  meta <- so@meta.data %>% mutate(Likely_CellType2 = case_when(
    Likely_CellType == "Cholangiocyte" & pooled_ident == "Ctrl" ~ "Epithelial_Malignant",
    Likely_CellType == "Malignant" ~ "Epithelial_Malignant",
    TRUE ~ Likely_CellType
  ))
  
  so@meta.data$Likely_CellType <- meta$Likely_CellType2
  
  # Filter by celltype
  Idents(so) <- so@meta.data$Likely_CellType
  so <- subset(so, idents = c('Endothelial','TAM_M1','B_cells',
                              'NK_cells','CD4','CD8','Fibroblasts','TAM_M2',
                              'Dendritic',
                              'Treg','Epithelial_Malignant'))
  
  # Subset by timepoint, store in list
  Idents(so) <- so@meta.data$pooled_ident
  
  so_list <- list(Ctrl = subset(so, idents = c("Ctrl")),
                  Five = subset(so, idents = c("Week_5")),
                  Eight = subset(so, idents = c("Week_8")))
  
  cellchat_list <- list()
  
  for (time in names(so_list)){
    # Load and process data
    data.input <- as.matrix(so_list[[time]]@assays$SCT@data)
    meta = so_list[[time]]@meta.data
    rownames(meta) = meta$Barcode
    cell.use = rownames(meta)
    
    meta = data.frame(labels = meta$Likely_CellType, row.names = meta$Barcode)
    
    # Create CellChat Object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    
    # Set ligand-receptor interaction database
    CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
    
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Infer cell-cell communication at signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # Calculate aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    cellchat_list[[time]] <- cellchat
  }
  
  rm(cellchat)
  
  # Change between ctrl - five_week, ctrl - eight_week, eight_week-five_week
  object.list <- list(time1_group = cellchat_list[[time1]], 
                      time2_group = cellchat_list[[time2]])
  
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  # Differential network plot of interaction number and strengths
  netVisual_diffInteraction(cellchat, weight.scale = T)

  library(ComplexHeatmap)
  
  i = 1
  # combining all the identified signaling pathways from different datasets 
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 14, height = 21)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 14, height = 21)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 14, height = 21, color.heatmap = "GnBu")
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 14, height = 21, color.heatmap = "GnBu")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
}