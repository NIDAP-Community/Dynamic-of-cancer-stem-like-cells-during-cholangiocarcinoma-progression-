cellchat_indv_wrapper <- function(so,
                                  timepoint,
                                  celltypes){
  
  library(CellChat)
  library(patchwork)
  options(stringsAsFactors = FALSE)
  
  Idents(so) <- so@meta.data$pooled_ident
  
  # First subset by timepoint
  # Ctrl, Five week, or Eight week depending on timepoint
  setwd("/rstudio-files/ccbr-data/users/Jing/CCBR_1119/CellChat/Ctrl")
  so <- subset(so, idents = timepoint)
  # so <- subset(so, idents = "Five")
  # so <- subset(so, idents = "Eight")
  
  # Then subset by relevant celltype
  Idents(so) <- so@meta.data$Likely_CellType
  
  so <- subset(so, idents = celltypes)
  
  # Load and process data
  data.input <- as.matrix(so@assays$SCT@data)
  meta = so@meta.data
  rownames(meta) = meta$Barcode
  cell.use = rownames(meta)
  
  meta = data.frame(labels = meta$Likely_CellType, row.names = meta$Barcode)
  
  # Create CellChat Object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  ## Append new metadata information to CellChat
  # cellchat <- addMeta(cellchat, meta = meta)
  # cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  # levels(cellchat@idents) # show factor levels of the cell labels
  # groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
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
  
  # project gene expression data onto Protein-Protein Interaction (PPI) network (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  #   could take into account the possibility of signal dropoff 
  #cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Infer cell-cell communication at signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  #### Results ####
  
  # Visualize aggregated cell-cell communication network
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  
  # Creating a plot
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  
  # Creating a plot
  net_view_weight <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat, label.size = 6, font.size = 20)
  
  plot(gg1)
}
