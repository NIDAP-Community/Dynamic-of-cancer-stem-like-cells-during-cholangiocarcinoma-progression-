library(Seurat)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(ggplot2)

localFilePaths <- list.files("~/Manuscript/Data/h5",full.names = TRUE)

obj.list <- lapply(localFilePaths, function(x) {
  return(Read10X_h5(x, use.names=TRUE))})

names(obj.list) <- c("Ctrl1","Ctrl2","Ctrl3","Eight_week1","Eight_week2",
                     "Eight_week3","Five_week1","Five_week2","Five_week3",
                     "Ctrl_Five_week1","Ctrl_Five_week2","Ctrl_Five_week3",
                     "Ctrl_Eight_week1","Ctrl_Eight_week2","Ctrl_Eight_week3")
obj.list <- obj.list[sort(names(obj.list))]

mincells = 3
mingenes = 200
organism = "Mouse"

mitoch = "^mt-"
cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
cc.genes$s.genes = str_to_title(cc.genes$s.genes)

seurat_object <- function(i) {
  
  so.nf <- so.orig.nf[[i]]
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  if ("Protein" %in% names(so.nf)){
    so.nf <- NormalizeData(so.nf, assay = "Protein", normalization.method = "CLR")
  }
  
  if ("HTO" %in% names(so.nf)){
    so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
  }
  
  so <- so.nf
  so.origcount = dim(so.nf)[2]
  
  #Start with filtering here:
  maxgenes = 2500
  complexity = 0.8
  minUMI = 500
  MAD_gene <- TRUE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 10
  
  MAD_mitoch <- TRUE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  if (MAD_gene == TRUE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    cat(paste0("Complexity Filter =",complexity,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
  }
  
  df.m <- melt(so@meta.data)
  df.m$filt <- "filt"
  df.m$filt <- as.factor(df.m$filt)
  df2.m <- melt(so.nf@meta.data)
  df2.m$filt <- "raw"
  df2.m$filt <- as.factor(df2.m$filt)
  
  v <- unique(df.m$variable)
  
  so2.list <- list(so,so.nf)
  
  return(so2.list)
}

#Create Seurat Object from original H5, splitting to multiple ones as necessary.  For splitting H5's only RNA slot is expected and supported.
so.orig.nf <- list()
for(i in seq_along(names(obj.list))){
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
  } else {
    k = names(obj.list[[i]])
    for(j in 1:length(k)){
      if(names(obj.list[[i]][j]) == "Gene Expression"){
        so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]][k][[j]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
      } else if(names(obj.list[[i]][j]) == "Antibody Capture"){
        Protein = rownames(obj.list[[i]][k][[j]])[!grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
        HTO = rownames(obj.list[[i]][k][[j]])[grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
        so.orig.nf[[i]][["Protein"]] <- CreateAssayObject(obj.list[[i]][k][[j]][Protein, colnames(so.orig.nf[[i]])])
        if(length(HTO)>0){
          so.orig.nf[[i]][['HTO']] <- CreateAssayObject(counts = obj.list[[i]][k][[j]][HTO, colnames(so.orig.nf[[i]])])}
      } else{
        print(paste(names(obj.list[[i]][j]),"found, not stored"))
      }
    }
  }
  names(so.orig.nf)[[i]] <- names(obj.list)[[i]]
}

so.list <- lapply(seq_along(so.orig.nf), seurat_object)

SO <- lapply(so.list,function(x) x[[1]])
names(SO) <- unlist(lapply(so.list, function(x) as.character(Seurat::Idents(x[[1]])[1])))

rm(so.list)
rm(so.f.list)
rm(so.orig.nf)

### PCA and Normalization ###
vars_to_regress <- c()
npcs = 30

# Linearly scale data without regressing anything.
scale_so <- function(so){
  so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst")
  all.genes <- rownames(so)
  so <- ScaleData(so,features=all.genes)
  return(so)
}

# Make PCA with SCTransform() and optional ScaleData, and do so with
# both regression (if user requests) and on all genes.
pca <- function(so) {
  # If user sets Linear Scaling toggle TRUE, also run ScaleData().
  # Use case: user has legacy project from Seurat 2 and wants to keep
  # methods consistent with pre-SCT Seurat.
  # Run SCTransform().
  if(is.null(vars_to_regress)){
    so <- so
  }
  else { 
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
  }
  # Make PCA using last transform run, which will always be that from
  # SCTransform().
  so <- RunPCA(object = so, npcs = npcs)
  slot(so,"commands") <- list()
  return(so)
}

# Do transformation with and without regression using SCTransform()
# and ScaleData().
so_scale <- lapply(SO, scale_so) 

### Combine and Renormalize ###
SO <- lapply(so_scale, pca)
rm(so_scale)
rm(obj.list)

## Add commentary on how this toggle works with add.only.var.genes.
conserve_memory <- TRUE

#initialize Citeseq functionality as false, 
#later the template will check for a Protein assay and run if it finds it
dat = vector()
integratedata = FALSE

for(i in 2:length(SO)){dat=c(dat,SO[[i]])}
SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
allgenes <- rownames(SO_merge)

rm(dat)
rm(SO)
gc()

Do_SCTransform = TRUE
vars_to_regress = c()

SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, conserve.memory = conserve_memory, return.only.var.genes = FALSE)

SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 5000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)

for (i in seq(0.2,1.2,0.2)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
}

