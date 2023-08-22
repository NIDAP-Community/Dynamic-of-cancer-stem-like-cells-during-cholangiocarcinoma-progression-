library(DescTools)
library(ggplot2)
library(ggpubr)
library(ggsignif)

# Load the Seurat object
so <- readRDS("Manuscript/Data/Seurat_rds/ccbr1119_malign.rds")

# Extract expression data
malign_expr <- data.frame(cell = so$Barcode, 
                          cluster = so$custom_cluster, 
                          Tm4sf1 = so$SCT@data["Tm4sf1",], 
                          Mki67 = so$SCT@data["Mki67",], 
                          Prom1 = so$SCT@data["Prom1",], 
                          Sox9 = so$SCT@data["Sox9",])

dunn_res <- list()

for(gene in c("Tm4sf1","Mki67","Prom1","Sox9")){
  kruskal_result <- kruskal.test(malign_expr[[gene]] ~ cluster, 
                                 data = malign_expr)
  
  print(kruskal_result)
  
  if(kruskal_result$p.value <= 0.05){
    # Post-Hoc
    dunn_out <- DunnTest(malign_expr[[gene]], malign_expr$cluster, 
                                 method = "BH")
    
    dunn_res[[gene]] <- dunn_out$pmat
  } else{
    print("non-significant Kruskal-Wallis p-value")
  }
}
