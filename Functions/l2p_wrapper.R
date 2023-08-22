l2p_wrapper <- function(wgcna_module){

  library(l2p)
  
  wcgna_mod <- read.csv("~/Manuscript/Data/Supplementary/ccbr1119_wcgna_modules.csv")
  
  module1 <- toupper(wcgna_mod$mod1[wcgna_mod$mod1 != "NA"])
  module2 <- toupper(wcgna_mod$mod2[wcgna_mod$mod2 != "NA"])
  module3 <- toupper(wcgna_mod$mod3[wcgna_mod$mod3 != "NA"])
  
  mod_func <- switch(wgcna_module,
                     "1" = l2p(module1),
                     "2" = l2p(module2),
                     "3" = l2p(module3)
  )
  
  mod_func <- mod_func[mod_func$category %in% c("GO",
                                                "REACTOME",
                                                "KEGG",
                                                "H"),]
  
  mod_func$pathway_db <- gsub(" ","_",toupper(
    paste(mod_func$category, mod_func$pathway_name, sep = " ")))
  
  mod_func <- mod_func[mod_func$fdr <= 0.05,]
  
  mod_func <- head(mod_func[order(mod_func$enrichment_score,
                                  decreasing = T),], 100)
  
  mod_func$description <- unlist(lapply(
    mod_func$pathway_id,l2pgetlongdesc))
  
  return(mod_func)
}