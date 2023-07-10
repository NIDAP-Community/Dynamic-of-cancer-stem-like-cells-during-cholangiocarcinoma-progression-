
survival_wrapper <- function(wcgna_mod){
  
  library(matrixStats)
  library(survival)
  library(survminer)
  library(tidyverse)

  df <- read.csv("Supplementary/Normalized_Counts_icca.csv")
  
  # From wcgna_wrapper function
  mod_genes <- read.csv("Supplementary/ccbr1119_wcgna_modules.csv")
  
  # Grab genes, remove NA
  mod_genes <- mod_genes[,wcgna_mod]
  mod_genes <- mod_genes[mod_genes != "NA"]
  
  # Convert to Uppercase for humans (approximately)
  mod_genes <- toupper(mod_genes)
  
  # Use genes found only in clinical dataset
  mod_genes <- mod_genes[mod_genes %in% df$Gene]
  
  #evaluate geometric or arithmentic mean
  #gl <- trimws(unlist(strsplit(c("GZMA,PRF1"), ",")), which=c("both")) # unpack the gene list provided by the user and remove white spaces 
  gl <- mod_genes   
  ind_gn <- match(gl, df$Gene)
  dfgn <- df[ind_gn,2:dim(df)[2]]
  
  dx <- colSums(dfgn)/2 #Calculate geometric mean with input log values 
  dg <- df[-ind_gn,]
  dy <- df[-ind_gn,2:dim(df)[2]]
  
  # evaluate correlation
  df_out <- as.data.frame(as.matrix(NA, nrow=dim(dg)[1], ncol=5))
  
  for (i in 1:dim(dg)[1]){
    cc <- cor.test(dx, as.numeric(dy[i,]), method="pearson")
    df_out[i,1] <- dg$Gene[i]
    df_out[i,2] <- round(cc$estimate,3)
    df_out[i,3] <- cc$p.value
    df_out[i,4] <- round(cc$conf.int[1],3)
    df_out[i,5] <- round(cc$conf.int[2],3)
  }
  colnames(df_out) <- c("Gene", "CorrCoef", "Pval","LowerCI", "UpperCI")
  
  ## Survival Plots ##
  df <- read.csv("Supplementary/Normalized_Counts_icca.csv")
  dhc <- read.csv("Supplementary/AndersenJBetal_Gastroenterology 2012_CCA104_ClinINFO.xlsb_Sheet5.csv")
  dg <- df_out
  
  # select only tumor samples
  ind_tp_grp <- grep("iCCA",dhc$Cancer)
  pt_samples <- dhc$SampleA[ind_tp_grp]
  dfpt <- df[,c("Gene",pt_samples)]
  
  # select the expression data for the gene list of choice
  gn <- dg %>% filter(CorrCoef < 0, Pval < 0.05) %>% pull(Gene) 
  ind <- match(gn, dfpt$Gene)  
  dfg <- as.matrix(dfpt[ind,2:dim(dfpt)[2]])
  dfg <- t(scale(t(dfg)))
  
  # classify samples above and below upper and lower quartiles of aggregated gene expression differences 
  md_thr <- rowMedians(dfg) # median across samples for the list of genes of choice
  diff_thr <- dfg-md_thr # difference between the expression and median threshold
  above_threshold <- colMedians(diff_thr)>0 # samples with difference above threshold
  below_threshold <- colMedians(diff_thr)<0 # samples with difference below threshold 
  label_above_thr <- "above median"
  label_below_thr <- "below median"
  
  samples_above_threshold <- colnames(dfpt[,-1])[above_threshold]
  samples_below_threshold <- colnames(dfpt[,-1])[below_threshold]
  
  ind_patients_above_threshold_hclinial <- match(samples_above_threshold,dhc$SampleA)
  ind_patients_below_threshold_hclinial <- match(samples_below_threshold,dhc$SampleA)
  
  # prepare the vectors necessary for survival analysis, if dhc vital is dead, use days to death, if not dead, use days to followup
  surv_days_all <- dhc$Days_elapsed
  ind_vec <- c(ind_patients_above_threshold_hclinial, ind_patients_below_threshold_hclinial)
  
  surv_days <- surv_days_all[ind_vec]
  
  cs <- dhc$"LD"[ind_vec]
  csn <- as.numeric(cs=="1") # if the vital_status is "Dead" attribute 1; else, attribute 0
  surv_censor <- csn
  cluster <- c(rep(0,length(ind_patients_above_threshold_hclinial)), rep(1, length(ind_patients_below_threshold_hclinial)))
  
  # prepare the dataframe necessary for survival analysis
  ds <- as.data.frame(matrix(NA, ncol=3, nrow=length(cluster)))
  ds[,1] <- surv_days
  ds[,2] <- surv_censor
  ds[,3] <- cluster
  colnames(ds) <- c("surv_days", "surv_censor", "cluster")
  
  # evaluate survival curves
  fit1 <- survfit(Surv(as.numeric(surv_days),surv_censor) ~ cluster, data = ds) # single variable
  
  # explicit measures of survival fit
  pval <- surv_pvalue(fit1,ds)
  sfit <- summary(fit1)$table
  print(pval)
  print(sfit)
  
  p <- survminer::ggsurvplot(fit1, #plot just the KM curves 
                             data = ds,
                             pval = TRUE, # p-value of the Log-Rank test comparing the groups 
                             conf.int = TRUE,  # 95% confidence interval
                             risk.table = TRUE, # display the risk table below the graph 
                             risk.table.col = "strata",  # control what is displayed in the table
                             title = unique(dhc$"project_id"),
                             surv.median.line = "hv", # display horizontal and vertical line for median survival
                             ggtheme = theme_gray(), # have a background theme
                             palette = c("red","black"),  # control the colors for the KM curves and associated data
                             ncensor.plot = TRUE, # add censor plot below the risk table to display the censor data
                             legend.labs = c(label_above_thr, label_below_thr)
  ) 
  print(p)           
  
  df_out <- as.data.frame(sfit)
  df_out <- cbind(as.data.frame(rownames(df_out)), df_out)
  colnames(df_out)[c(1,6:7)] <- c("cluster","rmean", "se_rmean")
  
  select_output <- 1
  if (select_output==1){
    df_out_r <- pval
  } else if (select_output==2){
    df_out_r <- df_out
  } else if (select_output==3){
    df_out_r <- ds
  } 
  
  return(df_out_r)
}
