# function
setwd("~/RRA-WGCNA/")
source("src/r-function.R")

# ======= RUN DEA & RRA ======= 
gse_list <- c("GSE76427", "GSE54236", "GSE39791","GSE14811", "GSE87630", "GSE136247", "GSE64041", "GSE45114")
multiple_limma <- GSE_manual(gse_list = gse_list)
save(multiple_limma, file = paste0(base_dir, "/Step2_RAW.RData"))

# rra & robust DEGs
if(!exists("multiple_limma"))
  load(paste0(base_dir, "/Step2_RAW.RData"))

robustdegs <- rra_analysis(m_list = multiple_limma, logfc = 0, fdr = 0.05, save_path = base_dir)
save(robustdegs, file = paste0(base_dir, "/Step2_RRA_DEG.RData"))
