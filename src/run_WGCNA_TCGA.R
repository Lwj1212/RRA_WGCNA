# function
setwd("~/RRA-WGCNA/")
source("src/r-function.R")

# load RData
if(!exists("robustdegs")){
  load(paste0(base_dir, "/Step2_RRA_DEG.RData"))
}

# weighted gene co-expression network costruction
mch_test <- purrr::quietly(mergecutheight_test)(pr_name = pr_name, robustdegs = robustdegs) 
mch_test_value <- mch_test$result$max_mch

network <- network_preprocessing(pr_name = pr_name, robustdegs = robustdegs, 
                                 mch = mch_test_value, time_stamp = time_stamp)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]
MEs0 <- moduleEigengenes(network[[1]], moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    
# find key hub gene
key_hub_gene <- find_key_modulegene(pr_name = pr_name, base_dir = base_dir,
                                    network = network, MEs = MEs, 
                                    select_clinical = NULL,
                                    mm = 0.85, gs = 0.2)

save(key_hub_gene, file = paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))


# GO / KEGG enrichment analysis (ORA)
# ora_go_kegg(gs = key_hub_gene$intra_module_gene$turquoise$GS,
#             geneName = key_hub_gene$intra_module_gene$turquoise$gene,
#             module_name = "Turquoise",
#             base_dir = base_dir)
