# function
setwd("~/RRA-WGCNA/")
source("src/r-function.R")

# load RData
if(!exists("key_hub_gene")){
  load(paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))
}

# Machine learning validation
module_candidate <- list()
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg

for(m_name in names(key_hub_gene$intra_module)){
  print(m_name)
  total_keyhub <- key_hub_gene$intra_module[[m_name]]
  
  # gene selection
  selected_gene <- gene_selection(base_dir = base_dir,
                                  total_keyhub_list = total_keyhub, 
                                  over_sampling = TRUE) %>% compact()
  
  final_candidate_gene <- ml_validation(base_dir = base_dir, selected_gene = selected_gene,
                                        over_sampling=TRUE, cv = FALSE, module_name = m_name) %>%
    auc_cutoff(sg_list = ., selected_gene = selected_gene, auc_cutoff = 0.7)

  # final write
  module_candidate[[m_name]] <- tibble(GENE_NAME = final_candidate_gene$auc_cutoff) %>%
    mutate(TRAIT = paste0(final_candidate_gene$trait, collapse = ";"),
           MODULE = m_name) %>%
    arrange(GENE_NAME)
}

# module size plot
module_size <- module_candidate %>% bind_rows() %>% 
  group_by(MODULE) %>% 
  summarise(cnt = n())
colnames(module_size) <- c("module", "module_size")
gene_module_size_plot(gene_module_key_groph = module_size, 
                      title = "After Machine Learning Validation, Module size",
                      save_path = paste(base_dir, "ML_LOG", sep = "/"),
                      prefix = "Step4_Module_size")


save(module_candidate, file = paste0(base_dir, "/Step4_MODULE_SELECTED_GENE.RData"))

# Furter analysis
module_candidate_analysis <- further_analysis(mc = module_candidate, ge = geneExpression, key_hub_gene = key_hub_gene,
                        base_dir = base_dir)

final_module_candidate <- filtering_combine(pr_name = pr_name, mc = module_candidate_analysis)

write_csv(final_module_candidate, file = paste0(base_dir, "/FINAL_CANDIDATE.csv"),na = "")
  
  

