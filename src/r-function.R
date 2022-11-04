# library ----
suppressMessages({
  # Bioconductor
  library(GEOquery) # bio
  library(limma) # bio
  library(DESeq2) # bio
  library(impute) #bio
  library(TCGAbiolinks) # bio
  library(WGCNA) # bio
  library(clusterProfiler) # bio
  library(enrichplot) # bio
  library(org.Hs.eg.db) # bio
  library(BiocParallel) # bio
  library(sva) # bio
  library(MultiAssayExperiment) # bio
  
  # Github
  library(ggVennDiagram) #github 
  library(CorLevelPlot) # github
  library(ELMER) # bio
  
  # CRAN
  library(glue)
  library(RobustRankAggreg) # cran
  library(pheatmap) # cran
  library(survival) # cran
  library(survminer) # cran
  library(rbioapi) # cran
  library(parallel) # base
  library(httr) # base
  library(RMariaDB) # cran
  library(glmnet)
  library(reticulate)
  library(tidyverse) # cran
  
  source_python("src/py-function.py")
  use_python("/usr/bin/python3")  
})


# global variable ----
time_stamp <- Sys.time() %>% str_split(pattern = " ") %>% 
  unlist() %>% .[1]
pr_name <- readline('enter cancer type : ')
base_dir <- paste("WGCNA_PIPELINE_RESULT", pr_name, time_stamp,sep = "/")
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# function ----
retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=5, sleep=0) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
      flog.fatal(msg)
      stop(msg)
    } else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors, 
                    capture.output(str(retval)))
      flog.error(msg)
      warning(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(retval)
}

#' Function that returns numeric values with 2 decimal numbers.
#' 
#' @param x input numeric value with N decimal numbers.
#' @return a numeric value with 2 decimal numbers.
#' @examples
#' aaa <- dec_two(8.31232)
dec_two <- function(x) {
  return (format(round(x, 2), nsmall = 2));
}

ora_go_kegg <- function(geneName, module_name, base_dir){
  # Create directory
  save_path <- paste0(base_dir, "/ANALYSIS/GO_KEGG_ORA/", module_name, "/")
  dir.create(path = save_path, showWarnings = FALSE, recursive = TRUE)
  
  # Symbol to Entrez
  geneList <- bitr(geneName, 
                   fromType="ALIAS", 
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, 
                   drop = FALSE) %>% 
    distinct(ALIAS, .keep_all = TRUE) %>% 
    pull(2)
  
  # GO ORA - BP
  ego_ora_bp <- enrichGO(geneList, 
                         OrgDb=org.Hs.eg.db, 
                         ont="BP", 
                         pAdjustMethod = "bonferroni",
                         pvalueCutoff=0.9,
                         qvalueCutoff=0.9,
                         readable=T) %>%
    mutate(., Description = paste0("(", ID ,") ", Description))
  
  p_ora_bp <- barplot(ego_ora_bp, title = paste0("GO:Biological Process(BP)-", module_name), 
                      showCategory = 20) + theme_linedraw(base_size = 30)
  
  # GO ORA CC
  ego_ora_cc <- enrichGO(geneList, 
                         OrgDb=org.Hs.eg.db, 
                         ont="CC", 
                         pAdjustMethod = "bonferroni",
                         pvalueCutoff=0.9,
                         qvalueCutoff=0.9,
                         readable=T) %>%
    mutate(., Description = paste0("(", ID ,") ", Description))
  p_ora_cc <- barplot(ego_ora_cc, title = paste0("GO:Cellular Component(CC)-", module_name), 
                      showCategory = 20) + theme_linedraw(base_size = 30)
  
  # GO ORA MF
  ego_ora_mf <- enrichGO(geneList, 
                         OrgDb=org.Hs.eg.db, 
                         ont="MF", 
                         pAdjustMethod = "bonferroni",
                         pvalueCutoff=0.9,
                         qvalueCutoff=0.9,
                         readable=T) %>%
    mutate(., Description = paste0("(", ID ,") ", Description))
  
  p_ora_mf <- barplot(ego_ora_mf, title = paste0("GO:Molecular Function(MF)-", module_name), 
                      showCategory = 20) + theme_linedraw(base_size = 30)
  
  
  
  # KEGG
  kk <- enrichKEGG(gene = geneList, organism = 'hsa', 
                   pvalueCutoff=0.9,
                   qvalueCutoff=0.9) 
  p_kk <- barplot(kk, title = paste0("KEGG -", module_name), 
                  showCategory = 20) + theme_linedraw(base_size = 30)
  
  # table save
  ego_ora_bp@result %>% write_csv(file = paste0(save_path, module_name, "-GO-BP.csv"))
  ego_ora_cc@result %>% write_csv(file = paste0(save_path, module_name, "-GO-CC.csv"))
  ego_ora_mf@result %>% write_csv(file = paste0(save_path, module_name, "-GO-MF.csv"))
  kk@result %>% write_csv(file = paste0(save_path, module_name, "-KEGG.csv"))
  
  # graph save
  ggsave(plot = p_ora_bp, filename = paste0(save_path, module_name,"_GO-BP.png"), width = 40, height = 25, dpi = 300)
  ggsave(plot = p_ora_cc, filename = paste0(save_path, module_name,"_GO-CC.png"), width = 40, height = 25, dpi = 300)
  ggsave(plot = p_ora_mf, filename = paste0(save_path, module_name,"_GO-MF.png"), width = 40, height = 25, dpi = 300)
  ggsave(plot = p_kk, filename = paste0(save_path, module_name,"_KEGG.png"), width = 40, height = 25, dpi = 300)
}

# STRING analysis
string_analysis <- function(mc, base_dir, score = 400){
  log_save <- paste(base_dir, "ANALYSIS/STRING", sep = "/")
  dir.create(paste(base_dir, "ANALYSIS/STRING", sep = "/"), showWarnings = FALSE, recursive = TRUE)
  
  
  suppressMessages({
    string_filtering <- lapply(X = names(mc), FUN = function(mc_name){
      g <- mc[[mc_name]]
      proteins_mapped <- rba_string_map_ids(ids = g$GENE_NAME, species = 9606)
      int_net <- rba_string_interactions_network(ids = proteins_mapped$stringId,
                                                 species = 9606,
                                                 required_score = score) # minimum confidence
      
      rba_string_network_image(ids = proteins_mapped$stringId,
                               image_format = "highres_image",
                               species = 9606,
                               save_image = paste0(getwd(),
                                                   "/", base_dir, "/ANALYSIS/STRING/", mc_name, "-PPI.png"),
                               required_score = score,
                               network_flavor = "confidence", 
                               hide_disconnected_nodes = TRUE)
      
      ee <- FALSE
      tryCatch(
        expr = {
          string_filtered_g <- c(int_net$preferredName_A, int_net$preferredName_B) %>% 
            unique() %>% 
            tibble(GENE_NAME = .) %>% 
            inner_join(x = g, y = ., by = "GENE_NAME")
        }, 
        error = {
          function(e){
            ee <<- TRUE
          }})
      
      if(ee)
        return(NULL)
      else
        return(string_filtered_g)
    })  
  })
  
  names(string_filtering) <- names(mc)
  
  return(string_filtering)
}

# Survival analysis
survival_analysis <- function(base_dir, geneExpression, mc){
  log_save <- paste(base_dir, "ANALYSIS/OS", sep = "/")
  dir.create(paste(base_dir, "ANALYSIS/OS", sep = "/"), showWarnings = FALSE, recursive = TRUE)
  
  survival_trait <- read_delim(paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2F",pr_name,"_survival.txt"),
                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
  
  suppressMessages({
    survival_filtering <- lapply(X = names(mc), function(mc_name){
      g <- mc[[mc_name]] %>% pull(1)
      suv_exp <- geneExpression %>% rownames_to_column(var = "sample") %>% 
        dplyr::select(sample, all_of(g)) %>% 
        filter(str_ends(sample, "-01"))
      
      deg_list <- suv_exp %>% colnames() %>% .[-1]
      suv_exp_trait <- inner_join(x = suv_exp, y = survival_trait, by = "sample")
      
      os_list <- list()
      dir.create(paste(log_save, mc_name, sep = "/"), showWarnings = FALSE, recursive = TRUE)
      for(gene_name in deg_list){
        # print(gene_name)
        exp_median <<- suv_exp_trait %>% 
          mutate(group = ifelse(.[[gene_name]] > median(.[[gene_name]]), "High-exp", "Low-exp")) %>% 
          mutate(group = factor(group))
        
        exp_median <<- exp_median %>%
          mutate(group = factor(group, levels = c("Low-exp", "High-exp")))
        
        sfit <- survfit(Surv(OS.time, OS) ~ group, data = exp_median)
        sdf <- survdiff(Surv(OS.time, OS) ~ group, data = exp_median)
        coxph_ <- summary(coxph(Surv(OS.time, OS) ~ group, data = exp_median))
        
        hazard_ratio <- coxph_$coefficients[[2]]
        # p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        
        ggsurvplot(sfit, title = paste0(gene_name, "-Expression High-Low(Median)"), pval = TRUE)
        ggsave(filename = paste0(log_save, "/", mc_name, "/", gene_name, "-OS.png"), device = "png")
        
        os_list[[gene_name]] <- tibble(GENE_NAME = gene_name, 
                                       HR = hazard_ratio,
                                       HR.P = coxph_$coefficients[[5]]) 
      }
      
      os_list %>% bind_rows() %>% 
        inner_join(x =  mc[[mc_name]], y = ., by = "GENE_NAME") %>% 
        return()
    })
  })
  
  names(survival_filtering) <- names(mc)
  return(survival_filtering)
}

# DNA methylation analysis
methylation_analysis <- function(pr_name, method_="both", base_dir){
  
  log_save <- paste(base_dir, "ANALYSIS/ELMER", sep = "/")
  dir.create(paste(base_dir, "ANALYSIS/ELMER", sep = "/"), showWarnings = FALSE, recursive = TRUE)
  
  if(file.exists(paste0(log_save, "/", pr_name, "_methylation_filtering.RData"))){
    load(paste0(log_save, "/", pr_name, "_methylation_filtering.RData"))
    
  } else {
    
    if(method_ == "both")
      method_ <- c("hypo", "hyper")
    
    if(!file.exists(paste0(getwd(), "/Data/", pr_name, "/", pr_name, "_RNA_hg19.rda"))){
      getTCGA(disease = pr_name, genome = "hg19")  
      
    } else {
      
      load(paste0(getwd(), "/Data/", pr_name, "/", pr_name, "_RNA_hg19.rda"))
      load(paste0(getwd(), "/Data/", pr_name, "/", pr_name, "_meth_hg19.rda"))
      geneExp <- assay(rna)
      Meth <- assay(met)
      rm(rna, met);gc()
    }
    
    distal.probes <- get.feature.probe(
      genome = "hg19", 
      met.platform = "450K", 
      rm.chr = paste0("chr",c("Y")) # Y는 male specific
      # rm.chr = paste0('chr', c(1:6, 10:22, "X","Y"))
    )
    
    mae <- createMAE(
      exp = geneExp, 
      met = Meth,
      save = TRUE,
      linearize.exp = TRUE,
      save.filename = paste0(log_save, "/", pr_name,"_mae.rda"),
      filter.probes = distal.probes,
      met.platform = "450K",
      genome = "hg19",
      TCGA = TRUE
    )
    
    sig.diff_list <- list()
    nearGenes_list <- list()
    methylation_filtering <- list()
    enriched.motif_list <- list()
    TF_list <- list()
    
    for(method in method_){
      sig.diff_list[[method]] <- get.diff.meth(data = mae, 
                                               group.col = "definition",
                                               group1 =  "Primary solid Tumor",
                                               group2 = "Solid Tissue Normal",
                                               minSubgroupFrac = 0.2, # if supervised mode set to 1
                                               sig.dif = 0.3,
                                               diff.dir = method, # Search for hypomethylated probes in group 1
                                               cores = 20, 
                                               dir.out = paste0(getwd(), "/", log_save), 
                                               pvalue = 0.05)
      
      nearGenes_list[[method]] <- GetNearGenes(data = mae, 
                                               probes = sig.diff_list[[method]]$probe, 
                                               numFlankingGenes = 30)
      
      methylation_filtering[[method]] <- get.pair(data = mae,
                                                  group.col = "definition",
                                                  group1 =  "Primary solid Tumor",
                                                  group2 = "Solid Tissue Normal",
                                                  nearGenes = nearGenes_list[[method]],
                                                  mode = "unsupervised",
                                                  permu.dir = paste0(getwd(), "/", log_save, "/permu"),
                                                  permu.size = 10000, # Please set to 100000 to get significant results
                                                  raw.pvalue = 0.05,   
                                                  Pe = 0.001, # Please set to 0.001 to get significant results
                                                  filter.probes = TRUE, # See preAssociationProbeFiltering function
                                                  filter.percentage = 0.05,
                                                  filter.portion = 0.3,
                                                  dir.out = log_save,
                                                  cores = 20,
                                                  label = method)
      
      # # Mehtylation enrichiment analsis
      # enriched.motif_list[[method]] <- get.enriched.motif(data = mae,
      #                                                     probes = pair_list[[method]]$Probe, 
      #                                                     dir.out = paste0(getwd(), "/", log_save), 
      #                                                     label = method,
      #                                                     min.incidence = 10,
      #                                                     lower.OR = 1.1)
      # 
      # # TF analyis
      # TF_list[[method]] <- get.TFs(data = mae, 
      #                              group.col = "definition",
      #                              group1 =  "Primary solid Tumor",
      #                              group2 = "Solid Tissue Normal",
      #                              mode = "unsupervised",
      #                              enriched.motif = enriched.motif_list[[method]],
      #                              dir.out = paste0(getwd(), "/", log_save), 
      #                              cores = 1, 
      #                              label = method)
    }
    
    save(methylation_filtering, file = paste0(log_save, "/", pr_name, "_methylation_filtering.RData"))
  }
  
  return(methylation_filtering)
}



# mapping function ====
symbol_mapping <- function(ge, col_name, platform_ann_df){
  thisGene <- rownames(ge) %>%
    lapply(FUN = function(value){
      platform_ann_df %>% 
        dplyr::filter(ID == value) %>% 
        pull(col_name)
    }) %>% 
    do.call(c, .)
  
  if(col_name != "gene_assignment"){
    
    return(thisGene)
    
  } else { # gene assignment
    
    SIZE_SPLIT_STRING <- 2
    GENE_SYMBOL_INDEX <- 2
    
    thisGene %>% 
      lapply(X = ., FUN = function(value){
        split_string <- strsplit(value, "//")
        if (length(split_string[[1]]) >= SIZE_SPLIT_STRING) {
          thisGeneSymbol_temp <- NULL
          thisGeneSymbol  <- NULL
          thisGeneSymbol_temp <- split_string[[1]][GENE_SYMBOL_INDEX]
          thisGeneSymbol <- gsub(" ", "", thisGeneSymbol_temp, fixed = TRUE)
          thisGeneSymbol %>% return()
        } else {
          value %>% return()
        }
      }) %>% do.call(c, .) %>% 
      return()
  }
}
biodbnet_db2db <- function(id){
  base_url <- "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/"
  json_url <- paste0(base_url, "biodbnetRestApi.json?")
  mapping_result <- list()
  id_chunk <- split(id, ceiling(seq_along(id)/1000))
  
  for(index in 1:length(id_chunk)){
    parameters <- list(method="db2db", 
                       inputValues=id_chunk[[index]],
                       input="genbanknucleotideaccession",
                       outputs=c("genesymbol"),
                       taxonId="9606",
                       format="row")
    
    mapping_result[[index]] <- rcurl_request(json_url, parameters)
  }
  
  mapping_result <- mapping_result %>% 
    bind_rows() %>% 
    filter(`Gene Symbol` != "-")
  
  return(mapping_result)
}


# GEO download function  ====
GSE_manual <- function(gse_list, pheno_edit = TRUE){
  # now_date <- Sys.time() %>% str_split(pattern = " ") %>% unlist() %>% .[1]
  dir.create(paste0(base_dir, "/DEA"), recursive = TRUE, showWarnings = FALSE)
  multiple_limma <- list()
  
  tryCatch(
    expr = {
      
      for(gse_name in gse_list){
        print(gse_name)
        if(file.exists(paste0(base_dir, "/DEA/", gse_name, "_limma.txt"))){
          multiple_limma[[gse_name]] <- read_delim(paste0(base_dir, "/DEA/", gse_name, "_limma.txt"), delim = "\t",
                                                   show_col_types = FALSE)
          next
        } else {
          # file load
          if(file.exists(paste0("GSE/", gse_name, ".RData"))){
            load(file = paste0("GSE/", gse_name, ".RData"))
          } else {
            dir.create("GSE/", showWarnings = FALSE, recursive = TRUE)
            gse_data <- retry(expr = getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                                              retrieveGeneSymbols = TRUE, 
                                                              verbose = TRUE),
                              maxErrors = 20
            )
            save(gse_data, file = paste0("GSE/", gse_name, ".RData"))
          }
          
          geneExpression <- gse_data$gene_expression
          boxplot(geneExpression[1:10, 1:10])
          
          # run log2-transformation
          log2_trans <- readline('Run log2-transformation? [y]/[n] : ')
          if(tolower(log2_trans) == "yes" | tolower(log2_trans) == "y"){
            log2_trans <- TRUE
          } else {
            log2_trans <- FALSE
          }
          
          if(log2_trans){
            geneExpression <- log2(geneExpression)
            boxplot(geneExpression[1:10, 1:10])
          }
          
          # phenotype selection
          if(file.exists(paste0("GSE/", gse_name, "_pheno.txt"))){
            pheno <- read.delim(file = paste0("GSE/", gse_name, "_pheno.txt"), sep = "\t")
          } else {
            pheno <- gse_data$pheno@data %>% mutate_all(tolower)  
            write.table(pheno, file = paste0("GSE/", gse_name, "_pheno.txt"), sep = "\t", row.names = TRUE)
          }
          
          pheno_check <- lapply(X = names(pheno), FUN = function(col_name){
            df <- pheno[col_name]
            df_el <- df %>% pull(1) %>% unique()
            
            if(length(df_el) < 2 | length(df_el) > 12){
              return(NULL) 
            } else {
              df_el <- c(length(df_el), df_el)
              df_el <- paste0(df_el, collapse = " / ")
              tibble(col_name = col_name, factor = df_el) %>% 
                return()
            }
          }) %>% compact() %>% bind_rows()
          View(pheno_check)
          
          
          
          # sample selection
          selected_pheno <- readline('enter phenotype : ')
          pheno[[selected_pheno]] <- str_replace_all(pheno[[selected_pheno]], pattern = " ", replacement = "_")
          
          # re pheno_check
          pheno_check <- lapply(X = names(pheno), FUN = function(col_name){
            df <- pheno[col_name]
            df_el <- df %>% pull(1) %>% unique()
            
            if(length(df_el) < 2 | length(df_el) > 12){
              return(NULL) 
            } else {
              df_el <- c(length(df_el), df_el)
              df_el <- paste0(df_el, collapse = " / ")
              tibble(col_name = col_name, factor = df_el) %>% 
                return()
            }
          }) %>% compact() %>% bind_rows()
          
          pheno_check %>% 
            dplyr::filter(col_name == selected_pheno) %>% 
            dplyr::pull(factor) %>% 
            print()
          
          # phenotype re-level
          relevel_check <- readline('Run re-level factor? [y]/[n] : ')
          if(tolower(relevel_check) == "yes" | tolower(relevel_check) == "y"){
            relevel_check <- TRUE
          } else {
            relevel_check <- FALSE
          }
          
          if(relevel_check){
            NT_v <- readline("enter NT : ") %>% 
              str_split(pattern = " ") %>% unlist()
            TP_v <- readline("enter TP : ") %>% 
              str_split(pattern = " ") %>% unlist()
            
            pheno[[selected_pheno]] <- forcats::fct_collapse(pheno[[selected_pheno]], 
                                                             NT = NT_v, TP = TP_v, other_level = 'NA') 
            pheno <- pheno %>% dplyr::filter(!!as.symbol(selected_pheno) == "NT" | !!as.symbol(selected_pheno) == "TP")
            pheno[[selected_pheno]] <- pheno[[selected_pheno]] %>% droplevels()
            
            pheno %>% group_by(!!as.symbol(selected_pheno)) %>% summarise(cnt = n()) %>% 
              print()
          }
          
          # selected pheno
          grp <- pheno %>% 
            dplyr::select(starts_with(selected_pheno)) %>% 
            pull(1) %>% 
            # filter(str_detect(`source_name_ch1`, "HCC")) %>%
            lapply(X = ., FUN = tolower) %>% 
            unlist() %>% 
            as.factor()
          # print(grp %>% unique())
          design <- model.matrix(~0 + grp)
          View(design)
          
          nt_tp_order <- readline("enter NT-TP order : ") %>% 
            str_split(pattern = " ") %>% unlist()
          colnames(design) <- nt_tp_order # 데이터에 맞추어 manual로 설정해야 됨
          
          # Limma
          
          geneExpression_filter <- geneExpression[, pheno %>% rownames()] # rowname
          limma_result <- run_limma(ge = geneExpression_filter, de = design)
          write_delim(x = limma_result, path = paste0(base_dir, "/DEA/", gse_name, "_limma.txt"), delim = "\t")
          multiple_limma[[gse_name]] <- limma_result
        }
      }
    },
    error = function(e) {
      print(e)
      return(multiple_limma)
    }
  )
  
  
  return(multiple_limma)
}

#' Function that reads in a URL to check and verifies if it exists (function taken from https://stackoverflow.com/a/12195574 )
#' 
#' @param url the URL of a webpage
#' @return the output of a webpage verification check
#' @examples y <- readUrl("http://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html")
readUrl <- function(url) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      
      readLines(con=url, warn=FALSE) 
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return("EMPTY_STRING")
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message(cond)
      # Choose a return value in case of warning
      return("EMPTY_STRING")
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message(paste("Processed URL:", url))
    }
  )    
  return(out)
}
rcurl_request <- function(service_url, parameters) {
  # Collapse all parameters into one string
  all_parameters <- paste(
    sapply(names(parameters), 
           FUN=function(param_name, parameters) {
             paste(param_name, paste(parameters[[param_name]], collapse=','), collapse='', sep='=')
           }, 
           parameters),
    collapse="&")
  
  # Paste base URL and parameters
  requested_url <- paste0(service_url, all_parameters)
  
  # Encode URL (in case there would be any space character for instance)
  # requested_url <- URLencode(requested_url)
  
  # Start request to service
  response <- GET(requested_url)
  
  raise <- content(response, as="text")
  #parse JSON
  new <- fromJSON(raise)
  
  return(new)
}
#' Function that reads in the GEO code of a dataset, and returns the gene expression dataframe.
#' 
#' @param datasetGeoCode the GEO code of a dataset.
#' @param retrieveGeneSymbols a boolean flag stating if the function should retrieve the gene symbols or not.
#' @param verbose a boolean flag stating if helping messages should be printed or not
#' @return a gene expression dataset.
#' @examples
#' geneExpressionDF1 <- getGeneExpressionFromGEO("GSE3268", FALSE, FALSE)
getGeneExpressionFromGEO <- function(datasetGeoCode, retrieveGeneSymbols, verbose = FALSE) {
  
  GSE_code <- datasetGeoCode
  
  
  # r <- NULL
  # attempt <- 1
  # while( is.null(r) && attempt <= 3 ) {
  #   attempt <- attempt + 1
  #   try(
  #     r <- some_function_that_may_fail()
  #   )
  # } 
  
  # check   URL
  checked_html_text <- "EMPTY_STRING"
  checked_html_text <- retry(expr =  xml2::read_html("https://ftp.ncbi.nlm.nih.gov/geo/series/"),
                             maxErrors=20, 
                             sleep=2)
  
  checked_html_text_url <- "EMPTY_STRING"
  url_to_check <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", datasetGeoCode)
  GSE_code_for_url <- GSE_code
  GSE_code_for_url <- substr(GSE_code_for_url,1,nchar(GSE_code_for_url)-3)
  GSE_code_for_url <- paste0(GSE_code_for_url, "nnn")
  complete_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", GSE_code_for_url, "/", GSE_code)
  
  checked_html_text_url <- lapply(complete_url, readUrl)
  
  
  
  #             
  if(all(checked_html_text == "EMPTY_STRING")) {
    
    cat("The web url https://ftp.ncbi.nlm.nih.gov/geo/series/ is unavailable right now. Please try again later. The function will stop here\n")
    return(NULL)
    
  } else if(all(checked_html_text_url == "EMPTY_STRING" | is.null(checked_html_text_url[[1]]) )) {
    
    cat("The web url ", complete_url," is unavailable right now (Error 404 webpage not found). The GEO code might be wrong. The function will stop here\n", sep="")
    return(NULL)        
    
  } else {
    
    gset <- retry(expr = GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE), maxErrors = 10, sleep = 2)
    
    thisGEOplatform <- toString((gset)[[1]]@annotation)
    
    if (length(gset) > 1) 
      idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    gene_expression <- as.data.frame(Biobase::exprs(gset))
    
    
    if(retrieveGeneSymbols == TRUE) {
      gene_expression$GeneSymbol <- ""
      
      
      # we retrieve the platform details
      platform_ann <- retry(expr = annotate::readGEOAnn(GEOAccNum = thisGEOplatform), maxErrors = 10, sleep = 2)
      platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
      
      if (verbose == TRUE) {
        print("sort(names(platform_ann_df))")
        print(sort(names(platform_ann_df)))
      }
      # "GENE"
      platformsWithGENEField <- c("GPL1528")
      
      # "gene_assignment
      platformsWithGene_assignmentField <- c("GPL11532", "GPL23126", "GPL6244", "GPL17586", "GPL5175")
      
      # "Gene Symbol"
      platformsWithGeneSpaceSymbolField <- c("GPL80", "GPL8300", "GPL80", "GPL96", "GPL570", "GPL571", "GPL3921")
      
      # "gene_symbol"
      platformsWithGene_SymbolField <- c("GPL20115", "GPL13667")
      
      # "HUGOname"
      platformsWithHUGOnamelField <- c("GPL5918")
      
      # "symbol"
      platformsWithSymbolField <- c("GPL1293", "GPL6102", "GPL6104", "GPL6883", "GPL6884", "GPL10558","GPL6947","GPL8177")
      
      # "GENE_SYMBOL
      platformsWith_GENE_SYMBOL_Field <- c("GPL13497", "GPL14550", "GPL17077", "GPL6480")
      
      # if symbol
      if(thisGEOplatform %in% c(platformsWithGENEField, platformsWithHUGOnamelField, 
                                platformsWithGene_assignmentField, platformsWithGeneSpaceSymbolField, 
                                platformsWithGene_SymbolField, platformsWithSymbolField, 
                                platformsWith_GENE_SYMBOL_Field)   ) {
        
        emptyGeneSymbol <- ""
        FIRST_GENE_EXPRESSION_INDEX <- 2
        
        if (verbose == TRUE)    
          cat("\n[start] loop for the association of the gene symbols to the probeset ID's\n", sep="")
        
        if(thisGEOplatform %in% platformsWithGENEField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "GENE", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithGeneSpaceSymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "Gene Symbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithGene_SymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "GeneSymbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithSymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "Symbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWith_GENE_SYMBOL_Field)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "GENE_SYMBOL", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithGene_assignmentField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "gene_assignment", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithHUGOnamelField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "HUGOname", platform_ann_df = platform_ann_df)
        gene_expression$GeneSymbol <- thisSymbol
        
        if (verbose == TRUE) 
          cat("\n [end] loop for the association of the gene symbols to the probeset ID's \n ", sep="")
        
      }  else{
        if (verbose == TRUE) 
          cat("\n\n[Impossible to perform gene symbol retrieval]\n")
        if (verbose == TRUE) 
          cat("We're sorry but the indicated platform (", thisGEOplatform, ") is not among the platforms included in this function.\nThe gene symbol retrieval cannot be performed.\nPlease visit the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", thisGEOplatform, " website for more information about this platform.\n\n", sep="")
        gene_expression$GeneSymbol <- NULL
      }
      
    }
    
    # gene expression duplicated probe remove using mdedin max
    # duplicate remove probe / gene
    gene_name <- gene_expression$GeneSymbol
    sample_name <- colnames(gene_expression)[1:(length(colnames(gene_expression)) - 1)]
    probe_MAD <- gene_expression %>% 
      select(-GeneSymbol) %>% 
      as.matrix() %>% 
      apply(.,1,mad) %>% 
      tibble(id = seq(length(.)), gene_name = gene_name, MAD = .) %>% 
      arrange(desc(MAD))
    probe_MAD_dup <- probe_MAD[which(!duplicated(probe_MAD$gene_name)),] %>% 
      filter(gene_name != "") %>% 
      arrange(id)
    gene_expression <- gene_expression[probe_MAD_dup$id, ]
    
    geneExpression_dup <- gene_expression %>% 
      select(-GeneSymbol) %>% 
      as.matrix()
    rownames(geneExpression_dup) <- gene_expression$GeneSymbol
    colnames(geneExpression_dup) <- colnames(gene_expression)[1:(length(gene_expression) - 1)]
    
    return(list(gene_expression = geneExpression_dup , 
                pheno = phenoData(gset)))
    
  }   }

# DEA & RRA function ====
run_limma <- function(ge, de){
  fit <- lmFit(ge,de)
  cont <- makeContrasts(TP-NT,levels=de) # 데이터에 맞추어 manual로 설정해야 됨
  fit.cont <- contrasts.fit(fit,cont)
  fit.cont <- eBayes(fit.cont)
  res <- topTable(fit.cont,number=Inf) 
  target <- res[res$P.Value < 0.05,] %>% 
    rownames_to_column() %>% as_tibble()
  
  return(target)
}

rra_extract <- function(ml, logfc = 0.0, fdr = 0.05){
  # combine deg
  combine_degs <- names(ml) %>% 
    lapply(X = ., FUN = function(list_name){
      tmp <- ml[[list_name]] %>% 
        dplyr::filter(adj.P.Val < fdr & (logFC >= logfc | logFC < -(logfc))) %>% 
        arrange(desc(logFC)) %>%
        dplyr::select(rowname, logFC)
      colnames(tmp) <- c("GENE", list_name)
      return(tmp)
    }) %>% 
    purrr::reduce(., full_join, by = "GENE") %>% 
    bind_cols(., apply(.[,-1], 1, mean, na.rm = TRUE) %>% 
                tibble(group = .)) 
  
  # up-regulated
  # updown_degs <- names(ml) %>% 
  #   lapply(X = ., FUN = function(list_name){
  #     ml[[list_name]] %>% 
  #       dplyr::filter(adj.P.Val < fdr & (logFC > logfc | logFC < -(logfc))) %>% 
  #       arrange(adj.P.Val) %>% 
  #       dplyr::pull(rowname) %>%
  #       return()
  #   }) 
  
  # up-regulated  
  up_degs <- names(ml) %>%
    lapply(X = ., FUN = function(list_name){
      ml[[list_name]] %>%
        dplyr::filter(adj.P.Val < fdr & (logFC >= logfc)) %>%
        arrange(adj.P.Val) %>%
        dplyr::pull(rowname) %>%
        return()
    })
  
  # down-regulated
  down_degs <- names(ml) %>%
    lapply(X = ., FUN = function(list_name){
      ml[[list_name]] %>%
        dplyr::filter(adj.P.Val < fdr & (logFC < -(logfc))) %>%
        arrange(adj.P.Val) %>%
        dplyr::pull(rowname) %>%
        return()
    })
  
  # Aggregate the inputs
  ## up-regulated RRA
  up_deg_rra <- aggregateRanks(glist = up_degs, method = "RRA") %>%
    as_tibble() %>%
    dplyr::filter(Score <= 0.01) %>% 
    arrange(Score)
  
  down_deg_rra <- aggregateRanks(glist = down_degs, method = "RRA") %>%
    as_tibble() %>%
    dplyr::filter(Score <= 0.01) %>% 
    arrange(Score)
  
  print(paste0("up-regulated : ", up_deg_rra %>% nrow()))
  print(paste0("down-regulated : ", down_deg_rra %>% nrow()))
  
  
  # run RRA
  # updown_deg_rra <- aggregateRanks(glist = updown_degs, method = "RRA") %>%
  #   as_tibble() %>%
  #   dplyr::filter(Score < 0.05) %>%
  #   arrange(Score)
  updown_deg_rra <- bind_rows(up_deg_rra, down_deg_rra) %>% 
    as_tibble() %>% 
    arrange(Score)
  
  # # 1 - combine deg, 2 - up_down-regulated RRA
  list(combine_degs = combine_degs, updown_rra = updown_deg_rra,
       up_rra = up_deg_rra$Name, down_rra = down_deg_rra$Name) %>% return()
}
rra_analysis <- function(m_list, logfc = 0, fdr = 0.05, save_path = getwd()){
  rra_result <- rra_extract(ml = m_list, logfc = logfc, fdr = fdr)
  combine_degs_rra <- rra_result$combine_degs %>% 
    dplyr::filter(GENE %in% rra_result$updown_rra$Name) %>% 
    arrange(desc(group))
  
  up_down_rra_gene <- bind_rows(head(combine_degs_rra, 20),
                                tail(combine_degs_rra, 20)) %>% 
    dplyr::select(-group)
  
  # heatmap
  combine_degs_m <- up_down_rra_gene[,-1] %>% as.matrix()
  rownames(combine_degs_m) <- up_down_rra_gene$GENE
  colnames(combine_degs_m) <- colnames(up_down_rra_gene)[2:length(colnames(up_down_rra_gene))]
  
  p <- pheatmap(combine_degs_m, 
                display_numbers = TRUE,
                number_color = "black",
                fontsize_number = 10,
                border_color = "black",
                cluster_rows = F,
                cluster_cols = F,
                cellwidth = 35,
                cellheight = 10
  )
  save_pheatmap(p, filename = paste0(save_path, "/Step2_merge_DEA.png"))
  
  # save
  save(up_down_rra_gene, file = paste0(save_path, "/Step2_merge_DEA.RData"))
  
  return(combine_degs_rra %>% 
           dplyr::pull(1))
}

save_pheatmap <- function(x, filename, width=480, height=960) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# WGCNA function ==== 
mergecutheight_test <- function(pr_name, robustdegs){
  network_list <- list()
  max_module_cnt <- 0
  max_module <- NULL
  
  for(mch_value in seq(from = 0.1, to = 0.3, by = 0.01)){
    print(paste0("mch_value : ", mch_value))
    network <- network_preprocessing(pr_name = pr_name, robustdegs = robustdegs, mch = mch_value, time_stamp = "test")
    module_cnt <- length(network[[2]]$colors %>% unique()) - 1
    print(module_cnt)
    
    if(max_module_cnt <= module_cnt){
      max_module <- mch_value
      max_module_cnt <- module_cnt
    }
    # network_list[[as.character(mch_value)]] <- list(network, module_cnt)
  }
  
  return(list(max_mch = max_module, max_module_cnt = max_module_cnt))
}
network_preprocessing <- function(pr_name, robustdegs, mch = 0.25, time_stamp){
  allowWGCNAThreads(nThreads = 20)
  
  suppressMessages({
    cancer_fpkm <- load_tcga_dataset(pkl_path = "PKL/", raw_path = "RAW_DATA/", cancer_type = pr_name) %>% 
      rownames_to_column(var = "id")
    
    gene_probe_mapping <- read_delim(file = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
                                     delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
      select(id, gene) %>% 
      inner_join(x = ., y = cancer_fpkm, by = "id") %>% 
      select(-id)
    
    dataFilt <- gene_probe_mapping %>% 
      mutate_if(is.numeric, .funs = function(value){
        2^value - 0.001 %>% return()
      }) %>% 
      mutate_if(is.numeric, .funs = function(value){
        ifelse(value < 0, 0, value) %>% return()
      }) %>% 
      mutate_if(is.numeric, .funs = function(value){
        log2(value + 1) %>% return()
      }) %>% 
      distinct(gene, .keep_all = TRUE) %>% 
      column_to_rownames(var = "gene") %>% 
      as.matrix()
    
    
    # row to col AND DESeq2 normalized log2(x+1)
    robustdeg_ge <- lapply(X = robustdegs, FUN = function(deg){
      error <- FALSE
      tryCatch(
        expr = {
          tmp <- dataFilt[deg, ]
        },
        error = function(e) {
          error <<- TRUE
        }
      )
      if(error){
        return(NULL)
      } else {
        tmp <- as.matrix(tmp) %>% t()
        rownames(tmp) <- deg
        return(tmp)
      }}) %>% do.call(rbind, .) %>% 
      t() %>% 
      as.data.frame()
    
    # split rowname
    rownames(robustdeg_ge) <- rownames(robustdeg_ge) %>% as_tibble() %>% 
      separate(col = value, into = c("A","B","C","D")) %>% 
      unite(col = id, sep = "-") %>% dplyr::pull(1)
    # gsub('.{1}$', '', .)
    
    # sample & gene filtering
    gsg <- goodSamplesGenes(robustdeg_ge, verbose = 3)
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(robustdeg_ge)[!gsg$goodGenes], collapse = ", ")));
      if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(robustdeg_ge)[!gsg$goodSamples], collapse = ", ")));
      # Remove the offending genes and samples from the data:
      robustdeg_ge <-  robustdeg_ge[gsg$goodSamples, gsg$goodGenes]
    }
    
    # power calculation
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(robustdeg_ge, powerVector = powers, verbose = 0) %>% 
      .$powerEstimate
    
    # network construct
    # net construct
    # time_stamp <- Sys.time()
    log_save <- paste(base_dir, "WGCNA_LOG", pr_name, time_stamp, sep = "/")
    dir.create(paste(base_dir, "WGCNA_LOG", pr_name, time_stamp, sep = "/"), showWarnings = FALSE, recursive = TRUE)
    net <- blockwiseModules(datExpr = robustdeg_ge, 
                            power = sft,
                            corType = "pearson",
                            TOMType = "unsigned", 
                            minModuleSize = 30,
                            # maxBlockSize = 500,
                            reassignThreshold = 0, 
                            mergeCutHeight = mch,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = paste0(log_save, "/TCGA-", pr_name),
                            verbose = 0)
    
  })
  
  
  return(list(deg = robustdeg_ge, network = net))
}
find_key_modulegene <- function(pr_name, base_dir, network, MEs, select_clinical=NULL, mm=0.85, gs=0.25, binarytocategory=FALSE){
  
  # variable
  expression_sample <- rownames(network[[1]])
  nGenes <- ncol(network[[1]])
  nSamples <- nrow(network[[1]])
  log_save <- paste(base_dir, "WGCNA_LOG/", sep = "/")
  dir.create(log_save, recursive = T, showWarnings = FALSE)
  
  # UCSCXena clinical
  clinical_trait <- read_delim(paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.",pr_name,".sampleMap%2F",pr_name,"_clinicalMatrix"),
                               delim = "\t", show_col_types = FALSE, progress = FALSE)
  survival_trait <- read_delim(paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2F",pr_name,"_survival.txt"),
                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
  
  # immune_trait <- read_delim("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Subtype_Immune_Model_Based.txt.gz",
  #                            delim = "\t", show_col_types = FALSE, progress = FALSE)
  # 
  # molecular_trait <- read_delim("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGASubtype.20170308.tsv.gz",
  #                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>%
  #   select(sampleID, Subtype_Integrative)
  
  clinical_trait <- left_join(x = clinical_trait, y = survival_trait, by = c("sampleID" = "sample"))
  # left_join(x = ., y = immune_trait, by = c("sampleID" = "sample")) %>% 
  # left_join(x = ., y = molecular_trait, by = "sampleID")
  
  default_clinical <- c('sample_type', 
                        # 'Subtype_Immune_Model_Based', 'Subtype_Integrative',
                        'pathologic_stage', 'pathologic_T', 'pathologic_N', 'pathologic_M')
  use_clinical <- c(default_clinical, select_clinical)
  
  # clinical trait preprocessing
  traitRows <- match(expression_sample, clinical_trait$sampleID)
  data_trait <- clinical_trait[traitRows, ] %>% 
    column_to_rownames(var = "sampleID") %>% 
    dplyr::select(all_of(use_clinical)) %>% 
    mutate(sample_type = ifelse(sample_type == "Primary Tumor", 1, 0),
           pathologic_stage = str_remove_all(string = pathologic_stage, pattern = "A|B|C"),
           pathologic_T = str_remove_all(string = pathologic_T, pattern = "a|b")
    ) %>%  # sample type 한정, pathogenic stage
    replace(is.na(.), 0) %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.character, as.numeric) %>% 
    mutate_if(is.character, function(value){value - 1})
  
  # remove trait element. less than 5%
  for(col_name in data_trait %>% colnames()){
    if(col_name == "age_at_initial_pathologic_diagnosis"){
      next
    }
    
    remove_trait <- data_trait[col_name] %>% 
      group_by_at(1) %>% 
      summarise(prop = n()) %>% 
      mutate(prop = prop / sum(prop)) %>% 
      filter(prop < 0.05) %>% 
      dplyr::pull(1)
    
    if(length(remove_trait) > 0){
      data_trait[[col_name]] <- ifelse(data_trait[[col_name]] %in% remove_trait, 0, data_trait[[col_name]]) 
    }
  }
  
  # binary to category
  if(binarytocategory){
    # category to binary
    dummy <- dummyVars(" ~ .", data=data_trait)
    data_trait <- data.frame(predict(dummy, newdata=data_trait))  
  }
  
  # module relation calculation 
  moduleTraitCor <-  WGCNA::cor(MEs, data_trait, use = "p")
  moduleTraitPvalue <-  corPvalueStudent(moduleTraitCor, nSamples)
  
  # module selection
  signModule <- lapply(X = 1:nrow(moduleTraitPvalue), FUN = function(row_index){
    trait_cnt <- moduleTraitPvalue[row_index, ] %>% .[. < 0.05] %>% length()
    row_name <- rownames(moduleTraitPvalue)[row_index] %>% 
      substring(., 3)
    if(str_detect(row_name, "grey"))
      return(NULL)
    else
      return(tibble(MM = row_name, signTrait = trait_cnt))
  }) %>% 
    bind_rows() %>% 
    arrange(desc(signTrait)) %>% 
    dplyr::filter(signTrait >= 2) %>% 
    dplyr::pull(1) %>% 
    # .[1:3] %>% 
    .[!is.na(.)]
  
  gene_module_key <- network[[2]]$colors %>% names() %>% 
    tibble(gene = .) %>% 
    bind_cols(., tibble(module = moduleColors)) %>% 
    filter(module %in% signModule) %>% 
    arrange(module)
  
  # Module membership
  MM <- as.data.frame(cor(network[[1]], MEs, use = "p")) 
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(MM), nSamples))
  MM <- MM %>% 
    rownames_to_column(var = "gene")
  
  # Gene significance
  # gene significance
  GS <- lapply(X =  data_trait %>% colnames(), FUN = function(s_type){
    trait <- data_trait[s_type]
    names(trait) <- s_type
    
    # gene significance
    geneTraitSignificance <- as.data.frame(cor(network[[1]], trait, use = "p")) %>%
      mutate_all(abs) %>%
      rownames_to_column("gene")
    
  }) %>% purrr::reduce(., inner_join, by = "gene") %>%
    column_to_rownames("gene") %>%
    bind_cols(., apply(., 1, median) %>% tibble(GS_median = .)) %>%
    # select(GS_median) %>%
    rownames_to_column("gene")
  
  # module 별 clinical trait 각각 calulation
  intra_module <- lapply(X = signModule, FUN = function(sig_module){
    module_gene <- gene_module_key %>% 
      filter(module == sig_module) %>% 
      dplyr::pull(gene)
    
    module_MM <- MM %>% 
      select(gene, MM = starts_with(paste0("ME", sig_module))) %>% 
      filter(gene %in% module_gene, abs(MM) > mm)
    module_MM_gene <- module_MM %>% dplyr::pull(gene)
    
    module_MM_GS_filtered <- lapply(X = data_trait %>% colnames(), FUN = function(t){
      GS %>% 
        select(gene, GS = t) %>% 
        filter(gene %in% module_MM_gene, abs(GS) > gs) %>% 
        arrange(desc(abs(GS)))
      
    })
    names(module_MM_GS_filtered) <- data_trait %>% colnames()
    
    return(module_MM_GS_filtered)
  })
  names(intra_module) <- signModule
  
  intra_module_size <- lapply(X = intra_module %>% names(), function(im_name){
    intra_module[[im_name]] %>% bind_rows() %>% 
      distinct(gene, .keep_all = TRUE) %>% 
      nrow() %>% 
      tibble(module = im_name, module_size = .) %>% 
      return()
  }) %>% bind_rows() %>% 
    arrange(desc(module_size))
  
  intra_module_bind <- lapply(X = intra_module %>% names(), function(im_name){
    intra_module[[im_name]] %>% bind_rows() %>% 
      distinct(gene, .keep_all = TRUE) %>%
      arrange(gene) %>% 
      return()
  })
  names(intra_module_bind) <- signModule
  
  # # module eigen gene intersection
  # total_keyhub <- list()
  # for(index in data_trait %>% colnames()){
  #   tmp <- lapply(X = intra_module, FUN = function(trait){
  #     trait[[index]]
  #   }) %>% bind_rows() %>% dplyr::pull(gene)
  #   
  #   if(length(tmp) > 1)
  #     total_keyhub[[index]] <- tmp
  # }
  # 
  # total_keyhub_merge <- names(total_keyhub) %>% lapply(X = ., FUN = function(lname){
  #   col_name <- lname %>% str_split(pattern = "\\.") %>% unlist %>% .[1]
  #   
  #   tibble(gene = total_keyhub[[lname]], col_name) %>% return()
  # }) %>% bind_rows() %>% 
  #   split(x = ., f = .$col_name) %>% 
  #   lapply(X = ., FUN = function(df){df %>% dplyr::pull(gene)})
  
  # total_keyhub <- total_keyhub_merge %>% unlist() %>% unname() %>% unique()
  
  # plot save
  sample_cluster_plot(network = network[[1]], clinical_trait = data_trait, save_path = log_save)
  module_cluster_plot(network = network, save_path = log_save)
  module_trait_plot(moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue,
                    data_trait = data_trait, MEs = MEs, save_path = log_save)
  gene_module_size_plot(gene_module_key_groph = intra_module_size, 
                        title = "After Intra Module analysis(MM-GS cut-off), Module Size",
                        prefix = "Step3_module_size",
                        save_path = log_save)
  # key_hub_intersection_plot(total_keyhub = total_keyhub_merge, save_path = log_save)
  
  
  
  # return(list(total_keyhub = total_keyhub, clinical_trait = data_trait, total_keyhub_merge = total_keyhub_merge, network = network))
  return(list(intra_module = intra_module, intra_module_gene = intra_module_bind, clinical_trait = data_trait, network = network))
  
}

sample_cluster_plot <- function(network, clinical_trait, save_path){
  sampleTree <-  hclust(dist(network[[1]]), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  # plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
  #      cex.axis = 1.5, cex.main = 2)
  
  traitColors <- numbers2colors(clinical_trait, signed = FALSE);
  
  png(paste0(save_path, "/sample_cluster.png"), width = 1000, height = 1000)
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(clinical_trait),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  
  
}
module_cluster_plot <- function(network, save_path){
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors <- labels2colors(network[[2]]$colors)
  # Plot the dendrogram and the module colors underneath
  
  png(paste0(save_path, "/Step3_module_cluster.png"), width = 1000, height = 1000)
  plotDendroAndColors(network[[2]]$dendrograms[[1]], mergedColors[network[[2]]$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)  
  dev.off()
}
module_trait_plot <- function(moduleTraitCor, moduleTraitPvalue, data_trait, MEs, save_path){
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  png(paste0(save_path, "/Step3_module_trait.png"), width = 1000, height = 1000)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(data_trait),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
}
gene_module_size_plot <- function(gene_module_key_groph, title, prefix, save_path){
  p <- ggpubr::ggbarplot(gene_module_key_groph, x = "module", y = "module_size",
                         color = "module",
                         fill = "module",
                         sort.val = "asc",
                         sort.by.groups = FALSE, 
                         x.text.angle = 90,          # Rotate vertically x axis texts
                         title = title,
                         ylab = "Module Size",
                         legend.title = "Module",
                         rotate = TRUE,
                         label = TRUE, lab.pos = "in", lab.col = "black",
                         ggtheme = theme_minimal()) + theme(legend.position="none")
  ggsave(p, filename = paste0(save_path, "/", prefix, ".png"), 
         width = 8, height = 8)
  return(p)
}
key_hub_intersection_plot <- function(total_keyhub, save_path){
  
  p <- ggVennDiagram::ggVennDiagram(x = total_keyhub, label_size = 7) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  
  ggsave(p, filename = paste0(save_path, "/Step3_keyhubgene_intersection.png"), dpi = 200, width = 30, height = 10)
  return(p)
}

# # plot - GS & MM correlation 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "GS for sample type",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Machine learning function ====
gene_selection <- function(base_dir, total_keyhub_list, over_sampling){
  
  log_save <- paste(base_dir, "ML_LOG/", sep = "/")
  dir.create(log_save, recursive = T, showWarnings = FALSE)
  
  total_keyhub_list_filter <- total_keyhub_list %>% lapply(X = ., FUN = function(df){
    if(nrow(df) <= 0)
      return(NULL)
    else{
      return(df)
    }
  }) %>% compact()
  
  trait_names <- total_keyhub_list_filter %>% names()
  gene_selection_list <- list()
  
  for(trait_name in trait_names){
    Y_col_name <- trait_name
    DF <- clinical_trait %>% 
      dplyr::select(all_of(Y_col_name)) %>%
      rownames_to_column(var = "sample") %>% 
      inner_join(x = ., y = geneExpression %>% 
                   rownames_to_column(var = "sample") %>% 
                   dplyr::select(sample, all_of(total_keyhub_list[[Y_col_name]]$gene)),
                 by = "sample") %>% 
      column_to_rownames(var = "sample")
    y_df <- DF %>% select(-all_of(total_keyhub_list[[Y_col_name]]$gene))
    x_df <- DF %>% select_at(all_of(total_keyhub_list[[Y_col_name]]$gene))
    if(ncol(x_df) <= 1){
      next
    }
    
    # gene selection
    write_csv(x_df, file = paste0(log_save, "x_df_temp.csv"))
    write_csv(y_df, file = paste0(log_save, "y_df_temp.csv"))
    
    system(glue("python3 src/py-lasso.py -b {base_dir} -x {x_path} -y {y_path} -o True", 
                base_dir = log_save,
                x_path = paste0(log_save, "x_df_temp.csv"), 
                y_path = paste0(log_save, "y_df_temp.csv"))
           )
    lasso_coef <- read_csv(file = paste0(log_save, "/lasso_result.csv"), show_col_types = FALSE) %>% pull(coef)
    # lasso_coef <- feature_selection_LASSO(x_df, y_df, over_sampling) # Invoke python ########## 실행시 crash 발생함
    lasso_selection_gene <- x_df %>% select(which(abs(lasso_coef) > 0)) %>% colnames()
    
    gene_selection_list[[trait_name]] <- lasso_selection_gene
  }
  
  # venn diagram
  # seleted_gene_intersection_plot(gene_selection_list, log_save)
  
  
  return(gene_selection_list)
  
}

ml_validation <- function(base_dir, selected_gene, over_sampling, cv, module_name){
  
  log_save <- paste(base_dir, "ML_LOG/", sep = "/")
  dir.create(log_save, recursive = T, showWarnings = FALSE)
  trait_names <- selected_gene %>% names()
  
  roc_auc_list <- list()
  for(trait_name in trait_names){
    Y_col_name <- trait_name
    DF <- clinical_trait %>% 
      select(Y_col_name) %>%
      rownames_to_column(var = "sample") %>% 
      inner_join(x = ., y = geneExpression %>% 
                   rownames_to_column(var = "sample") %>% 
                   select(sample, all_of(selected_gene[[Y_col_name]])),
                 by = "sample") %>% 
      column_to_rownames(var = "sample")
    
    write_csv(DF, file = paste0(log_save, "df_temp.csv"))
    
    
    # gene selection
    if(!cv){
      # roc_auc_list[[trait_name]] <- roc_acu_calculator(DF, Y_col_name, log_save, over_sampling = over_sampling, module_name = module_name)
      system(glue("python3 src/py-roc.py -b {base_dir} -d {df_path} -l {log_save} -f {f_name} -m {m_name} -o True", 
                  base_dir = log_save,
                  log_save = log_save,
                  df_path = paste0(log_save, "df_temp.csv"), 
                  f_name = Y_col_name,
                  m_name = module_name)
      )
      roc_auc_list[[trait_name]] <- jsonlite::read_json(paste0(log_save, 'ml_validation_result.json'))
      
    } else {
      # roc_auc_list[[trait_name]] <- roc_acu_calculator_cv(DF, Y_col_name, log_save, over_sampling = over_sampling, module_name = module_name)
      system(glue("python3 src/py-roc.py -b {base_dir} -d {df_path} -l {log_save} -f {f_name} -m {m_name} -o True", 
                  base_dir = log_save,
                  log_save = log_save,
                  df_path = paste0(log_save, "df_temp.csv"), 
                  f_name = Y_col_name,
                  m_name = module_name)
      )
      roc_auc_list[[trait_name]] <- jsonlite::read_json(paste0(log_save, 'ml_validation_result.json'))
      
    }
  }
  
  return(roc_aut_result = roc_auc_list)
  
}
auc_cutoff <- function(sg_list, selected_gene, auc_cutoff = 0.75){
  clinical_name <- sg_list %>% names()
  auc_cutoff_gene <- lapply(X = clinical_name, FUN = function(list_name){
    value <- sg_list[[list_name]]
    
    if(length(value) >= 2){
      if(value$micro >= auc_cutoff){
        return(list_name)
      }
    } else {
      if(value >= auc_cutoff){
        return(list_name)
      }
    }
  }) %>% unlist() %>% 
    selected_gene[.] %>% 
    unlist() %>% 
    unique()
  
  return(list(auc_cutoff = auc_cutoff_gene, trait = clinical_name))
}

seleted_gene_intersection_plot <- function(gene_selection_list, save_path){
  
  p <- ggVennDiagram::ggVennDiagram(x = gene_selection_list, label_size = 7) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  
  ggsave(p, filename = paste0(save_path, "/Step4_seleted_gene_intersection.png"), dpi = 200, width = 30, height = 10)
  return(p)
}


# STRING function ====
retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=5, sleep=0) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
      flog.fatal(msg)
      stop(msg)
    } else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors, 
                    capture.output(str(retval)))
      flog.error(msg)
      warning(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(retval)
}
#' Function that returns STRING network
#' @param hub_gene input character vector
#' @return network dataframe
#' @examples
#' aaa <- string_network(hub_gene = my_gene)
string_network <- function(hub_gene){
  # URLs
  string_api_url <- "https://version-11-5.string-db.org/api"
  output_format <- "tsv-no-header"
  method <- "network"
  request_url <- paste(string_api_url, output_format, method, sep = "/")
  
  # post payload
  params <- list(
    identifiers = paste0(hub_gene, collapse = "%0d"),
    species = "9606",
    caller_identity = "www.hallym.ac.kr"
  )
  
  # output
  network_colname <- c('stringId_A','stringId_B','preferredName_A','preferredName_B','ncbiTaxonId',
                       'score','nscore','fscore','pscore','ascore','escore','dscore','tscore')
  
  ppi_network <- httr::POST(request_url, body = params, encode = "form") %>%  
    httr::content(encoding = "UTF-8") 
  colnames(ppi_network) <- network_colname
  ppi_network <- ppi_network %>% arrange(preferredName_A) %>% distinct_all()
  
  return(ppi_network)
}

# Filtering ----
# DGIdb
dgidb_interaction_parallel <- function(gene_name){
  base_url <- "https://dgidb.org"
  request_url <- paste0(base_url, "/api/v2/interactions.json?")
  result_list <- list()
  
  # chunk id
  id_chunk <- split(gene_name, ceiling(seq_along(gene_name)/200))
  
  mclapply(X = 1:length(id_chunk), FUN = function(index){
    # print(index)
    payload <-  list(genes = paste0(id_chunk[[index]], collapse = ","),
                     fda_approved_drug="true")
    
    # output
    dgidb_result <- POST(request_url, body = payload, encode = "form", config = httr::config(connecttimeout = 60)) %>%  
      httr::content(encoding = "UTF-8") 
    
    lapply(X = dgidb_result$matchedTerms, FUN = function(dgidb_element){
      gene_category <- dgidb_element$geneCategories %>% 
        sapply(X = ., FUN = function(value) {value$name}) %>% 
        paste0(collapse = ",")
      
      interaction <- dgidb_element$interactions %>% 
        sapply(X = ., FUN = function(value){
          drug_name <- value$drugName
          score <- value$score
          types <- value$interactionTypes %>% unlist() %>% paste0(collapse = "&")
          
          paste0(c(drug_name, score, types), collapse = ";") %>% 
            as_tibble() %>% 
            return()
          
          # return(drug_name)  
        }) %>% unlist() %>% 
        paste0(., collapse = "&")
      
      tibble(
        gene = dgidb_element$geneName,
        DGI_GENE_CATEGORY = gene_category, 
        `DGI(DRUG_NAME;SCORE;TYPE)` = interaction,
        DGI_COUNT = length(dgidb_element$interactions)
      )  %>% return()
      
    }) %>% return()
  }, mc.cores = 3) %>% bind_rows() %>% return()
}

# Textmining
textmining_extract <- function(pr_name){
  con <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", 
                        port = 3306, user = "root", password = "sempre813!", dbname = "Textmining")
  
  tm <- tbl(con, pr_name) %>% collect()
  
  return(tm)  
}

# NT-TP DEA
run_deseq_normal <- function(pr_name, base_dir, batch_removal){
  register(MulticoreParam(20))
  suppressMessages({
    if((!file.exists(paste0(base_dir, "/GDCdata/", pr_name, "_normal.RData"))) | 
       (!file.exists(paste0(base_dir, "/GDCdata/", pr_name, "_RnaseqSE_normal.RData")))){
      query <- GDCquery(project = paste0("TCGA-", pr_name), 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        experimental.strategy = "RNA-Seq",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
                        legacy = TRUE)
      
      GDCdownload(query, directory = paste0(base_dir, "/GDCdata"))
      RnaseqSE <- GDCprepare(query)
      
      save(RnaseqSE, file = paste0(base_dir, "/GDCdata/", pr_name, "_RnaseqSE_normal.RData"))
      
      Rnaseq_CorOutliers <- assay(RnaseqSE, "raw_count") # to matrix
      
      # normalization of genes, # quantile filter of genes
      dataNorm <- TCGAanalyze_Normalization(tabDF = Rnaseq_CorOutliers, geneInfo =  geneInfo)
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)
      
      save(dataFilt, file = paste0(base_dir, "/GDCdata/", pr_name, "_normal.RData"))
    } else {
      load(paste0(base_dir, "/GDCdata/", pr_name, "_RnaseqSE_normal.RData"))
      load(paste0(base_dir, "/GDCdata/", pr_name, "_normal.RData"))
    }
    
    
    # selection of normal samples "NT"
    samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("NT")) %>% 
      as_tibble() %>% 
      mutate(group = 0) %>% 
      dplyr::rename(sample = value)
    
    # selection of tumor samples "TP"
    samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("TP")) %>% 
      as_tibble() %>% 
      mutate(group = 1) %>% 
      dplyr::rename(sample = value)
    
    metadata <- bind_rows(samplesNT, samplesTP) %>% 
      mutate(group = ifelse(group == 0, "NT", "TP"))
    metadata$group <- factor(metadata$group, levels = c("NT", "TP"))
    
    tcga_se <- DESeqDataSetFromMatrix(countData = dataFilt, colData = metadata, design = ~ group)
    tcga_se$group <- relevel(tcga_se$group, ref = "NT")
    tcga_deseq <- DESeq(tcga_se, parallel = TRUE)
    
    # Hidden Batch effect remove
    if(batch_removal){
      # hidden batch removal
      dat <- counts(tcga_deseq, normalized=TRUE)
      mod <- model.matrix(~ group, colData(tcga_deseq))
      mod0 <- model.matrix(~ 1, colData(tcga_deseq))
      
      # run sva
      svseq <- svaseq(dat, mod, mod0, n.sv=2)
      tcga_se_batch <-tcga_se
      
      tcga_se_batch$SV1 <- svseq$sv[,1]
      tcga_se_batch$SV2 <- svseq$sv[,2]
      design(tcga_se_batch) <- ~ SV1 + SV2 + group
      tcga_deseq <- DESeq(tcga_se_batch, parallel = TRUE)
    }
    
    tcga_deseq_result <- results(tcga_deseq, 
                                 alpha = 0.9999)
    tcga_deseq_result_tidy <- results(tcga_deseq, 
                                      alpha = 0.9999,
                                      tidy = TRUE)
    
    # volcano plot
    p <- EnhancedVolcano(tcga_deseq_result,
                         lab = rownames(tcga_deseq_result),
                         x = 'log2FoldChange',
                         y = 'padj',
                         title = 'Primary Tumor versus Solid Tissue Normal',
                         pCutoff = 0.05,
                         FCcutoff = 1.5,
                         pointSize = 3.0,
                         labSize = 6.0)
    
    ggsave(plot = p, filename = paste0(base_dir, "/GDCdata/", pr_name, "_volcano/", pr_name, "_DESEQ2_normal_volcano.png"), height = 8, width = 12, dpi = 70)    
  })  
  
  return(tcga_deseq_result_tidy)
}

# protein atlas
protein_atlas <- function(gene_name){
  suppressMessages({
    protein_atlas_url <- "https://www.proteinatlas.org/"
    gene_list_mapping <- mapIds(org.Hs.eg.db,
                                keys=gene_name, 
                                column="ENSEMBL",
                                keytype="SYMBOL",
                                multiVals="first") %>% 
      tibble(gene = names(.), ENSEMBL = .) %>% 
      mutate(PROTEIN_ATLAS = ifelse(!is.na(ENSEMBL), paste0(protein_atlas_url, ENSEMBL, "-", gene), "")) %>% 
      dplyr::select(gene, PROTEIN_ATLAS)
  })
  return(gene_list_mapping)
}

further_analysis <- function(mc, ge, key_hub_gene, base_dir){
  # STRING analysis
  string_filtering <- string_analysis(mc, base_dir) %>% 
    compact() %>% 
    lapply(X = ., FUN = function(df){
      df %>% mutate(STRING = TRUE)
    }) %>% bind_rows()
  
  # Survival analysis
  survival_filtering <- survival_analysis(base_dir = base_dir, 
                                          geneExpression = ge, 
                                          mc = mc) %>% 
    bind_rows()
  
  # DNA metylation analysis
  methylation_filtering <- methylation_analysis(pr_name = pr_name, base_dir = base_dir)
  methylation_filtering_ <- lapply(X = names(methylation_filtering), FUN = function(method){
    GENE_NAME <- methylation_filtering[[method]]$Symbol
    P.adjust <- methylation_filtering[[method]]$Pe
    
    tmp_df <- tibble(GENE_NAME = GENE_NAME, Methylated = method, P.adjust = P.adjust) %>% 
      distinct(GENE_NAME, .keep_all = TRUE)
    colnames(tmp_df) <- c("GENE_NAME", paste0('Methylated_', method),
                          paste0('P.adjust_', method))
    return(tmp_df)
  })
  names(methylation_filtering_) <- names(methylation_filtering)
  
  # analysis merge
  module_candidate_analysis <- mc %>% bind_rows() %>% 
    left_join(x = ., y = key_hub_gene$intra_module_gene %>% bind_rows(), by = c("GENE_NAME" = "gene")) %>% 
    left_join(x = ., y = string_filtering, by = c("GENE_NAME", "TRAIT", "MODULE")) %>% 
    left_join(x = ., y = survival_filtering, by = c("GENE_NAME", "TRAIT", "MODULE")) %>% 
    left_join(x = ., y = methylation_filtering_$hypo, by = "GENE_NAME") %>% 
    left_join(x = ., y = methylation_filtering_$hyper, by = "GENE_NAME") %>% 
    unite(Methylated, Methylated_hypo, Methylated_hyper, sep = "/", na.rm = TRUE) %>% 
    unite(Methylated.Padj, P.adjust_hypo, P.adjust_hyper, sep = "/", na.rm = TRUE) %>% 
    rename(GENE_SIGNIFICANCE = GS)
  
  # Enrichment analysis
  lapply(X = names(mc), FUN = function(module_name){
    ora_go_kegg(geneName = mc[[module_name]]$GENE_NAME,
                module_name = module_name,
                base_dir = base_dir)
  })
  
  return(module_candidate_analysis)
}

filtering_combine <- function(pr_name, mc){
  # filtering merge
  # nt_tp <- run_deseq_normal(pr_name = pr_name, base_dir = base_dir, batch_removal = FALSE) %>% 
  #   select('row', 'log2FoldChange', 'pvalue')
  # colnames(nt_tp) <- c('gene', 'NT-TP_Log2FC', 'NT-TP.Padj')
  tm <- textmining_extract(pr_name = pr_name)
  colnames(tm) <- c('gene', 'TM.type', 'TM.Support', 'TM.Confidence', 'TM.Lift', 'TM.Count')
  
  dgi <- dgidb_interaction_parallel(mc$GENE_NAME)
  colnames(dgi) <- c('gene', 'DGI.Gene Category', 'DGI.DrugName;Score;Type', 'DGI.Count')
  
  pdb <- symbol2pdb(mc$GENE_NAME)
  colnames(pdb) <- c('gene', 'PDB.Id', 'PDB.Count')
  
  oncokb <- oncokb_allcuratedGenes()
  colnames(oncokb) <- c('gene', 'OncoKB.Is Oncogene', 'OncoKB.Is TSG', 
                        'OncoKB.Highest Level of Evidence(sensitivity)',
                        'OncoKB.Highest Level of Evidence(resistance)',
                        'OncoKB.Background')
  
  protein <- protein_atlas(gene_name = mc$GENE_NAME)
  
  mc %>% 
    # left_join(x = ., y = nt_tp, by = c('GENE_NAME' = 'gene')) %>% 
    left_join(x = ., y = tm %>% select(-`TM.type`), by = c('GENE_NAME' = 'gene')) %>% 
    left_join(x = ., y = dgi, by = c('GENE_NAME' = 'gene')) %>% 
    left_join(x = ., y = pdb, by = c('GENE_NAME' = 'gene')) %>% 
    left_join(x = ., y = oncokb, by = c('GENE_NAME' = 'gene')) %>% 
    left_join(x = ., y = protein, by = c('GENE_NAME' = 'gene')) %>% 
    return()
}
