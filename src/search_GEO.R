setwd("~/RRA-WGCNA/")

library(GEOmetadb)
library(tidyverse)
library(reshape2)

con <- dbConnect(SQLite(), "RAW_DATA/GEOmetadb.sqlite") # 0607 version
geo_tables <- dbListTables(con)

# sample count
sample_count <- tbl(con, "gse_gsm") %>% 
  collect() %>% 
  group_by(gse) %>% 
  summarise(sample_count = n())

#gpl_count
gpl_count <- tbl(con, "gse_gpl") %>% 
  collect() %>% 
  group_by(gse) %>% 
  summarise(gpl_count = n())

gpl <- tbl(con, "gse_gpl") %>% 
  collect() 

gpl_id <- lapply(X = gpl %>%  pull(gse) %>% unique(),
       FUN = function(value){
         tmp <- gpl %>% filter(gse == value)
         return(tmp %>% pull(gpl) %>% paste0(collapse = ";"))
       })

# query
query <- paste("SELECT gse, title, summary, type, pubmed_id FROM gse", 
               "WHERE (title LIKE '%hepatocellular%' OR summary LIKE '%hepatocellular%'",
               "OR title LIKE '%liver%' OR summary LIKE '%liver%')",
               "AND (title LIKE '%cancer%' OR summary LIKE '%cancer%' OR title LIKE '%carcinoma%' OR summary LIKE '%carcinoma%' )", 
               "AND type = 'Expression profiling by array'",
               "AND gse IN (SELECT DISTINCT series_id FROM gsm WHERE organism_ch1 == 'Homo sapiens')",
               sep = " ")

rs <- dbGetQuery(conn = con, statement = query) %>% 
  as_tibble()

# 1st filter -> sample_count
rs_sample <- inner_join(rs, sample_count, by = "gse") %>% 
  filter(sample_count >= 30) %>% 
  arrange(desc(sample_count))

# 2nd filter -> gpl count
rs_gpl <- inner_join(rs_sample, gpl_count, by = "gse") %>% 
  filter(gpl_count <= 4) %>% 
  # filter(!is.na(pubmed_id)) %>% 
  arrange(gpl_count)

# Write with add GEO LINKS
rs_gpl %>% 
  mutate(GEO_LINK = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse)) %>% 
  write_delim("GEO_2nd_filtered_HCC.csv", delim = ",")
