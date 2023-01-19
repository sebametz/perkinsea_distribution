##############################################################################################################
# Script: get_tables.R
# Created: 21/10/2022
# Author: Metz Sebastian
# Objective: get the Perkinsea tables, abundance, metadata and taxonomy
# Modified: 21/10/2022
# NOTE: for data processioning please refer to the manuscript
##############################################################################################################
# Libreries ----
library(tidyverse) # for data table manipulation
library(vroom) # for acelerate the table upload

# 1) get abundance table = ASVs as row and samples as Columns
# Extract the ASVs from two different tables from EukBank, one with Perkinsea ASVs and the other off 
# unknown Alveolata therefore we need to check both tables for the ASVs that we have and create only one table

# uploading ASVs from the Phylogenetic classification
taxonomy <- vroom("data/CLASSIFICATION/taxonomy_per_query.tsv", delim = "\t")

# uploading Perkinsea table
perkinsea <- vroom("data/EUKBANK_TABLES/EukBank_2020_10_27_Perkinsea/eukbank_18SV4_asv.table", num_threads = 8)

# uploading unknown Alveolata table
alveolata <- vroom("data/EUKBANK_TABLES/EukBank_2020_10_27_unknown_alveolata/eukbank_18SV4_asv.table")

# creating a new table
# amplicon format of the tables don't have the ;size=XX as we have in the taxonomy file in name columns
perkinsea_new_table <- tibble(amplicon = str_extract(taxonomy$name, "\\w+")) %>%
  left_join(add_case(perkinsea, alveolata), by = "amplicon") %>% # merge Perkinsea and alveolata and filter ASVs
  select(amplicon, where( ~ is.numeric(.x) && sum(.x) != 0)) # select the samples with Perkinsea reads

# 2) get the metadata file 
eukbank_meta <- read_tsv("data/EUKBANK_TABLES/eukbank_18SV4_asv.metadata_28062021")
perkinsea_meta <- filter(eukbank_meta, sample %in% colnames(perkinsea_new_table))

# there are 758 samples that don't have metadata, most of them are repeated samples, 
# therefore we will delete them

perkinsea_new_table_filtered <- perkinsea_new_table[,colnames(perkinsea_new_table) %in% c("amplicon", perkinsea_meta$sample)]

# write tables
write_tsv(perkinsea_new_table_filtered, "data/PERKINSEA_TABLES/Perkinsea_abundance.table")
write_tsv(perkinsea_meta, "data/PERKINSEA_TABLES/Perkinsea_metadata.table")

# 3) prepare taxonomy table
ncol <- max(str_count(taxonomy$taxopath, ";")) + 1
colnm <- str_c("TAX", 1:ncol)
  
taxonomyF <- taxonomy %>%
  separate(col = taxopath, sep = ";", into = colnm, remove = F)

colnames(taxonomyF)[1] <- "amplicon"
write_tsv(taxonomyF, "data/PERKINSEA_TABLES/Perkinsea_taxonomy.table")

# The table of samples were then manually revised and the samples were 
# classified according to their environmental type

# 4) get rarefied tables <- from the eukbank database
table <- read_tsv("data/EUKBANK_TABLES/Perkinsea_18SV4_rareffied_10K.table")
