# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
source("revigo.R")


# Set paths --------------------------------------------------------------------
data_dir <- file.path("results", "scRNAseq", "01_DGEA")
res_dir <- file.path("results", "scRNAseq", "02_ORA")


# Input data -------------------------------------------------------------------
## DGEA results 
degs_files <- list.files(data_dir, pattern = "*.degs.txt", full.names = TRUE)

degs <- map(degs_files, read_tsv, id = "dataset") %>% 
  list_rbind() %>% 
  mutate(
    dataset = str_remove(basename(dataset), ".degs.txt"),
    status = case_when(
      (p_val_adj <= 0.05) & (avg_log2FC >= 1) ~ "Up",
      (p_val_adj <= 0.05) & (avg_log2FC <= (-1)) ~ "Down",
      .default = "Stable")) 


# Filter DGEA results ----------------------------------------------------------
degs_flt <- degs %>% 
  filter(status != "Stable") %>% 
  dplyr::select(gene, status, celltype, dataset)

## Save combined table with DEGs from all datasets
write_tsv(degs_flt, file.path(res_dir, "degs_combined_table.txt"))


# Group DEGs by cell type (referred to as subtype in SCovid DB) and dataset ----
## Make nested column for dataset:tissue:UpDown groupping 
degs_by_celltype <- degs_flt %>% 
  # Rename mitochondrial genes to be found in org.Hs.eg.db
  mutate(gene = str_remove(gene, "^MT-")) %>% 
  group_by(celltype, dataset, status) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(celltype_dataset_status = str_c(celltype, dataset, status, sep = "//"))


# Functional enrichment --------------------------------------------------------
## Function for GO BP Over-Representation Analysis (ORA)
run_enrichGO <- function(data){
  data$go_bp = vector(mode = "list", length = nrow(data))
  for (i in 1:nrow(data)){
    data$go_bp[[i]] = as.data.frame(enrichGO(gene = data$data[[i]]$gene, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP"))
  }
  return(data)
}

## Run ORA
go_by_celltype <- run_enrichGO(degs_by_celltype)


## Save RDS object with ORA results 
saveRDS(go_by_celltype, file.path(res_dir, "go_by_celltype.rds"))


# Filter GO BP terms -----------------------------------------------------------
## Parse ORA results - add dataset, cell type, and Up/Down status info
go_by_celltype_df <- go_by_celltype$go_bp %>% 
  map2(go_by_celltype$celltype_dataset_status, ., ~ add_column(.y, celltype_dataset_status = .x)) %>% 
  list_rbind() %>% 
  separate(celltype_dataset_status, into = c("celltype", "dataset", "status"), sep = "//")


## Lipid-related words
lipid_terms <- c("lipid", "fat", "triglyceride", "cholesterol")
nonlipid_terms <- c("fate", "sulfat")

go_lipid <- go_by_celltype_df %>%
  # Keep only significant GO BP terms
  dplyr::filter(qvalue <= 0.05) %>% 
  # Keep only GO BP terms related to lipids metabolism
  dplyr::filter(str_detect(Description, paste(lipid_terms, collapse = "|"))) %>%
  # Remove some GO BP terms
  dplyr::filter(!str_detect(Description, paste(nonlipid_terms, collapse = "|")))


# Run Revigo - discard redundant GO terms --------------------------------------
go_lipid_revigo <- go_lipid %>% 
  dplyr::select(ID, qvalue) %>% 
  arrange(qvalue) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  revigo(cutoff = "small") 


# Discard redundant GO terms from all ORA results ------------------------------
revigo_flt <- go_lipid_revigo %>% 
  dplyr::filter(Representative == "null")

go_lipid_revigo_flt <- go_lipid %>% 
  dplyr::filter(ID %in% revigo_flt$TermID)

## Save Revigo-filtered ORA results
write_tsv(go_lipid_revigo_flt, file.path(res_dir, "go_lipid_revigo_flt.txt"))


# Substitute discarded GO terms with representative ones -----------------------
go_lipid_map <- go_lipid %>% 
  dplyr::select(ID, Description_revigo = Description) %>% 
  distinct()

go_lipid_revigo_imp <- go_lipid_revigo %>% 
  mutate(TermID_revigo = if_else(Representative == "null", 
                                 TermID,
                                 paste0("GO:", str_pad(Representative, width = 7, side = "left", pad = "0")))) %>% 
  left_join(go_lipid_map, by = c("TermID_revigo" = "ID"))

go_lipid_revigo_imp_res <- go_lipid %>% 
  left_join(dplyr::select(go_lipid_revigo_imp, TermID, TermID_revigo, Description_revigo), by = c("ID" = "TermID"))

## Save filtered ORA results 
write_tsv(go_lipid_revigo_imp_res, file.path(res_dir, "go_lipid_revigo_imp_res.txt"))

