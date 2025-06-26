# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
# install.packages("./KEGG.db_1.0.tar.gz", repos = NULL, type = "source")
library(KEGG.db)
source("revigo.R")


# Set paths --------------------------------------------------------------------
data_dir <- file.path("results", "scRNAseq", "01_DGEA")
res_dir <- file.path("results", "scRNAseq", "02_ORA")


# Input data -------------------------------------------------------------------
## Datasets meta
ds_meta <- read_tsv("sc_dataset_meta.txt")

### Exclude datasets of organoids/cells infected with SARS-CoV-2 in vitro
datasets_2exclude <- c(
  "GSE159556", "GSE182298", "GSE167747", "GSE178404", "GSE208034", "GSE166766",
  "GSE156760", "GSE151878")


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
  dplyr::select(gene, status, celltype, dataset) %>% 
  # Rename mitochondrial genes to be found in org.Hs.eg.db
  mutate(gene = str_remove(gene, "^MT-")) 
  

## Convert SYMBOL to ENTREZ IDs ------------------------------------------------
degs_entrez <- bitr(degs_flt$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(SYMBOL)) 

degs_flt <- degs_flt %>% 
  left_join(degs_entrez, by = c("gene" ="SYMBOL"))


## Save combined table with DEGs from all datasets
write_tsv(degs_flt, file.path(res_dir, "degs_combined_table.txt"))


# Group DEGs by cell type (referred to as subtype in SCovid DB) and dataset ----
## Make nested column for dataset:tissue:UpDown groupping 
degs_by_celltype <- degs_flt %>% 
  group_by(celltype, dataset, status) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(celltype_dataset_status = str_c(celltype, dataset, status, sep = "//"))


# Functional enrichment --------------------------------------------------------
## Function for GO BP Over-Representation Analysis (ORA)
run_enrichGO <- function(data){
  go_bp = as.data.frame(enrichGO(gene = data$gene, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP"))
  return(go_bp)
}


## Function for KEGG Over-Representation Analysis (ORA)
run_enrichKEGG <- function(data){
  gene_list = data %>% drop_na(ENTREZID) %>% pull(ENTREZID)
  kegg = enrichKEGG(gene = gene_list, organism = "hsa", use_internal_data = TRUE)
  # Convert EntrezIDs to gene symbols
  if (!is.null(kegg)) {
    kegg = setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  }
  kegg = as.data.frame(kegg)
  return(kegg)
}


## Run ORA ---------------------------------------------------------------------
ora_by_celltype$go_bp <- map(degs_by_celltype$data, run_enrichGO)
ora_by_celltype$kegg <- map(ora_by_celltype$data, run_enrichKEGG)


## Save RDS object with ORA results --------------------------------------------
saveRDS(ora_by_celltype, file.path(res_dir, "ora_by_celltype.rds"))


# Filter GO BP terms -----------------------------------------------------------
## Parse ORA results - add dataset, cell type, and Up/Down status info
go_by_celltype_df <- ora_by_celltype$go_bp %>% 
  map2(ora_by_celltype$celltype_dataset_status, ., ~ add_column(.y, celltype_dataset_status = .x)) %>% 
  list_rbind() %>% 
  separate(celltype_dataset_status, into = c("celltype", "dataset", "status"), sep = "//")


## Lipid-related words
lipid_terms <- c("lipid", "fat", "triglyceride", "cholesterol")
nonlipid_terms <- c("fate", "sulfat")

go_lipid <- go_by_celltype_df %>%
  # Remove 'in vitro' datasets 
  dplyr::filter(str_detect(dataset, paste(datasets_2exclude, collapse = "|"), negate = TRUE)) %>% 
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


# Filter KEGG pathways ---------------------------------------------------------
## Parse ORA results - add dataset, cell type, and Up/Down status info
kegg_by_celltype_df <- ora_by_celltype$kegg %>% 
  map2(ora_by_celltype$celltype_dataset_status, ., ~ add_column(.y, celltype_dataset_status = .x)) %>% 
  list_rbind() %>% 
  separate(celltype_dataset_status, into = c("celltype", "dataset", "status"), sep = "//")


kegg_lipid <- kegg_by_celltype_df %>%
  # Remove 'in vitro' datasets 
  dplyr::filter(str_detect(dataset, paste(datasets_2exclude, collapse = "|"), negate = TRUE)) %>% 
  # Keep only significant KEGG pathways
  dplyr::filter(qvalue <= 0.05) %>% 
  # Keep only KEGG pathways related to lipids metabolism
  dplyr::filter(
    (subcategory == "Lipid metabolism") | 
      (str_detect(str_to_lower(Description), paste(lipid_terms, collapse = "|")))) %>% 
  # Remove some GO BP terms
  dplyr::filter(
    category %in% c("Metabolism", "Organismal Systems"),
    !str_detect(Description, paste(nonlipid_terms, collapse = "|")))


## Save KEGG ORA results
write_tsv(kegg_lipid, file.path(res_dir, "kegg_lipid.txt"))


# Retrieve DE genes in GO terms and KEGG pathways ------------------------------
## GO BP
go_lipid_revigo_imp_2genetable <- go_lipid_revigo_imp_res %>% 
  left_join(ds_meta, by = "dataset") %>% 
  dplyr::select(ID, Description, geneID, status, celltype, Tissue, dataset, `Dataset name`) 

write_tsv(go_lipid_revigo_imp_2genetable, file.path(res_dir, "go_lipid_revigo_imp_degs.txt"))


## KEGG
kegg_lipid_2genetable <- kegg_lipid %>% 
  left_join(ds_meta, by = "dataset") %>% 
  dplyr::select(ID, Description, geneID, status, celltype, Tissue, dataset, `Dataset name`) %>% 
  separate_rows(geneID, sep = "/")

write_tsv(kegg_lipid_2genetable, file.path(res_dir, "kegg_lipid_degs.txt"))


