# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)


# Set paths --------------------------------------------------------------------
data_dir <- file.path("SCovid_DBv2.0", "raw_counts")
res_dir <- file.path("results", "scRNAseq", "01_DGEA")


# Input data -------------------------------------------------------------------
## Raw counts
raw_counts_files <- list.files(data_dir, pattern = "raw_counts", full.names = TRUE)

## Cells metadata
meta_files <- list.files(data_dir, pattern = "metadata", full.names = TRUE)

dataset_ids <- raw_counts_files %>% 
  basename() %>% 
  str_remove("_rawdata_raw_counts.txt.gz")


# DGEA -------------------------------------------------------------------------
run_sc_dgea <- function(ds_id, raw_counts_file, meta_file) {
  print(paste("Start", ds_id))
  
  # Read raw counts
  mtx = read.table(raw_counts_file, row.names = 1, header = TRUE)
  
  # Read cells metadata
  meta = read.csv(meta_file, row.names = 1)
  rownames(meta) = rownames(meta) %>% 
    str_replace_all("-", ".")
  
  # Reorder cells in metadata to match counts matrix colnames
  meta = meta[colnames(mtx),]
  
  # Create Seurat object
  sobj = CreateSeuratObject(
    counts = mtx,
    project = ds_id,
    meta.data = meta)
  
  # Normalize counts
  sobj = NormalizeData(sobj)
  
  # Rename cells based on cell type (referred to as subtype in SCovid DB) and disease status
  sobj$group = paste(sobj$subtype, sobj$Status, sep = ":")
  Idents(sobj) = "group"
  
  # Create table with contrasts (cell type X: COVID-19 vs Healthy)
  groups2compare = tibble(
    ident.1 = paste(unique(sobj$subtype), "COVID-19", sep = ":"),
    ident.2 = paste(unique(sobj$subtype), "Healthy", sep = ":"))
  
  # Run Wilcoxon tests to find differentially expressed genes within each cell type in dataset
  degs_list = pmap(
    groups2compare,
    safely(FindMarkers), # In case of inability to perform a test
    object = sobj,
    test.use = "wilcox")
  
  # Keep dataframes with DGEA results 
  degs_list_res <- degs_list %>% 
    map("result")
  
  # Create a list of cell types used in DGEA
  used_celltypes = str_remove(groups2compare$ident.1, ":.+$")[!map_vec(degs_list_res, is.null)]

  # Add cell type info and combine all results into one table
  degs_df = degs_list_res %>%
    # Remove empty DGEA results (if test produced an error)
    discard(is.null) %>% 
    # Add cell type info
    map2(used_celltypes, \(x, y) mutate(x, celltype = y)) %>% 
    # Move gene name from rownames to a column named "gene"
    map(rownames_to_column, "gene") %>% 
    # Combine all results into one table
    list_rbind()
  
  # Save DGEA results 
  write_tsv(degs_df, file.path(res_dir, paste0(ds_id, ".degs.txt")))
}


## Run DGEA
tibble(
  ds_id = dataset_ids, 
  raw_counts_file = raw_counts_files, 
  meta_file = meta_files) %>% 
  pwalk(run_sc_dgea)
