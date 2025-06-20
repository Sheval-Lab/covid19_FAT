# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(ragg)


# Set paths --------------------------------------------------------------------
data_dir <- file.path("results", "scRNAseq", "02_ORA")
res_dir <- file.path("results", "scRNAseq", "03_plots")


# Set ggplot2 theme ------------------------------------------------------------
theme_set(
  theme_bw() +
    theme(
      legend.justification = c("right"),
      strip.text = element_text(size = rel(1.1)),
      axis.text = element_text(color = "black", size = rel(0.9)),
      axis.text.x = element_text(angle = 35, hjust = 1)))


# Input data -------------------------------------------------------------------
## Datasets meta
ds_meta <- read_tsv("sc_dataset_meta.txt")

## GO ORA results: redundant terms are substituted with representative ones
go_lipid_revigo_imp_res <- read_tsv(file.path(data_dir, "go_lipid_revigo_imp_res.txt"))

## Combine ORA results with datasets metadata
go_lipid_revigo_imp_2plot <- go_lipid_revigo_imp_res %>% 
  left_join(ds_meta, by = "dataset") %>% 
  # Add labels for Up/Down-regulated GO terms
  mutate(
    symbol = if_else(status == "Up", "\u25D7", "\u25D6"),
    color = if_else(status == "Up", "#f4a261", "#0a9396")) 


# Semicircle plots -------------------------------------------------------------
## By tissue -------------------------------------------------------------------
tissues <- c("Lung", "Airway", "Nose", "Liver", "Blood", "Lymph node", "Brain", "Heart", "Pancreas", "Spleen", "Kidney", "Urine", "Placenta")

go_lipid_revigo_imp_2plot_by_tissue <- go_lipid_revigo_imp_2plot %>% 
  distinct(Description_revigo, dataset, status, Tissue, symbol, color) %>% 
  # At least one cell type in the dataset has Up/Down-regulated GO term
  count(Description_revigo, status, Tissue, symbol, color) %>% 
  # distinct(Description_revigo, status, Tissue, symbol, color, n) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description_revigo = fct_rev(as.factor(Description_revigo)))

min_size_t <- min(go_lipid_revigo_imp_2plot_by_tissue$n)
max_size_t <- max(go_lipid_revigo_imp_2plot_by_tissue$n)


go_lipid_revigo_imp_2plot_by_tissue %>%
  ggplot(aes(x = Tissue, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ndatasets",
    breaks = log2(c(min_size_t, c(min_size_t:max_size_t)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_t, c(min_size_t:max_size_t)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 
  
ggsave(file.path(res_dir, "GOBP_by_tissue_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 12, height = 10, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_tissue_lipid_keywords.HC.svg"), width = 12, height = 10, units = "cm", scale = 1/0.6)


## By dataset - only Lung and Airway -------------------------------------------
datasets_order <- c("Delorey TM (Lung)", "Izar B (Lung)", "Xu (Lung)", "Kropski JA (Lung)", "Bost P (BALF)", "Liao (BALF)", "Delorey TM (Trachea)", "Eddins (Endotracheal aspirates)", "Ravindra (Bronchial epithelium)", "Misharin (Bronchial epithelium)")

go_lipid_revigo_imp_2plot_by_dataset <- go_lipid_revigo_imp_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway")) %>% 
  mutate(
    `Dataset name` = str_remove(`Dataset name`, "et al. "),
    `Dataset name` = factor(`Dataset name`, levels = datasets_order)) %>% 
  distinct(Description_revigo, dataset, `Dataset name`, celltype, status, Tissue, symbol, color) %>% 
  count(Description_revigo, dataset, `Dataset name`, status, Tissue, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description_revigo = fct_reorder(Description_revigo, n))

min_size_ds <- min(go_lipid_revigo_imp_2plot_by_dataset$n)
max_size_ds <- max(go_lipid_revigo_imp_2plot_by_dataset$n)


go_lipid_revigo_imp_2plot_by_dataset %>% 
  ggplot(aes(x = `Dataset name`, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  facet_wrap(vars(Tissue), scales = "free_x") +
  ggh4x::force_panelsizes(cols = c(6, 4)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ncell types",
    breaks = log2(c(min_size_ds, c(min_size_ds:max_size_ds)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_ds, c(min_size_ds:max_size_ds)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 

ggsave(file.path(res_dir, "GOBP_by_dataset_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 14, height = 10, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_dataset_lipid_keywords.HC.svg"), width = 14, height = 10, units = "cm", scale = 1/0.6)


## By cell type - only Lung and Airway ----------------------------------------
go_lipid_revigo_imp_2plot_by_celltype <- go_lipid_revigo_imp_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway")) %>% 
  mutate(celltype_tidy = str_remove_all(celltype, "^[A-Z0-9]*_|_[A-Z0-9]*$")) %>% 
  distinct(Description_revigo, dataset, celltype_tidy, status, Tissue, symbol, color) %>% 
  count(Description_revigo, celltype_tidy, status, Tissue, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description_revigo = str_wrap(Description_revigo, width = 40),
    Description_revigo = fct_reorder(Description_revigo, n))

min_size_ct <- min(go_lipid_revigo_imp_2plot_by_celltype$n)
max_size_ct <- max(go_lipid_revigo_imp_2plot_by_celltype$n)


go_lipid_revigo_imp_2plot_by_celltype %>% 
  ggplot(aes(x = celltype_tidy, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  facet_wrap(vars(Tissue), scales = "free_x") +
  ggh4x::force_panelsizes(cols = c(3, 2)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ndatasets",
    breaks = log2(c(min_size_ct, c(min_size_ct:max_size_ct)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_ct, c(min_size_ct:max_size_ct)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 

ggsave(file.path(res_dir, "GOBP_by_celltype_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 20, height = 10, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_celltype_lipid_keywords.HC.svg"), width = 20, height = 10, units = "cm", scale = 1/0.6)

