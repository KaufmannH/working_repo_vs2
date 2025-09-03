
# using fused and cell type annotated seurat object and calc module score

library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(readxl)


# add gene sets
gene_set_list <- list()
gene_set_list[["mm"]] <- list()
gene_set_list[["mm"]]

# add lps
lps_data <- read_excel("droplet_analysis/data/Bhatt_2012_data.xlsx")
lps_data <- lps_data |> rename(group = `...4`, Gene = 'Probe...1' )
gene_set_list[["mm"]][["LPS up early"]] <- lps_data |> dplyr::filter(group <= 5) |> pull(Gene) |> unique() |> sort()
gene_set_list[["mm"]][["LPS up late"]] <- lps_data |> dplyr::filter(group > 5) |> pull(Gene) |> unique() |> sort()



# add the spern DNA condensation signature for basal level
sperm_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_sperm_dna_condensation.xlsx')
sperm_list <- sperm_table |> pull(Symbol) |> unique() |> sort()
gene_set_list[["mm"]][["Sperm"]] <- sperm_list


# add the meiosis signature for basal level
meiosis_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_G1_M_transition_meiosis.xlsx')
meiosis_list <- meiosis_table |> pull(Symbol) |> unique() |> sort()
gene_set_list[["mm"]][["Meiosis"]] <- meiosis_list


# add the oocyte maturation signature for basal level
oocyte_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_oocyte_maturation.xlsx')
oocyte_list <- oocyte_table |> pull(Symbol) |> unique() |> sort()
gene_set_list[["mm"]][["Oocyte"]] <- oocyte_list


# load the mouse/human conversions
m2h <- read.delim("droplet_analysis/data/20200307_ensembl/mouse.txt",  header = FALSE, sep = "\t")
colnames(m2h) <- c("mouse", "human")
h2m <- read.delim("droplet_analysis/data/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
colnames(h2m) <- c("human", "mouse")
head(h2m)


# add autimmune pleiotropic loci
autimmune_df <- read_excel('droplet_analysis/data/Marquez_pleiotropic_genes.xlsx')
autoimmune_list <- autimmune_df |> pull(Gene) |> unique() |> sort()
autoimmune_list
# conversion from human to mouse
mouse_auto <- h2m |> filter(human %in% autoimmune_list) |> pull(mouse)
gene_set_list[["mm"]][["Autoiummunity"]] <- mouse_auto





# load seurat object
# all cell types
main_seurat <- readRDS("droplet_analysis/data/fused_seurat_droplet_3m.rds")
# immune cell types
main_seurat <- readRDS("droplet_analysis/data/immune_seurat_droplet_3m.rds")
#


# set the og assay for the moduel score analysis
DefaultAssay(main_seurat) <- "originalexp" # data was already normalized
main_seurat

# rename cell types
cell_type_map <- c(
  "leukocyte_macrophage"    = "Leukocyte macrophage",
  "leukocyte_progenitor"    = "Leukocyte progenitor",
  "lymphocyte_B"            = "Lymphocyte B",
  "lymphocyte_T"            = "Lymphocyte T",
  "leukocyte_monocyte"      = "Leukocyte monocyte",
  "myeloid"                 = "Myeloid",
  "lymphocyte"              = "Lymphocyte",
  "lymphocyte_plasma"       = "Lymphocyte plasma",
  "lymphocyte_NK"           = "Lymphocyte NK",
  "leukocyte_dc"            = "Leukocyte DC",
  "lymphocyte_NKT"          = "Lymphocyte NKT",
  "leukocyte_neutrophil"    = "Leukocyte neutrophil",
  "leukocyte_basophil"      = "Leukocyte basophil",
  "leukocyte_granulocyte"   = "Leukocyte granulocyte",
  "granulocyte_progenitor"  = "Granulocyte progenitor",
  "leukocyte"               = "Leukocyte",
  "lymphocyte_thymocyte"    = "Lymphocyte thymocyte"
)
main_seurat@meta.data$cell_type <- cell_type_map[main_seurat@meta.data$cell_type]
#change idents for plot
Idents(object = main_seurat) <- "cell_type"


# use addmodule score to get the info of gene tag
for (set in names(gene_set_list[["mm"]])) {
  features <- list(c(gene_set_list[["mm"]][[set]]))
  main_seurat <- AddModuleScore(object = main_seurat, features = features, ctrl = 100, name = set)
  score <- colnames(main_seurat@meta.data) == paste0(set, 1)
  colnames(main_seurat@meta.data)[score] <- paste0("Module score ", set)
}

#violin plots per cell type
plts_vln <- list()
for (set in names(gene_set_list[["mm"]])) {
  id_feat <- paste0("Module score ", set)
  plts_vln[[set]] <- VlnPlot(main_seurat, features = id_feat, pt.size = 0) + 
    NoLegend() + 
    coord_flip() + 
    geom_hline(yintercept = 0) +
    ylim(c(-1,2)) +
    theme(
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold", hjust = 0))
}

pdf("droplet_analysis/cell_level/plots/module_score.pdf", width = 5, height = 6)
for (id in sort(names(gene_set_list[["mm"]]))) {
  plot(plts_vln[[id]])
}
dev.off()





# old


# get hvgs and lvgs
df_main <- df_main |> filter(cell_type == "Macrophage")
clusters <- unique(df_main$cluster_id)
clusters

# not for now
hvg_df <- data.frame(cluster = character(), gene = character(), stringsAsFactors = FALSE)
lvg_df <- data.frame(cluster = character(), gene = character(), stringsAsFactors = FALSE)
for (c in clusters) {
  subset_cluster <- df_main[df_main$cluster_id == c, ]
  hvg_genes <- subset_cluster$gene[subset_cluster$gene_set == "HVG"]
  hvg_genes <- unique(hvg_genes)
  if (length(hvg_genes) > 0) {
    hvg_df <- rbind(hvg_df, data.frame(cluster = c, gene = hvg_genes))
  }

  lvg_genes <- subset_cluster$gene[subset_cluster$gene_set == "LVG"]
  lvg_genes <- unique(lvg_genes)
  if (length(lvg_genes) > 0) {
    lvg_df <- rbind(lvg_df, data.frame(cluster = c, gene = lvg_genes))
  }
}
#_____

## from here i need to do hvg tagging 
# get mcas

main_seurat <- subset(main_seurat, subset = cell_ontology_class == "lung macrophage")
head(main_seurat@meta.data)



output_dir <- "droplet_analysis/cell_level/plots/negativ_control"
for (cat in unique(df_agg_macs$category)) {
  df_subset <- df_agg_macs %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cluster_id), y = gmean, fill = as.factor(seurat_clusters))) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.2,  aes(colour = as.factor(cluster_id))) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = paste("Expression in", cat),
         x = "Cluster ID", y = "Gene Expression") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  filename <- paste0(output_dir, "/violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}





# og flow
cleaned_long_df <- long_df %>%
   filter(values == TRUE) 


# just macrophages
unique(cleaned_long_df$cell_type)

df_agg_macs <- cleaned_long_df %>% 
  filter(cell_type == "Macrophage") %>%
  group_by(category, cluster_id, gene) #%>%
  #summarise(gmean = mean(gmean, na.rm = TRUE), .groups = "drop")



output_dir <- "droplet_analysis/cell_level/plots/negative_contol"
for (cat in unique(df_agg_macs$category)) {
  df_subset <- df_agg_macs %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cluster_id), y = gmean, fill = as.factor(cluster_id))) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.2,  aes(colour = as.factor(cluster_id))) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = paste("Expression in", cat),
         x = "Cluster ID", y = "Gene Expression") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  filename <- paste0(output_dir, "/violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}
