
library(Seurat)
library(dplyr)

classic_df <- readRDS("data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.rds")
immune_df <- readRDS("data/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.rds")

meta_data_classic <- classic_df@meta.data
head(meta_data_classic)

unique(meta_data_classic$Age_bin)
unique(meta_data_classic$Tissue)



meta_data_immune <- immune_df@meta.data
head(meta_data_immune)

unique(meta_data_immune$Age_bin)
unique(meta_data_immune$Tissue)
unique(meta_data_immune$granular)

no_macs <- meta_data_immune |>
    filter(granular == c("Immune (macrophage I)", "Immune (macrophage II)", "Immune (macrophage III)")) |>
    group_by(tissue) |>
    count()
no_macs

dim(meta_data_immune)

