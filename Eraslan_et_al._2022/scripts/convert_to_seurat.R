library(Seurat)
library(zellkonverter)


sce <- readH5AD("data/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad")
seurat_obj <- as.Seurat(sce, counts = "X", data = "X")
saveRDS(seurat_obj, "data/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.rds")
