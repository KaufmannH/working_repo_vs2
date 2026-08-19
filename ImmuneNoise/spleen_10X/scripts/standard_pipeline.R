
library(tidyverse)
library(Seurat)
library(ggplot2)

# 2 months male mice (3 individulas i think)
# need old seurat version to make it compatible 

expression_matrix <- ReadMtx(
  mtx = "data/raw/matrix.mtx.gz", 
  features = "data/raw/features.tsv.gz",
  cells = "data/raw/barcodes.tsv.gz")

so <- CreateSeuratObject(counts = expression_matrix)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")


# percent mt genes filtered out
plot <- so@meta.data %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  geom_hline(yintercept = 250) +
  geom_vline(xintercept = 1000) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw()

ggsave('plots/pct_mt_genes.png', plot)


# num genes
plot <- ggplot(so@meta.data, aes(x=nFeature_RNA)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  geom_vline(xintercept = 500) +
  theme_bw()

ggsave('plots/num_genes.png', plot)


# num counts
plot <- ggplot(so@meta.data, aes(x=nCount_RNA)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  geom_vline(xintercept = 2500) +
  theme_bw()

ggsave('plots/num_counts.png', plot)




# qc filtering

# cutoffs like in regnoise pipeline script
so[["qc_pass"]] <- so@meta.data$nFeature_RNA >= 500 & 
  so@meta.data$nCount_RNA >= 2500 &
  so@meta.data$percent.mt <= 5

table(so@meta.data$qc_pass)

so <- subset(so, subset = qc_pass)

saveRDS(so, 'data/so_filtered_stringent.rds')






# normalization (also like in regnoise pipeline) -> not needed for regnoise pipeline tho

so <- SCTransform(so, assay = "RNA", verbose = FALSE)
so <- RunPCA(so, verbose = FALSE, npcs = 100)

plot <- ElbowPlot(object = so, ndims = 100) +theme_classic()
ggsave('plots/elbow_plot_pcs.png', plot)


so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
plot <- DimPlot(object = so, reduction = "umap", pt.size = 0.3) + NoLegend()
ggsave('plots/umap.png', plot)


# clustering
cluster_res <- 1.3
dim_neighbors <- 30
so <- FindNeighbors(so, reduction = "pca", dims = 1:dim_neighbors, verbose = FALSE) %>%
    FindClusters(resolution = cluster_res, algorithm = 4, verbose = FALSE)


plot <- DimPlot(object = so,   reduction = "umap",    pt.size = 0.3,  label = FALSE,   group.by = "seurat_clusters")
ggsave('plots/umap_clusters.png', plot)






# add metadata cols for pipeline
so <- readRDS('data/so_filtered_stringent.rds')
so <- SeuratObject::RenameAssays(so, RNA = "originalexp")

DefaultAssay(object = so) <- "originalexp"

so@meta.data <- so@meta.data |>
  mutate(tissue = 'Spleen',
          age = '2m',
          sex = 'male',
          cell_ontology_class = 'Spleenocytes',
          mouse.id = 'unknown')
head(so@meta.data)

# test
#so_full <- subset(so, subset = nFeature_originalexp > 2 & nCount_originalexp > 2 & nCount_originalexp < 10000 & percent.mt < 90)


saveRDS(so, 'data/prepped_for_pipeline/so_prepped.rds')
