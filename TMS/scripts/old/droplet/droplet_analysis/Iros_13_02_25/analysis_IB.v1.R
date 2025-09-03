
library(tidyverse)
library(Seurat)


#~~~~~~~~~~~~~~~~~
#Prepare Gene Sets
##################

##Note: in InnateDB_cytokine_signalling.csv some symbols do not seem mouse gene symbol (human, or protein id)

gsets <- list()
gsets[["mm"]] <- list()

d <- read_tsv("sets/Bhatt_2012_lps_data.txt")
gsets[["mm"]][["LPS_up_early"]] <- d %>% dplyr::filter(Cluster <= 5) %>% pull(Gene) %>% unique() %>% sort()
gsets[["mm"]][["LPS_up_late"]] <- d %>% dplyr::filter(Cluster > 5) %>% pull(Gene) %>% unique() %>% sort()

d <- read_delim("sets/Housekeeping_TranscriptsMouse.csv", delim=";")
gsets[["mm"]][["HK"]] <- d %>% pull(Genes) %>% unique() %>% sort()

d <- read_delim("sets/Housekeeping_MostStable_Mouse.csv", delim=";")
gsets[["mm"]][["HK_stable"]] <- d %>% pull(Gene.name) %>% unique() %>% sort()

d <- read_delim("sets/InnateDB_cytokine_signalling.csv", delim=",")
gsets[["mm"]][["Cytokine_signaling"]] <- d %>% pull(name) %>% unique() %>% sort()

d <- read_delim("sets/InnateDB_genes_TLR.csv", delim=",")
gsets[["mm"]][["TLR_signaling"]] <- d %>% pull(name) %>% sort() %>% unique()

d <- read_delim("sets/mouse_innate_genes.csv", delim=",")
gsets[["mm"]][["Innate_immunity"]] <- d %>% pull("Gene Symbol") %>% sort() %>% unique()

d <- read_delim("sets/lineage_specific_markers.csv", delim=",")
cts <- d$cell_type %>% unique()
for (ct in cts) {
  id <- paste0("Lineage_", ct)
  gsets[["mm"]][[id]] <- d %>% dplyr::filter(cell_type == ct) %>% pull("gene") %>% sort() %>% unique()
}
#DC unique 
w <- !gsets[["mm"]]$Lineage_DC %in% c(gsets[["mm"]]$Lineage_Macrophage, gsets[["mm"]]$Lineage_Monocyte)
gsets[["mm"]][["Lineage_DC_uniq"]] <- gsets[["mm"]]$Lineage_DC[w]
#Macro unique 
w <- !gsets[["mm"]]$Lineage_Macrophage %in% c(gsets[["mm"]]$Lineage_DC, gsets[["mm"]]$Lineage_Monocyte)
gsets[["mm"]][["Lineage_Macrophage_uniq"]] <- gsets[["mm"]]$Lineage_Macrophage[w]
#Mono unique 
w <- !gsets[["mm"]]$Lineage_Monocyte %in% c(gsets[["mm"]]$Lineage_Macrophage, gsets[["mm"]]$Lineage_DC)
gsets[["mm"]][["Lineage_Monocyte_uniq"]] <- gsets[["mm"]]$Lineage_Monocyte[w]

#Convert to human gene symbols

m2h <- read_tsv("20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets[["hs"]] <- list()
for (id in names(gsets[["mm"]])) {
  genes_hs <- m2h %>% dplyr::filter(gene %in% gsets[["mm"]][[id]]) %>% pull(ortholog) %>% sort() %>% unique()
  gsets[["hs"]][[id]] <- genes_hs
}

# Also add HK from Lin et al.
## https://academic.oup.com/gigascience/article/8/9/giz106/5570567

d <- read_tsv("sets/scHK_human.symbols.txt", col_names = FALSE)
gsets[["hs"]][["HK_Lin"]] <- d %>% pull(1) %>% unique() %>% sort()


#~~~~~~~~~~~~~
#Granja et al.
##############

so <- readRDS("mpm_data/scRNA-Healthy-Hematopoiesis-191120.subset.so.rds")
so <- NormalizeData(so, verbose = TRUE)

#re-order cell types
lvls_ord <- c("01_HSC", "02_Early.Eryth", "03_Late.Eryth",
              "06_CLP.1", "15_CLP.2", "22_CD4.M", "20_CD4.N1", "21_CD4.N2", 
              "19_CD8.N", "23_CD8.EM", "24_CD8.CM", "16_Pre.B", "17_B",
              "18_Plasma", "25_NK", "05_CMP.LMPP", "07_GMP", "08_GMP.Neut",
              "09_pDC", "10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono", 
              "04_Early.Baso", "14_Unk", "26_Unk")
so@meta.data$orig.BioClassification <- factor(so@meta.data$orig.BioClassification, levels = lvls_ord)
Idents(so) <- so@meta.data$orig.BioClassification

#module score for each gene set

for (id in names(gsets[["hs"]])) {
  foi <- list(c(gsets[["hs"]][[id]]))
  so <- AddModuleScore(object = so, features = foi, ctrl = 100, name = id)
  w <- colnames(so@meta.data) == paste0(id, 1)
  colnames(so@meta.data)[w] <- paste0("score.", id)
}

#violin plots per cell type

plts_vln <- list()
for (id in names(gsets[["hs"]])) {
  id_feat <- paste0("score.", id)
  plts_vln[[id]] <- VlnPlot(so, features = id_feat, pt.size = 0) + 
    NoLegend() + 
    coord_flip() + 
    geom_hline(yintercept = 0) +
    ylim(c(-0.5,0.5))
}

pdf("analysis_IB.v1.moduleScore.vlnPlots.pdf", width = 5, height = 6)
for (id in sort(names(gsets[["hs"]]))) {
  plot(plts_vln[[id]])
}
dev.off()

#aggregated expression for each gene set

for (id in names(gsets[["hs"]])) {
  foi <- gsets[["hs"]][[id]]
  w <- rownames(so@assays$RNA@data) %in% foi
  aggr_expr <- colMeans(so@assays$RNA@data[w,])
  col_id <- paste0("aggrExpr.", id)
  so@meta.data[[col_id]] <- aggr_expr
}

plts_vln <- list()
for (id in names(gsets[["hs"]])) {
  id_feat <- paste0("aggrExpr.", id)
  plts_vln[[id]] <- VlnPlot(so, features = id_feat, pt.size = 0) + 
    NoLegend() + 
    coord_flip() + 
    ylim(c(0,0.8))
}

pdf("analysis_IB.v1.aggrExpr.vlnPlots.pdf", width = 5, height = 6)
for (id in sort(names(gsets[["hs"]]))) {
  plot(plts_vln[[id]])
}
dev.off()

