

library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(tibble)
library(grid)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)



df_main_all <- readRDS("facs_analysis/data/df_with_markers.rds")
colnames(df_main_all)



df_main <- df_main_filtered %>%
    group_by(cluster_id) %>%
    mutate(expression_bin = as.factor(ntile(gmean, 6))) %>%
    rename(variability_class = gene_set) %>%
  mutate(category = paste(variability_class, expression_bin, sep = "_"))
head(df_main)


lps_lvg1 <- df_main |>
    filter(lps_stim == TRUE) |>
    filter(category == "LVG_1") |>
    pull(gene) |>
    unique()
length(lps_lvg1)


lps_lvg_other <- df_main |>
    filter(lps_stim == TRUE) |>
    filter(category %in% c("LVG_2", "LVG_3", "LVG_4", "LVG_5", "LVG_6")) |>
    pull(gene) |>
    unique()
length(lps_lvg_other)


# enrichment LPS LVG 1

gene_entrez <- mapIds(org.Mm.eg.db, keys = lps_lvg1, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  

#GO Enrichment
go_res <- enrichGO(gene = gene_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")
dotplot_go <- dotplot(go_res, showCategory = 10)
ggsave('facs_analysis/plots/Test_9/go_lps_lvg1.png', plot = dotplot_go)



# enrichment LPS LVG rest 

gene_entrez <- mapIds(org.Mm.eg.db, keys = lps_lvg_other, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  

#GO Enrichment
go_res <- enrichGO(gene = gene_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")
dotplot_go <- dotplot(go_res, showCategory = 10)
ggsave('facs_analysis/plots/Test_9/go_lps_lvg_rest.png', plot = dotplot_go)
