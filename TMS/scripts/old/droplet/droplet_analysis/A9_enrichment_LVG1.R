# taking df from A8 (TMS, 3m, lung, macs, mos, dcs) isolated the LVG lowest expression bin


library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ggplot2)


# clusters in lung with high var in LVG 1 : strat_LVG1.csv
# clusters in all tissues with high or low var in LVG 1: 

gene_df <- read.csv('data/strat_LVG1.csv')
gene_list <- gene_df$gene
head(gene_list)

gene_entrez <- mapIds(org.Mm.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  

# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez, organism = "mmu") 
dotplot_kegg <- dotplot(kegg_res, showCategory = 10)
ggsave('plots/3_m/Test_9/dotplot_kegg.png', plot = dotplot_kegg)


#GO Enrichment
go_res <- enrichGO(gene = gene_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")
dotplot_go <- dotplot(go_res, showCategory = 10)
ggsave('plots/3_m/Test_9/dotplot_go.png', plot = dotplot_go)




# what is the difference in gene composition in the LVG 1 cluster between the variable clusters and the non variable clusters? 

# df taken from stratification
df_main <- read.csv('data/strat_LVG1_high_and_low.csv')
high_var_cluster <- c("Lung_3m_male_14", "Lung_3m_male_12", "Lung_3m_male_17", "Kidney_3m_female_7", "Kidney_3m_male_11", "Marrow_3m_female_15","Marrow_3m_male_8")


# classify the genes for high var cluster or low var cluster
high_var_LVG_1_clusters <- df_main |>
  filter(cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 
unique(high_var_LVG_1_clusters$cluster_id)
high_var_LVG_1_gene_list <- unique(pull(high_var_LVG_1_clusters))

low_var_LVG_1_clusters <- df_main |>
  filter(!cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 
unique(low_var_LVG_1_clusters$cluster_id)
low_var_LVG_1_gene_list <- unique(pull(low_var_LVG_1_clusters))

special_high_var_genes <- setdiff(high_var_LVG_1_gene_list, low_var_LVG_1_gene_list) # the ones specific to high
special_low_var_genes <- setdiff(low_var_LVG_1_gene_list, high_var_LVG_1_gene_list) # the ones specific to low
special_low_var_genes
# there are no exclusive genes


## Enrichment for all organs
gene_list <- df_main$gene

gene_entrez <- mapIds(org.Mm.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  

# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez, organism = "mmu") 
dotplot_kegg <- dotplot(kegg_res, showCategory = 10)
ggsave('plots/3_m/Test_9/dotplot_kegg_all_organs.png', plot = dotplot_kegg, width = 6, height = 4, dpi = 300)


#GO Enrichment
go_res <- enrichGO(gene = gene_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")
dotplot_go <- dotplot(go_res, showCategory = 10)
ggsave('plots/3_m/Test_9/dotplot_go_all_organs.png', plot = dotplot_go)
