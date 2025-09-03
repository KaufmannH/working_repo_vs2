# From A5 we got a list of lowly expressed genes (lower quartile) and low variability and the ones involved in immune response are selected. 
# Now we want to know what pathways they are enriched in.

# replaced by A9

library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ggplot2)


c("Tnf", "Il6", "Ccl2", "Cxcl10", "Mapk1", "Pik3ca", "Cd40", "Stat1", "Nfkb1")
gene_df <- read.csv('data/low_low_gene_top_100.csv')
gene_list <- gene_df$gene
head(gene_list)

gene_entrez <- mapIds(org.Mm.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  # Remove missing values

# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez, organism = "mmu") 
dotplot_kegg <- dotplot(kegg_res, showCategory = 10)
ggsave('plots/3_m/dotplot_kegg.png', plot = dotplot_kegg)


#GO Enrichment
go_res <- enrichGO(gene = gene_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")
dotplot_go <- dotplot(go_res, showCategory = 10)
ggsave('plots/3_m/dotplot_go.png', plot = dotplot_go)
