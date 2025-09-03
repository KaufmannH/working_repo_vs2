# Are innate immune response genes expressed at a lower level than other lineage specific genes? 

library(Seurat)
library(SeuratData)
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggridges)
library(clusterProfiler)



# get immune cell clusters 
cell_type_list_selected <- read_xlsx('facs_analysis/cell_type_list_selected1.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))
print(innate_cell_list)

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

expression_matrix_all <- read_csv('facs_analysis/combined_data.csv')
print(head(expression_matrix_all))

# list of all genes in TMS
all_genes_list <- expression_matrix_all %>% 
    distinct(gene) %>% 
    pull(gene)


#TODO:  not filtered for immune clusters
expression_matrix <- expression_matrix_all  %>% 
    filter(cell_ontology_class %in% innate_cell_list) %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    relocate(gene, .after = 3)  %>% 
    relocate(cell_ontology_class, .after = 5) 


print((expression_matrix)[50:60,])
print(unique(expression_matrix$cluster_id))

# where do i select for innate clusters again? 

# general analysis of data
clu = 'Aorta_18m_female_1' # 15


gene_counts_cluster <- expression_matrix_all %>%
  group_by(cluster_id) %>%
  summarize(num_genes = n_distinct(gene))
print(gene_counts_cluster)


cell_counts_cluster <- expression_matrix %>%
  group_by(cluster_id) %>%
  summarize(num_cells = n_distinct(gene))
print(gene_counts_cluster)


#DGE: one vs all approach (no DESEQ2)
gsea_results <- map_dfr(unique(expression_matrix$cluster_id), function(clu) {
  

  gene_clu <- expression_matrix %>%
    filter(cluster_id == clu) %>%
    select(cluster_id, gene, gmean)
  print(length(gene_clu$gene))

  gene_other <- expression_matrix %>%
    filter(cluster_id != clu) %>%
    mutate(cluster_id = "other_clusters") %>%
    select(cluster_id, gene, gmean)
  print(gene_other$gene)


  cells <- bind_rows(gene_clu, gene_other)
  #print(colnames(cells))
  #print(intersect(gene_clu$gene, gene_other$gene))

  logfc_data <- gene_clu %>%
    full_join(gene_other, by = "gene") %>%
    mutate(logFC = log2(gmean.x + 1) - log2(gmean.y + 1)) 
 
   print(head(sum(!is.na(logfc_data$gmean.x) & !is.na(logfc_data$gmean.y))))
  # no overlap between the genes of the clusters

  geneList <- logfc_data$logFC
  names(geneList) <- logfc_data$gene
  geneList <- sort(geneList, decreasing = TRUE)

  gsea <- GSEA(geneList, TERM2GENE = cells)
  as.data.frame(gsea) %>%
    mutate(cluster_id = clu)
})

print(head(gsea_results))




# GSEA
# only gsea version
gsea <- expression_matrix %>%
  group_by(cluster_id) %>%
  summarise(
    geneList = list(sort(setNames(as.numeric(gmean), as.character(gene)), decreasing = TRUE))
  ) %>%
  mutate(
    gsea = (map(geneList, ~ as.data.frame(GSEA(.x, TERM2GENE = cells))))) %>%
    dplyr::select(cluster_id, gsea) %>%
    unnest(gsea) 

print(gsea$geneList)
print(head(gsea))

# there are ties in the matrix meaning that there are genes that have the same expression and are therefore ranked the same.
# this can lead to varying results each time i run the analysis

# there are no significanlty overrepresented genes in 


