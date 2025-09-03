# Are hvgs enriched in innate immune response genes in innate immune cells?

library(Seurat)
library(readxl)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(clusterProfiler)
library(ggplot2)



# get immune cell clusters 
cell_type_list_selected <- read_xlsx('cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))


#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

expression_matrix_all <- read_excel("combined_data.xlsx")
# list of all genes in TMS
all_genes_list <- expression_matrix_all  %>% 
    distinct(gene) %>% 
    pull(gene)

expression_matrix <- expression_matrix  %>% 
    filter(manual_final %in% innate_cell_list) %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    relocate(gene, .after = 3) 
print((expression_matrix)[50:60,])




# list of hvgs
hvg_list <- expression_matrix  %>% 
    filter(hvg == TRUE) %>% 
    distinct() %>% 
    pull(gene)
print((hvg_list))


# df for hvgs: either innate response gene or not 
hvg_df <- expression_matrix  %>% 
    mutate(selection = ifelse(inner_gene, "innate response gene", "other")) %>% 
    filter(hvg == TRUE) %>% 
    dplyr::select(selection, gene)
print((hvg_df))


# list for innate innune genes
innate_list <- expression_matrix %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    filter(inner_gene == TRUE)  %>% 
    pull(gene)
print(innate_list)



#  shows which metabolic or signaling pathways are potentially impacted by the conditions tested
# overrepresentation analysis (hypergeometric test)
result <- enricher(
  gene = hvg_list, # list of overexpressed genes
  universe = all_genes_list, # use tms as background (all that could be measured)
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = hvg_df
)

# View results
print(result)

# No, there is no enrichment.


