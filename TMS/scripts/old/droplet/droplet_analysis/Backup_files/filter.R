library(Seurat)
library(readxl)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
#library(clusterProfiler)
#library(org.Mm.eg.db)
#library(fgsea)
library(clusterProfiler)
#library(enrichplot)
library(ggplot2)



droplet <- readRDS('SeuratProject.Rds')
print(tail(droplet@meta.data))


# get immune cell clusters 
cell_type_list_selected <- read_xlsx('cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))
innate_clusters <- droplet@meta.data %>%
    filter(manual_final %in% innate_cell_list) %>%#  c("macrophage","monocyte" ,"neutrophil", "dendritic", "killer", "eosinophil", "basophil", "mast", "innate")
    #filter(str_detect(manual_final, paste(c("macrophage","monocyte", "neutrophil", "dendritic", "killer", "eosinophil", "basophil", "mast"), collapse = "|"))) %>%
    rownames()
print(innate_clusters)

droplet_innate <- subset(droplet, cells = innate_clusters)
print(head(droplet_innate))

# filter for the hvg clusters
#hvg_clusters_df <- droplet_innate@meta.data %>%
  
  
  # filter(perc.hvg > 0.9) %>%  #filter(hvg == TRUE) 
  # rownames()
   

#hvg_droplet_innate <- subset(droplet, cells = hvg_clusters_df)
#print(head(hvg_droplet_innate))


#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)


# add column to mark the genes as innate 
#droplet_innate@metadata <- droplet_innate@metadata %>% 
  #  mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|"))))
   # print(sum(gene_metadata$inner_gene))

print((droplet_innate))
print(head(droplet_innate[["RNA"]]$counts))

# UMAP
# relabel the data since it is already QCd, but needed to rescale the data

#droplet_innate[["RNA"]]$data <- droplet_innate[["RNA"]]$counts
#droplet_innate[["RNA"]]$scale.data <- droplet_innate[["RNA"]]$counts
droplet_innate <- FindVariableFeatures(droplet_innate)

droplet_innate <- NormalizeData(droplet_innate)
droplet_innate <- ScaleData(droplet_innate)
droplet_innate <- RunPCA(droplet_innate)
droplet_innate <- FindNeighbors(object = droplet_innate, dims = 1:30)
droplet_innate <- FindClusters(object = droplet_innate)

droplet_innate <- RunUMAP(object = droplet_innate, dims = 1:30)
png("umap_plot.png", width = 1000, height = 1000)
DimPlot(droplet_innate, reduction = "umap")
dev.off()


print(head(droplet_innate@meta.data))


droplet_innate@meta.data[droplet_innate@meta.data$hvg == TRUE, ]



# genes are expressed in multiple clusters: group by cluster
expression_matrix <- read_excel("combined_data.xlsx")
print(head(expression_matrix$gene))

# order df
expression_matrix <- expression_matrix  %>% 
    #mutate(cluster_id = rownames(droplet_innate@meta.data)) %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    #group_by(cluster_id)  %>% 
    #arrange(desc(gmean), .by_group = TRUE)  %>% 
    relocate(gene, .after = 3) 

print((expression_matrix)[50:60,])


# list of all genes in TMS
gene_expression_data <- fread("droplet/full_results_10_22.tsv", header = TRUE, sep = "\t")
all_genes_list <- gene_expression_data  %>% 
    distinct(gene) %>% 
    pull(gene)
print((all_genes_list))


# list of hvgs
hvg_list <- expression_matrix  %>% 
    filter(hvg == TRUE) %>% 
    distinct() %>% 
    pull(gene)
print((hvg_list))


# df for hvgs: either a or b
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




