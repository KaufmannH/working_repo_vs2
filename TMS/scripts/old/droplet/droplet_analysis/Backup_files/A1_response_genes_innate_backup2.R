# Are hvgs enriched in innate repsonse genes in innate immune cells? 

# hvg are compered per cluster to all other genes in the cluster

library(Seurat)
library(SeuratData)
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggridges)
library(clusterProfiler)
library(httpgd)


# get immune cell clusters 
cell_type_list_selected <- read_xlsx('droplet_analysis/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))


#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

expression_matrix_all <- read_csv("combined_data.csv")
print(expression_matrix_all)

# process dataframe
expression_matrix <- expression_matrix_all  %>% 
    filter(manual_final %in% innate_cell_list) %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    relocate(gene, .after = 3) %>% 
    relocate(manual_final, .after = 5) %>% 
    select (-1) %>%
    filter(inner_gene ==TRUE)
print(unique(expression_matrix$cluster_id))


# add fake gene to check
row <- data.frame(cluster_id = "Fat_18m_female_12",   gene = 'test_gene',
res_var = NA, inner_gene = TRUE, manual_final = 'test', gmean = 350, cluster = NA, hvg = TRUE, mean_resvar_bs = NA,
perc.hvg = NA, quant_low = NA, quant_high = NA, tissue.x = NA, age = NA, tissue.y = NA, tis_co = NA, tis_fa = NA, 
fraction.cell_ont = NA, fraction.free_ann = NA,   manual_ann.fa = NA,  manual_ann.co = NA)          
expression_matrix <- rbind(expression_matrix, row)
maxgene <- max(expression_matrix$gmean)
print(expression_matrix[expression_matrix$gmean == maxgene, ])

#optional: select specific clusters
#expression_matrix <- expression_matrix %>%
 # filter(cluster_id %in% c('Fat_18m_female_12', "Marrow_1m_male_14", "Trachea_3m_male_14" ))
print(unique(expression_matrix$cluster_id))
print(head(expression_matrix))

# additional control of the data
gene_counts_cluster <- expression_matrix %>%
  group_by(cluster_id) %>%
  summarize(num_genes = n_distinct(gene))
print(gene_counts_cluster)


#cluster_list = unique(expression_matrix$cluster_id)
cluster_reference <- data.frame(cluster_id = unique(expression_matrix$cluster_id), 
                                cell_type = expression_matrix$manual_final[match(unique(expression_matrix$cluster_id), expression_matrix$cluster_id)])
cluster_reference$has_hvg <- FALSE
cluster_reference$p_value <- NA
gsea <- list()


#DGE: one vs all approach 
gsea_func <- function(cluster_reference, dataframe) {
  print(paste("staring cluster: ", cluster_reference[2]))
  # select non hvg genes of cluster
  gene_clu <- expression_matrix %>%
    filter(cluster_id == cluster_reference[1]) cluster_reference_summary
    filter(hvg == TRUE) %>%
    mutate(hvg = 'hvg') %>%
    select(cluster_id, hvg, gene, gmean)
  #print(head(gene_clu))
  #n_occur <- data.frame(table(gene_clu$gene))

  # select non hvg genes of cluster
  gene_other <- expression_matrix %>%
    filter(cluster_id == cluster_reference[1]) %>%
    filter(hvg != TRUE) %>%
    mutate(hvg = 'not hvg') %>%
    mutate(cluster_id = "other_clusters") %>%
    select(cluster_id, hvg, gene, gmean) %>%
    #calculate average of all other cluster gmeans
    group_by(gene) %>%
    mutate(gmean = mean(gmean, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE)
  #n_occur <- data.frame(table(gene_other$gene))
  #print(n_occur)

  # bind dfs with non hvgs and hvgs
  bound <- bind_rows(gene_clu, gene_other)
  cells <- bound %>%
    select(hvg, gene) %>%
    rename(Term = hvg, Gene = gene)
 
  # set gene list for GSEA
  geneList <- (bound$gmean)
  names(geneList) <- bound$gene
  geneList <- sort(geneList, decreasing = TRUE)
  #gsea <- GSEA(geneList, TERM2GENE = cells)
  gsea[[as.character(cluster_reference[1])]] <<- GSEA(geneList, TERM2GENE = cells)

}

gsea_results <- apply(cluster_reference, 1, gsea_func) # 2 would be column wise


number_clusters <- seq_along(gsea_results)

#TODO: does not work automatically to plot all clusters together
gsea_plot_func <- function(number_clusters, gsea_results, cluster_reference) {
  
  if (nrow(gsea_results[[number_clusters]]@result) != 0){
    print(gsea_results[[number_clusters]])
    cluster_reference$has_hvg[number_clusters] <<- TRUE
    cluster_reference$p_value[number_clusters] <<- gsea_results[[number_clusters]]@result$pvalue
    #print(cluster_reference$p_value[number_clusters])
    file_name <- paste0("plots/gsea_plot_cluster_", number_clusters, ".jpeg")
    jpeg(file = file_name, width = 800, height = 600)
    gseaplot(gsea_results[[number_clusters]], geneSetID = c("hvg"))
    dev.off()
  } 
}

# manual plots
file_name <- paste0("plots/gsea_plot_cluster_", 20, ".jpeg")
jpeg(file = file_name, width = 800, height = 600)
gseaplot(gsea_results[[20]], geneSetID = c("hvg"))
dev.off()

# clusters with enriched genes are plotted
gsea_plots <- lapply(number_clusters, gsea_plot_func, gsea_results = gsea_results, cluster_reference = cluster_reference)

#print(sum(cluster_reference$has_hvg == TRUE))
print(cluster_reference[cluster_reference['has_hvg'] == TRUE, ])

# Results: 192 clusters have enriched hvgs within the innate response genes

#visualization
cluster_reference_summary <- cluster_reference[!is.na(cluster_reference$p_value), ]
cluster_reference_summary$tissue <- str_extract(cluster_reference_summary$cluster_id, "^[^m_]+")

cluster_reference_summary <- cluster_reference_summary %>%
  mutate(age = str_split_i(cluster_id, "_",2)) 



#unique(cluster_reference_summary$age)

head(cluster_reference_summary)



ggplot(cluster_reference_summary, aes(x = tissue, y = -log10(p_value))) +
  geom_point(aes(color = cluster_id), size = 3) +
  #scale_y_log10()+ 
  labs(
    title = "P-Values by Cell Type for Each Cluster",
    x = "Cell Type",
    y = "-log10(p-Value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )
ggsave(filename = paste0("plots/gsea_plot_cluster.jpeg"), width = 8, height = 6, dpi = 300)


# yes, 12 tissues in immune cells have clusters where innate immune genes are significantly enriched in hvgs
