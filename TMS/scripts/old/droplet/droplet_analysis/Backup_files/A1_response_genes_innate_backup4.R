# Are hvgs overrepresented or overexpressed in innate immune cells? 
# 1. Data preparation
# 2. GSEA
# 3. Reporting

# not used, issues


library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggridges)
library(clusterProfiler)



## 1. data preparation
# get immune cell clusters (manually selected)
cell_type_list_selected <- read_xlsx('droplet_analysis/data/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

df_main_all <- read_csv("droplet_analysis/data/combined_data.csv")
head(df_main_all)

# process dataframe
df_main <- df_main_all  %>% 
    filter(manual_final %in% innate_cell_list) %>% 
    mutate(inner_gene = str_detect(gene, paste(mouse_innate_genes, collapse = "|"))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    relocate(gene, .after = 3) %>% 
    relocate(manual_final, .after = 5) %>% 
    select (-1)
unique(df_main$cluster_id)

# optinal: add fake gene to check 
#row <- data.frame(cluster_id = "Fat_18m_female_12",   gene = 'test_gene',
#res_var = NA, inner_gene = TRUE, manual_final = 'test', gmean = 350, cluster = NA, hvg = TRUE, mean_resvar_bs = NA,
#perc.hvg = NA, quant_low = NA, quant_high = NA, tissue.x = NA, age = NA, tissue.y = NA, tis_co = NA, tis_fa = NA, 
#fraction.cell_ont = NA, fraction.free_ann = NA,   manual_ann.fa = NA,  manual_ann.co = NA)          
#df_main <- rbind(df_main, row)
#maxgene <- max(df_main$gmean)
#print(df_main[df_main$gmean == maxgene, ])

#optional: select specific clusters
#df_main <- df_main %>%
 # filter(cluster_id %in% c('Fat_18m_female_12', "Marrow_1m_male_14", "Trachea_3m_male_14" ))

# additional control 
gene_counts_cluster <- df_main %>%
  group_by(cluster_id) %>%
  summarize(num_genes = n_distinct(gene))
head(gene_counts_cluster)







## 2. GSEA

gsea_func <- function(clusters, dataframe) {
  print(paste("starting cluster: ", clusters))
  cluster_reference <<- rbind(cluster_reference, data.frame(cluster_id = clusters))
  # select non hvg genes of cluster
  gene_clu <- df_main %>%
    filter(cluster_id == clusters) %>%
    filter(hvg == TRUE) %>%
    mutate(hvg = 'hvg') %>%
    select(cluster_id, hvg, gene, gmean)

  # select non hvg genes of cluster
  gene_other <- df_main %>%
    filter(cluster_id == clusters) %>%
    filter(hvg != TRUE) %>%
    mutate(hvg = 'not hvg') %>%
    mutate(cluster_id = "other_clusters") %>%
    select(cluster_id, hvg, gene, gmean) %>%
    #calculate average of all other cluster gmeans
    group_by(gene) %>%
    mutate(gmean = mean(gmean, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE)
  
  # bind dfs with non hvgs and hvgs
  bound <- bind_rows(gene_clu, gene_other)
  cells <- bound %>%
    select(hvg, gene) %>%
    rename(Term = hvg, Gene = gene)
 
  # set gene list for GSEA
  geneList <- (bound$gmean)
  names(geneList) <- bound$gene
  geneList <- sort(geneList, decreasing = TRUE)
  gsea <- GSEA(geneList, TERM2GENE = cells)
}

cluster_list = unique(df_main$cluster_id)
cluster_reference <- data.frame(cluster_id = integer()) # TODO: why integer??
gsea_results <- lapply(cluster_list, gsea_func, dataframe = df_main)
print(gsea_results)
print(cluster_reference)


## 3. reporting 
print((gsea_results[[337]]@result$ID))
cluster_reference$has_hvg <- FALSE
cluster_reference$p_value <- NA
print(cluster_reference$has_hvg[1])
number_clusters <- seq_along(gsea_results)


gsea_plot_func <- function(number_clusters, gsea_results, cluster_reference) {
  if (nrow(gsea_results[[number_clusters]]@result) != 0){
    cluster_reference$has_hvg[number_clusters] <<- TRUE
    cluster_reference$p_value[number_clusters] <<- gsea_results[[number_clusters]]@result$pvalue
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

# only clusters with enriched genes are plotted
gsea_plots <- lapply(number_clusters, gsea_plot_func, gsea_results = gsea_results, cluster_reference = cluster_reference)
print(sum(cluster_reference$has_hvg == TRUE))

# GSEA: would expect them to have hvs enriched. kinda sanity check. 

#Error in object@geneSets[[geneSetID]] : subscript out of bounds
