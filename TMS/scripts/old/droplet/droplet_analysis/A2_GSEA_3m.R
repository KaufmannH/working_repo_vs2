# Are hvgs expressed at a higher level than other in innate repsonse genes in innate immune cells? 
# only takes 3 months old mice into account -> healty homeostatic conditions

# 1. Data preparation
# 2. GSEA
# 3. Reporting

library(readxl)
library(patchwork)
library(tidyverse)
library(clusterProfiler)


# colors
palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")



## 1. Data preparation

# get immune cell clusters (manually selected)
cell_type_list_selected <- read_xlsx('droplet_analysis/data/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))
innate_cell_list

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

df_main_all <- read_csv("droplet_analysis/data/combined_data.csv")
colnames(df_main_all)

# process dataframe
df_main <- df_main_all  %>% 
    # get only innate immune cells
    filter(manual_final %in% innate_cell_list) %>%  
    # tag innate immune genes
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(gene, inner_gene, manual_final, .after = 4) %>% 
    select (-1) %>%
    filter(inner_gene == TRUE) %>%
    mutate(age = str_extract(age, "\\d+m")) %>%  
    mutate(age = as.numeric(str_remove(age, "m"))) %>% # change if age is already preped
    filter(age == 3)
head(df_main)


# optional: add fake gene to check
#row <- data.frame(cluster_id = "Fat_18m_female_12",   gene = 'test_gene',
#res_var = NA, inner_gene = TRUE, manual_final = 'test', gmean = 350, cluster = NA, hvg = TRUE, mean_resvar_bs = NA,
#perc.hvg = NA, quant_low = NA, quant_high = NA, tissue.x = NA, age = NA, tissue.y = NA, tis_co = NA, tis_fa = NA, 
#fraction.cell_ont = NA, fraction.free_ann = NA,   manual_ann.fa = NA,  manual_ann.co = NA)          
#df_main <- rbind(df_main, row)
#maxgene <- max(df_main$gmean)
#print(df_main[df_main$gmean == maxgene, ])

#optional: select specific clusters for testing
#df_main <- df_main %>%
 # filter(cluster_id %in% c('Fat_18m_female_12', "Marrow_1m_male_14", "Trachea_3m_male_14" ))
#print(unique(df_main$cluster_id))
#print(head(df_main))

# additional control of the data
gene_counts_cluster <- df_main %>%
  group_by(cluster_id) %>%
  summarize(num_genes = n_distinct(gene))
gene_counts_cluster
# the same age, sex and tissue have the same total n of genes bc hvg analysis was run on those. 
# that means if there could be genes listed in the subclusters that were not detected in a cell there

number_clusters <- df_main %>%
  distinct(cluster_id) %>%
  pull()
length(number_clusters)



# 2. GSEA

# make df where gsea results are stored
cluster_reference <- data.frame(cluster_id = unique(df_main$cluster_id), 
                                cell_type = df_main$manual_final[match(unique(df_main$cluster_id), df_main$cluster_id)])
cluster_reference$has_hvg <- FALSE
cluster_reference$p_value <- NA
gsea <- list()


# generate input for gsea   
prepared_data <- df_main %>%
  group_by(cluster_id, tissue) %>%
  summarise(
    TERM2GENE = list(data.frame(
      Term = hvg,
      Gene = gene
    )),
    geneList = list({
      gmeans <- gmean
      names(gmeans) <- gene
      sort(gmeans, decreasing = TRUE)
    }),
    .groups = "drop"
  )
prepared_data

# do gsea 
finisehd_gsea_data <- prepared_data %>%
  mutate(
    gsea_result = map2(geneList, TERM2GENE, ~ GSEA(.x, TERM2GENE = .y))
  )


# get p values and enrichment score 
# data structure: prepared_data$gsea_result[[1]]@result
gsea_summary <- finisehd_gsea_data %>%
  mutate(
    p_values = map(gsea_result, ~ .x@result$pvalue),  
    p.adjust = map(gsea_result, ~ .x@result$p.adjust),  
    enrichment_score = map(gsea_result, ~ .x@result$enrichmentScore),
    NES = map(gsea_result, ~ .x@result$NES)) %>%
  unnest(cols = c(p_values, p.adjust, enrichment_score, NES)) 
gsea_summary

# number of enriched clusters
cluster_enriched_list <- gsea_summary %>%
  distinct(cluster_id) %>%
  pull()
length(cluster_enriched_list)

# TODO: genes that positively contribute to gene list, @result$rank is just number




## 3. Reporting 

# dotplot per tissue
ggplot(gsea_summary, aes(x = tissue, y = -log10(p_values))) +
  geom_point(aes(color = cluster_id), size = 3) +
  labs(
    x = "Tissues",
    y = "-log10(p-Value)") +
  theme_minimal() +
  theme(
    legend.position = "none")
ggsave(filename = paste0("droplet_analysis/plots/3_m/dotplot_tissue.jpeg"), width = 8, height = 6, dpi = 300)



# histogramm of enrichment score
enrich_score <- gsea_summary %>%
  arrange(enrichment_score) 


hist <- ggplot(enrich_score, aes(x = enrichment_score)) +
  geom_histogram(binwidth = 0.01, fill = "#D0BFAD",color = "black", size = 0.2, alpha = 0.7) +
  geom_density(color = palette[(1)], size = 1) +  
  labs(
    x = "Enrichment Score",
    y = "Number of clusters") +
  theme_classic()
ggsave("droplet_analysis/plots/3_m/hist_enrichment_score.png", plot = hist, width = 8, height = 6)

#only the ones that are overexpressed not the underexpressed ones




# for plots: incorporate in main df if cluster is enriched or not
cluster_enriched_df <- df_main %>%
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(age = as.numeric(str_remove(age, "m"))) %>% # change if age is already preped
  select(cluster_id, manual_final, tissue, age) %>%
  mutate(cluster_enriched = cluster_id %in% cluster_enriched_list) %>%
  distinct()
sum(cluster_enriched_df$cluster_enriched)


perc_cell_type <- cluster_enriched_df %>%
  group_by(manual_final) %>%
  reframe(
    tissue = tissue,
    cluster_id = cluster_id,
    age = age,
    cluster_enriched = cluster_enriched,
    num_clusters = n_distinct(cluster_id),                         
    percent = sum(cluster_enriched)/n()*100 )
  
print(perc_cell_type, n= 50)

perc_tissue <- cluster_enriched_df %>%
  group_by(tissue) %>%
  reframe(
    tissue = tissue,
    cluster_id = cluster_id,
    age = age,
    cluster_enriched = cluster_enriched,
    num_clusters = n_distinct(cluster_id),                         
    percent = sum(cluster_enriched)/n()*100 )
  
print(perc_tissue, n= 50)



# bar plot 

# CELL TYPE
# barplot cell type
cell_type_plot <- ggplot(perc_cell_type, aes(x = manual_final, y = percent, fill = as.factor(cluster_enriched))) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_fill_manual(values = palette[c(5, 6)]) +
  labs(
    x = "Cell type",
    y = "Percentage",
    fill = "Is enriched cluster") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.margin = margin(20, 20, 20, 20),  
    axis.title.x = element_text(margin = margin(t = 10))) 


# barplot cell type number of clusters
cell_type_numb_clusters_plot <- ggplot(perc_cell_type, aes(x = manual_final, y = num_clusters)) +
  geom_bar(stat = "identity") +  #  position = "fill" makes it a stacked bar plot
  scale_y_continuous(limits = c(0, 150)) + 
  scale_fill_manual(values = palette[c(2)]) +
  labs(
    y = "Number of clusters per group") +
  theme_classic() +
  theme(
    plot.margin = margin(20, 20, 40, 20),  
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank() , 
    axis.line.x = element_blank())


bar_cell_type_combo_plot <- (cell_type_numb_clusters_plot / cell_type_plot) + plot_layout(heights = c(1, 1))
ggsave("droplet_analysis/plots/3_m/barchart_cell_type_combo.png", plot = bar_cell_type_combo_plot)




# TISSUE
# barplot tissue
tissue_plot <- ggplot(perc_tissue, aes(x = tissue, y = percent, fill = as.factor(cluster_enriched))) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_fill_manual(values = palette[c(4, 5)]) +
  labs(
    x = "Tissue",
    y = "Percentage",
    fill = "Is enriched cluster") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.margin = margin(20, 20, 20, 20),  
    axis.title.x = element_text(margin = margin(t = 10))) 


# barplot tissue number of clusters
tissue_numb_clusters_plot <- ggplot(perc_tissue, aes(x = tissue, y = num_clusters)) +
  geom_bar(stat = "identity") +  #  position = "fill" makes it a stacked bar plot
  scale_y_continuous(limits = c(0, 150)) + 
  scale_fill_manual(values = palette[c(2)]) +
  labs(
    y = "Number of clusters per group") +
  theme_classic() +
  theme(
    plot.margin = margin(20, 20, 40, 20),  
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank() , 
    axis.line.x = element_blank())


bar_tissue_combo_plot <- (tissue_numb_clusters_plot / tissue_plot) + plot_layout(heights = c(1, 1))
ggsave("droplet_analysis/plots/3_m/barchart_tissue_combo.png", plot = bar_tissue_combo_plot)


