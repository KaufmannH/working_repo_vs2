# First broad analysis: 

# in progress

# Are HVGs expressed ingreasingly over time?



library(readxl)
#library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(purrr)
#library(patchwork)
#library(ggridges)

# over time. for each age find 

# get immune cell clusters (manually selected)
cell_type_list_selected <- read_xlsx('droplet_analysis/data/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

df_main_all <- read.csv("droplet_analysis/data/combined_data.csv")

# list of all genes in TMS
all_genes_list <- df_main_all  %>% 
    distinct(gene) %>% 
    pull(gene)

df_main <- df_main_all  %>% 
    filter(manual_final %in% innate_cell_list) %>% 
    mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    relocate(inner_gene, .after = 4) %>% 
    relocate(gene, .after = 3) 
print((df_main)[50:60,])


# df for hvgs: either innate response gene or not 
hvg_df <- df_main  %>% 
    #mutate(selection = ifelse(inner_gene, "innate response gene", "other")) %>% 
    #filter(cluster_id == 'Bladder_18m_male_13')  
    group_by(hvg) %>%
    mutate(mean_hvg_expression = mean(gmean))


head(hvg_df$mean_hvg_expression)



# plot

# Create the ridge plot

ggplot(hvg_df, aes(x = df_main$gmean, y = as.factor(df_main$cluster_id), fill = df_main$age)) +
  geom_density_ridges(scale = 2, alpha = 0.8) +
  facet_wrap(~ df_main$gene, scales = "free") +  # Facet by gene, with independent scales
  labs(x = "Expression Level", y = "Age Group") +
  theme_ridges() +
  theme(legend.position = "none")



# violin plot 1
hvg_df$cluster_hvg <- interaction(hvg_df$cluster_id, hvg_df$hvg)

violin_hvg <- ggplot(hvg_df, aes(x = cluster_hvg, y = gmean, fill = hvg)) +
  geom_violin(scale = "width", alpha = 0.8) +
  facet_wrap(~ cluster_id, nrow = 8) + 
  labs(x = "Cluster & HVG Status", y = "Expression Level") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis_d()
ggsave("droplet_analysis/plots/violin_hvg.png", violin_hvg)



# violin 2


# number of tissues
hvg_list_tissue <- hvg_df %>%
  distinct(tissue) %>%
  pull()
hvg_list_tissue



plot_violin <- function(df){

  violin_hvg <- ggplot(hvg_df, aes(x = cluster_hvg, y = mean_hvg_expression, fill = hvg)) +
    geom_violin(scale = "width", alpha = 0.8) +
    labs(
      x = "Cluster & HVG Status",
      y = "Expression Level"
    ) +
    theme_classic() +
    theme(
      legend.position = "top",  
      strip.text = element_text(face = "bold"), 
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D")) 
}

result = map2(hvg_list_tissue, hvg_df, ~ plot_violin(.y, title = .x))

ggsave("droplet_analysis/plots/violin_hvg.png", violin_hvg)


top_genes <- df_main %>%
  group_by(gene, cluster_id) %>%
  summarize(mean_expression = mean(gmean, na.rm = TRUE)) %>%
  arrange(desc(mean_expression)) %>%
  slice_head(n = 20)  

# Plot the top genes
ggplot(top_genes, aes(x = reorder(gene, mean_expression), y = mean_expression)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Gene", y = "Mean Expression Level", title = "Top 20 Most Expressed Genes") +
  theme_minimal()




# Calculate mean expression for each gene in each cluster
top_genes_per_cluster <- df_main %>%
  group_by(cluster_id, gene) %>%
  summarize(mean_expression = mean(gmean, na.rm = TRUE)) %>%
  arrange(cluster_id, desc(mean_expression)) %>%
  group_by(cluster_id) %>%
  slice_head(n = 5)  

# Join the top genes back with the original data
df_top_genes <- df_main %>%
  semi_join(top_genes_per_cluster, by = c("cluster_id", "gene"))

# Create the violin plot
ggplot(df_top_genes, aes(x = gene, y = gmean)) +
  geom_violin(scale = "width", alpha = 0.8) +
  facet_wrap(~ cluster_id, scales = "free", nrow = 6) +  
  labs(x = "Top Genes per Cluster", y = "Expression Level") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() 

print(df_top_genes)