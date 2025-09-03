# Are innate immune genes enriched in LVGs in innate immune cells?
# lung, 3m, only monocytes. macrophages, DCs

# wrong

# 1. Data loading
# 2. Data preparation
# 3. Hypergeometric test
# 4. Reporting

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(ComplexHeatmap)


# colors
palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")



## 1. Data preparation

# get immune cell clusters (manually selected)
cell_type_list_selected <- read_xlsx('data/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))
innate_cell_list

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('data/mouse_innate_genes.csv')
mouse_innate_genes_raw <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

# load lineage specific genes
file_paths <- list.files(path = "data/lineage_spec_markers", pattern = "\\.tsv$", full.names = TRUE)

data_list <- lapply(file_paths, function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  cell_type <- sub("\\.tsv$", "", basename(file))  
  df$cell_type <- cell_type 
  return(df)})
markers_bound <- bind_rows(data_list) 
head(markers_bound)
markers_df <- markers_bound %>%
  mutate(official.gene.symbol = str_to_lower(official.gene.symbol),           
         official.gene.symbol = str_replace(official.gene.symbol, "^(.)", str_to_upper)) 
unique(markers_df$organ)

markers_df <- markers_df %>%
    filter(grepl("Mm", species)) %>%
    select(official.gene.symbol, cell_type, organ) %>%
    rename(gene = official.gene.symbol) %>%
    filter(cell_type %in% c("leukocyte_dc", "leukocyte_macrophage", "leukocyte_monocyte")) %>%
    filter(organ %in% c("Lungs")) %>%
    mutate(cell_type = case_when(
      str_detect(cell_type, "leukocyte_dc") ~ "DC",
      str_detect(cell_type, "leukocyte_macrophage") ~ "Macrophage",
      str_detect(cell_type, "leukocyte_monocyte") ~ "Monocyte",
      TRUE ~ as.character(cell_type))) %>%
    select(-organ)
head(markers_df)


df_main_all <- read.csv("data/combined_data.csv")
colnames(df_main_all)

# 
housekeeping_df <- read.csv('data/Housekeeping_TranscriptsMouse.csv', sep = ";") 

head(housekeeping_df)

housekeeping_list <- housekeeping_df %>%
    pull(Genes)
housekeeping_list


# filter for age and cell types
df_main_filtered <- df_main_all  %>% 
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(age = as.numeric(str_remove(age, "m"))) %>% 
  filter(age == 3) %>%
  filter(tissue == "Lung") %>%
  filter(manual_final %in% c("leukocyte_dc", "leukocyte_macrophage", "leukocyte_monocyte")) %>%
  mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
  #filter(inner_gene == TRUE) %>%
  mutate(cell_type = manual_final) %>% 
  mutate(cell_type = case_when(
      str_detect(cell_type, "leukocyte_dc") ~ "DC",
      str_detect(cell_type, "leukocyte_macrophage") ~ "Macrophage",
      str_detect(cell_type, "leukocyte_monocyte") ~ "Monocyte",
      TRUE ~ as.character(cell_type))) %>%
  select(cluster_id, cell_type, gene, res_var, gmean, hvg, lvg, num_cells_per_cluster, inner_gene)
head(df_main_filtered)
dim(df_main_filtered)

# tag genes
df_main <- df_main_filtered %>%
    # tag house keeping genes
    mutate(housekeeping_gene = (str_detect(gene, paste(housekeeping_list, collapse = "|")))) %>% 
    mutate(mouse_id = str_extract(cluster_id, ".*(?=_[^_]+$)")) 
head(df_main)
unique(df_main$cluster_id)


# join dfs to tag markers
df_main_markers <- df_main %>%
  left_join(markers_df %>% mutate(is_marker = TRUE), by = c("gene", "cell_type")) %>%
    mutate(is_marker = ifelse(is.na(is_marker), FALSE, TRUE))
head(df_main_markers)
sum(df_main_markers$is_marker)

marker_summary <- df_main_markers %>% 
  group_by(cell_type, is_marker) %>%
  summarise(count = n(), .groups = 'drop')
marker_summary
  
 df_main_with_markers <- df_main_markers %>%
  mutate(
    gene_set = case_when(
      hvg ~ "HVG",  
      lvg ~ "LVG", 
      TRUE ~ "Other")) %>%
  mutate(
    gene_type = case_when(
      inner_gene ~ "Immune Response Gene", 
      is_marker ~ "Lineage Specific Gene",
      inner_gene & is_marker ~ "Both",
      TRUE ~ "Other")) %>%
  select(gmean, gene, cell_type, gene_set, gene_type, res_var, cluster_id)
head(df_main_with_markers, n = 300)

# overlap between lineage specific an response genes: 0
sum(df_main_with_markers$gene_type == "Both")


test <- df_main_with_markers %>%
  group_by(cluster_id, gene_set) %>%
  summarise(gene_count = n(), .groups = "drop")
test





# 2. Hypergeometric test

all_clusters <- unique(df_main_with_markers$cluster_id)
hypergeometric_test <- function(cluster_id, df, mouse_innate_genes_raw) {
  # Filter data for the specific cluster
  df_cluster <- df[df$cluster_id == cluster_id, ]

  # get list of all genes
  all_genes_list <- df %>%
      distinct(gene) %>%
      pull(gene) 

  # get mouse innate genes in dataset
  mouse_innate_genes <- intersect(mouse_innate_genes_raw, all_genes_list)

  # list of lvgs
  lvg_list <- df  %>% 
      filter(gene_set == "LVG") %>% 
      distinct(gene) %>% 
      pull(gene)

  # list of overlap between hvg and innate genes
  overlap_list <- intersect(mouse_innate_genes, lvg_list)


  # signature (overall pool) - category of interest = failures
  n <- length(all_genes_list) - length(lvg_list) # 34
  # category of interest
  m <- length(lvg_list) # 14 602
  # number of draws (536?)
  k <- length(mouse_innate_genes)  # 542
  # overlap 
  q <- 0:(length(overlap_list) - 1) 
  q2 <- length(overlap_list) # 537

  p_value <- phyper(q2 - 1, m, n, k, lower.tail = FALSE)
  
  return(p_value)
}


p_values <- sapply(all_clusters, function(cluster_id) {
  hypergeometric_test(cluster_id, df_main_with_markers, mouse_innate_genes_raw)
})
p_values

p_values_df <- data.frame(
  cluster_id = all_clusters,
  p_value = p_values,
  significant = p_value < 0.05
)
p_values_df






# probability phyper
# P(X > q) 
prob_number <- phyper(q2, m, n, k, lower.tail = FALSE)
prob_number <- formatC(prob_number, format = "e", digits = 2)
prob_number
p rob <- phyper(q, m, n, k, lower.tail = FALSE)
p_value <- phyper(q2, m, n, k, lower.tail = FALSE)
data_p <- data.frame(q = q, cum_prob = prob)
paste(m, n, k, q)

p_value # 0.9921773



# 3. Reporting 
plot_p <- ggplot(data_p, aes(x = q, y = prob)) +
  geom_line(color = "#D6EFC3") +
  geom_point(color = "#D6EFC3") +
  xlab("Number of HVGs that are innate immune genes") +
  ylab("Cumulative Probability") +
  xlim(0, 500) +
  theme_classic() 

observed_q <- length(overlap_list)  
plot_p <- plot_p +
  geom_vline(xintercept = observed_q, linetype = "dashed", color = "#3B052C") +
  annotate("text", x = observed_q, y = 0.5, label = paste0("q = ", observed_q, "\nP = ", prob_number), color = "#3B052C", hjust = -0.1)

ggsave("plots/3_m/phyper_lvg.png", plot = plot_p, width = 8, height = 6)



# density dhyper
densities <- dhyper(q, m, n, k) 
data_d <- data.frame(q = q, density = densities)

plot_d <- ggplot(data_d, aes(x = q, y = density)) +
  geom_line(color = "#D6EFC3") +
  geom_point(color = "#D6EFC3") +
  xlab("Number of HVGs that are innate immune genes") +
  ylab("Probability Density") +
  xlim(0, 450) +
  ylim(0, 0.04) +
  theme_classic()

plot_d <- plot_d +
  geom_vline(xintercept = observed_q, linetype = "dashed", color = "#3B052C") +
  annotate("text", x = observed_q, y = 0.02, label = paste0("q = ", observed_q, "\nP = ", prob_number), color = "#3B052C", hjust = -0.1)

ggsave("plots/3_m/dhyper_lvg.png", plot = plot_d, width = 14, height = 6)


# test is how many clusters each HVG is (could be that the cells are in so many different clusters that the clusters were diluted)

# histogram 
hvg_cluster <- df_main %>%
  filter(hvg == TRUE) %>%
  group_by(gene) %>%
  summarise(number_clusters = n_distinct(cluster_id), .groups = "drop") %>%
  arrange(desc(number_clusters))

hvg_cluster

histogram_data <- hvg_cluster %>%
  group_by(number_clusters) %>%
  summarise(gene_count = n(), .groups = "drop")


histogram <- ggplot(histogram_data, aes(x = number_clusters, y = gene_count)) +
  geom_col(fill = "#70916B", color = "black", size = 0.1) +
  theme_classic()+
  labs(
    x = "Number of Clusters containing the LVG",
    y = "Number of LVGs"
  )
 ggsave("plots/3_m/hist_lvg.png", plot = histogram, width = 8, height = 6)








