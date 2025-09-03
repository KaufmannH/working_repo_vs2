# Are innate immune genes enriched in LVGs in innate immune cells?
# lung, 3m, only monocytes. macrophages, DCs

# test for LPS response genes


# 1. Data loading
# 2. Hypergeometric test
# 3. Reporting



library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(purrr)


# colors
palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")


# 1. Data loading
df_main_with_markers_expr_level <- readRDS("data/df_with_markers.rds")

check <- df_main_with_markers_expr_level %>%
  group_by(gene_set) %>%
  summarise(count = n(), .groups = 'drop')
check
colnames(df_main_with_markers_expr_level)


# remove NA cluster
df_main_with_markers_expr_level <- subset(df_main_with_markers_expr_level, !is.na(cluster_id))


lps_response_genes <- df_main_with_markers_expr_level %>% filter(lps_stim == TRUE) %>% distinct(gene)
lps_response_genes

# 2. Hypergeometric test

all_clusters <- unique(df_main_with_markers_expr_level$cluster_id)
hypergeometric_test <- function(cluster_id, df, lps_response_genes) {
  # Filter data for the specific cluster
  df_cluster <- df[df$cluster_id == cluster_id, ]

  # get list of all genes
  all_genes_list <- df %>%
      distinct(gene) %>%
      pull(gene) 

  # get mouse innate genes in dataset
  mouse_innate_genes <- intersect(lps_response_genes, all_genes_list)

  # list of lvgs
  lvg_list <- df_cluster  %>% 
      filter(gene_set == "HVG") %>% 
      distinct(gene) %>% 
      pull(gene)



  # list of overlap between lvg and innate genes
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


p_values_result <- sapply(all_clusters, function(cluster_id) {
  hypergeometric_test(cluster_id, df_main_with_markers_expr_level, mouse_innate_genes_raw)
})
p_values_result

p_values_df <- data.frame(
  cluster_id = all_clusters,
  p_value = p_values_result,
  significant = p_values_result < 0.05
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








