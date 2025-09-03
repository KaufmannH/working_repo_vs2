# Are hvgs enriched in innate immune response genes in innate immune cells?

# 1. Data preparation
# 2. Hypergeometric test
# 3. Reporting

library(readxl)
library(tidyverse)

palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")


# 1. Data preparation

# get immune cell clusters 
cell_type_list_selected <- read_xlsx('droplet_analysis/data/cell_type_list_selected.xlsx')
innate_cell_list <- as.vector(na.omit(cell_type_list_selected$innate_cell))

# load dataset
df_main_all <- read.csv("droplet_analysis/data/combined_data.csv")

#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('mouse_innate_genes.csv')
mouse_innate_genes_raw <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)
length(mouse_innate_genes_raw)

# get only immune cells, tag innate immune genes
df_main <- df_main_all  %>% 
    filter(manual_final %in% innate_cell_list) # %>% 
    #mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes_raw, collapse = "|")))) %>% 
    #relocate(gene, inner_gene, .after = 4)

# get list of all genes
all_genes_list <- df_main %>%
    distinct(gene) %>%
    pull(gene) 

# get mouse innate genes in dataset
mouse_innate_genes <- intersect(mouse_innate_genes_raw, all_genes_list)

# list of hvgs
hvg_list <- df_main  %>% 
    filter(hvg == TRUE) %>% 
    distinct(gene) %>% 
    pull(gene)
print((hvg_list))

# list of overlap between hvg and innate genes
overlap_list <- intersect(mouse_innate_genes, hvg_list)
overlap_list


# signature (overall pool) - category of interest = failures
n <- length(all_genes_list) - length(hvg_list)  # 10375
# category of interest
m <- length(hvg_list) # 8038
# number of draws
k <- length(mouse_innate_genes) #606
# overlap 
q <- 0:(length(overlap_list) - 1) # 434
q2 <- length(overlap_list) 


# 2. Hypergeometric test

# probability phyper
# P(X > q) 
prob_number <- phyper(q2, m, n, k, lower.tail = FALSE)
prob_number <- formatC(prob_number, format = "e", digits = 2)
prob_number
prob <- phyper(q, m, n, k, lower.tail = FALSE)
data_p <- data.frame(q = q, cum_prob = prob)
paste(m, n, k, q)



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

ggsave("droplet_analysis/plots/all_ages/phyper.png", plot = plot_p, width = 8, height = 6)



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

ggsave("droplet_analysis/plots/all_ages/dhyper.png", plot = plot_d, width = 14, height = 6)


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
    x = "Number of Clusters containing the HVG",
    y = "Number of HVGs"
  )
 ggsave("droplet_analysis/plots/all_ages/hist_hvg.png", plot = histogram, width = 8, height = 6)








