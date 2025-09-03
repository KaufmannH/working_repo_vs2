# taking the TMS df (from P1_assemble_TMS_df.) filter down to 3m lung macs, mos and dcs and statify it


library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# filter and prepare expression df (from P1_assemble_TMS_df.)

df_main_all <- read.csv("data/combined_data.csv") 
colnames(df_main_all)

num_cells <- df_main_all %>%
  group_by(cluster_id, num_cells_per_cluster) %>%
  select (cluster_id, num_cells_per_cluster) %>%
  unique()
tail(num_cells)
sum(num_cells$num_cells_per_cluster, na.rm = TRUE) # one cluster without annotation (non expressed innate genes)

# filter for age and cell types
df_main_filtered <- df_main_all  %>% 
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(age = as.numeric(str_remove(age, "m"))) %>% 
  filter(age == 3) %>%
  #filter(tissue == "Lung") %>% # until here 3460 cells
  filter(manual_final %in% c("leukocyte_dc", "leukocyte_macrophage", "leukocyte_monocyte")) %>%
  mutate(cell_type = manual_final) %>% 
  mutate(cell_type = case_when(
      str_detect(cell_type, "leukocyte_dc") ~ "DC",
      str_detect(cell_type, "leukocyte_macrophage") ~ "Macrophage",
      str_detect(cell_type, "leukocyte_monocyte") ~ "Monocyte",
      TRUE ~ as.character(cell_type))) %>%
  select(tissue, cluster_id, cell_type, gene, res_var, gmean, hvg, lvg, num_cells_per_cluster)
head(df_main_filtered)
sum(df_main_filtered$inner_gene)
dim(df_main_filtered)

unique(df_main_filtered$tissue)


# all genes list
all_genes <- unique(df_main_filtered$gene)

num_cells <- df_main_filtered %>%
  group_by(cluster_id, num_cells_per_cluster) %>%
  select (cluster_id, num_cells_per_cluster) %>%
  unique()
sum(num_cells$num_cells_per_cluster, na.rm = TRUE)

# add the spermatogenesis signature for basal level

sperm_table <- read_excel('data/GO_spermatogenesis_signature.xlsx')
sperm_list <- sperm_table |> pull(Symbol)

df_main_filtered <- df_main_filtered %>%
  mutate(sperm_gene = gene %in% sperm_list) 



# add the meiosis signature for basal level

meiosis_table <- read_excel('data/GO_G1_M_transition_meiosis.xlsx')
meiosis_list <- meiosis_table |> pull(Symbol)

df_main_filtered <- df_main_filtered %>%
  mutate(meiosis_gene = gene %in% meiosis_list) 



#--------
# put the genes according to expression level into 6 bins

df_main <- df_main_filtered %>%
    group_by(cluster_id) %>%
    mutate(expression_bin = as.factor(ntile(gmean, 6))) %>%
    mutate(variability_class = case_when(
    hvg ~ "HVG",
    lvg ~ "LVG",
    TRUE ~ "Other"  )) %>%
  mutate(category = paste(variability_class, expression_bin, sep = "_"))
head(df_main)



# show expression categories (just to check)
summary_table <- df_main |>
  group_by( expression_bin)
    


strat_violin <- ggplot(summary_table, aes(x = expression_bin, y = gmean, fill = expression_bin)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Category",
       y = "Log-normalized expression",
       fill = "Category") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_8/strat_violin_0_all.png', plot = strat_violin,  width = 10, height = 14 )



summary_table <- df_main |>
  group_by(category)

#my_colors <- c("Marker Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
#summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


strat_violin <- ggplot(summary_table, aes(x = category, y = gmean, fill = category)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Category",
       y = "Log-normalized expression",
       fill = "Category") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_8/strat_violin_1_all.png', plot = strat_violin,  width = 10, height = 14 )

# check number of cells in one cluster
c <- df_main_filtered |>
  distinct(cluster_id, num_cells_per_cluster) |>
  filter(num_cells_per_cluster < 100)
c





# get the LPS response genes, replace by proper gene set later on

# load LPS stimulated genes in monocytes
lps_df <- read_excel("data/Bhatt_2012_data.xlsx")
colnames(lps_df)[colnames(lps_df) %in% c("Probe...1", "Probe...2", "Promoter Class", "...4")] <- c("gene", "gene_id", "promoter_class", "group")
lps_df_early <- lps_df %>% filter(group <= 5)
lps_df_late <- lps_df %>% filter(group > 5)
lps_gene_list_early <- lps_df_early %>% pull(gene)
lps_gene_list_late <- lps_df_late %>% pull(gene)


df_main <- df_main %>% 
 mutate(lps_stim_early = str_detect(gene, paste(lps_gene_list_early, collapse = "|"))) %>%
 mutate(lps_stim_late = str_detect(gene, paste(lps_gene_list_late, collapse = "|"))) %>%
 mutate(lps_stim = lps_stim_early | lps_stim_late)  
 
 
df_lps <- df_main %>%
 filter(lps_stim == TRUE)
head(df_lps)


summary_table <- df_main |>
    group_by(category)

#my_colors <- c("Marker Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
#summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


strat_violin <- ggplot(summary_table, aes(x = category, y = gmean, fill = category)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Category",
       y = "Log-normalized expression",
       fill = "Category") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_8/strat_violin_2_all.png', plot = strat_violin,  width = 10, height = 14 )




# how many genes are in each category?

summary_table <- df_main |>
    group_by(cluster_id, category) |>
    summarise(cell_type, gene_count = n_distinct(gene)) 
summary_table

strat_violin <- ggplot(summary_table, aes(x = category, y = gene_count, fill = category)) +
  geom_col() +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Category",
       y = "Number of genes per category",
       fill = "Category") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_8/strat_bar_3_all.png', plot = strat_violin,  width = 10, height = 14 )



# correlation of LVG 1 genes across clusters


df_lvg1 <- df_main %>%
  filter(category == "LVG_1")

cluster_genes <- df_lvg1 %>%
  group_by(cluster_id) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

cluster_names <- cluster_genes$cluster_id
n <- length(cluster_names)
shared_genes_matrix <- matrix(0, nrow = n, ncol = n,
                              dimnames = list(cluster_names, cluster_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    genes_i <- cluster_genes$genes[[i]]
    genes_j <- cluster_genes$genes[[j]]
    shared_genes_matrix[i, j] <- length(intersect(genes_i, genes_j))
  }
}

col_fun <- colorRamp2(c(0, max(shared_genes_matrix)), c("white", "darkred"))

png("plots/3_m/Test_8/lvg1_corr_all.png", width = 1200, height = 1000)

Heatmap(shared_genes_matrix,
        name = "Shared genes",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(shared_genes_matrix[i, j], x, y, gp = gpar(fontsize = 8))
        })

dev.off()

# there is a higher correlation between the tissues and sexes

# is the correlation of lvg 1 genes that are LPSs stronger?


df_lvg1_lps <- df_lps %>%
  filter(category == "LVG_1") 

cluster_genes <- df_lvg1_lps %>%
  group_by(cluster_id) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

cluster_names <- cluster_genes$cluster_id
n <- length(cluster_names)
shared_genes_matrix <- matrix(0, nrow = n, ncol = n,
                              dimnames = list(cluster_names, cluster_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    genes_i <- cluster_genes$genes[[i]]
    genes_j <- cluster_genes$genes[[j]]
    shared_genes_matrix[i, j] <- length(intersect(genes_i, genes_j))
  }
}

col_fun <- colorRamp2(c(0, max(shared_genes_matrix)), c("white", "darkred"))

png("plots/3_m/Test_8/lvg1_corr_lps.png", width = 1200, height = 1000)

Heatmap(shared_genes_matrix,
        name = "Shared genes",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(shared_genes_matrix[i, j], x, y, gp = gpar(fontsize = 8))
        })

dev.off()










# there are clusters with high LVG_1 expression:
high_var_cluster <- c("Lung_3m_male_14", "Lung_3m_male_12", "Lung_3m_male_17") #just lung
high_var_cluster <- c("Lung_3m_male_14", "Lung_3m_male_12", "Lung_3m_male_17", "Kidney_3m_female_7", "Kidney_3m_male_11", "Marrow_3m_female_15","Marrow_3m_male_8")

# all genes select these special clusters
selected_LVG_1 <- df_main |>
  filter(cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 

#write.csv(selected_LVG_1, 'data/strat_LVG1.csv') # just high var
#write.csv(df_main, 'data/strat_LVG1_high_and_low.csv') # all in LVG 1


# lps genes select special clusters
selected_LVG_1_lps <- df_lps |>
  filter(cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 

print(selected_LVG_1_lps, n = 100)


# is there a big overlap of highly expressed genes in LVG 1 cluster and LPS response genes? 
intersect <- intersect(unique(selected_LVG_1_lps$gene), unique(selected_LVG_1$gene))

overlap <- length(intersect)
all <- length(unique(selected_LVG_1$gene))
perc <- overlap/all
all
perc
# 5% Lung
# 5.2% all




# is it more than exprected statistically
# hypergeometric test

hyper_df <- df_main |>
  filter(cluster_id %in% high_var_cluster)

all_clusters <- unique(hyper_df$cluster_id)
query_genes <- intersect(unique(lps_df$gene), unique(df_main$gene))


hypergeometric_test <- function(cluster_id, df, lps_response_genes) {
  # Filter data for the specific cluster
  df_cluster <- df[df$cluster_id == cluster_id, ]

  # get list of all genes
  all_genes_list <- df %>%
      distinct(gene) %>%
      pull(gene) 

  # get mouse innate genes in dataset
  mouse_innate_genes <- intersect(lps_response_genes, all_genes_list)

  # list of lvg_1s
  lvg_list <- df_cluster  %>% 
      filter(category == "Other_6") %>% 
      distinct(gene) %>% 
      pull(gene)



  # list of overlap between lvg and innate genes
  overlap_list <- intersect(mouse_innate_genes, lvg_list)


  # signature (overall pool) - category of interest = failures
  n <- length(all_genes_list) - length(lvg_list)
  # category of interest
  m <- length(lvg_list) 
  # number of draws
  k <- length(mouse_innate_genes) 
  # overlap 
  q <- 0:(length(overlap_list) - 1) 
  q2 <- length(overlap_list) 

  p_value <- phyper(q2 - 1, m, n, k, lower.tail = FALSE)
  
  return(p_value)
}


p_values_result <- sapply(all_clusters, function(cluster_id) {
  hypergeometric_test(cluster_id, hyper_df, query_genes)})
p_values_result

p_values_df <- data.frame(
  cluster_id = all_clusters,
  p_value = p_values_result,
  significant = p_values_result < 0.05
)
p_values_df


# significant (lung)
# siginificant (all)



# negative control
# is there a big overlap of highly expressed genes in LVG 1 cluster and house keeping genes? 

hk_df <- read_lines("data/scHK_human.symbols.txt")


df_hk <- df_main %>% 
 mutate(hk = str_detect(gene, paste(hk_list, collapse = "|"))) %>%
 filter(hk == TRUE)
head(df_hk)

# all genes select these special clusters
selected_LVG_1 <- df_main |>
  filter(cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 
unique(selected_LVG_1$gene)
#write.csv(selected_LVG_1, 'data/strat_LVG1.csv')

# lps genes select special clusters
selected_LVG_1_hk <- df_hk |>
  filter(cluster_id %in% high_var_cluster) |>
  filter(category == "LVG_1") 

print(selected_LVG_1_hk, n = 100)

intersect <- intersect(unique(selected_LVG_1_hk$gene), unique(selected_LVG_1$gene))

overlap <- length(intersect)
all <- length(unique(selected_LVG_1_hk$gene))
perc <- overlap/all
all
perc
# 100%


