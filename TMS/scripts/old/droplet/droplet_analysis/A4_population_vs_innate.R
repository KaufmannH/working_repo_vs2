# What is the expression level of immune response genes in comparison to lineage specific genes?
# only in Lung, 3m, monocytes, macrophages, DCs


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
mouse_innate_genes <- mouse_innate_genes_df %>% 
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
    select(-organ) %>%
    unique() %>%
    mutate(is_marker = TRUE)
head(markers_df)

df_main_all <- read.csv("data/combined_data.csv")
colnames(df_main_all)

# 
housekeeping_df <- read.csv('data/Housekeeping_TranscriptsMouse.csv', sep = ";") 
housekeeping_list <- housekeeping_df %>%
    pull(Genes)
housekeeping_list

housekeeping_stable_df <- read.csv('data/MostStable_Mouse.csv', sep = ";") 
housekeeping_stable_list <- housekeeping_stable_df %>%
    pull(Gene.name)
housekeeping_stable_list

# once they send the data:
#immune_dict <- readRDS("data/Immunedict/ref_data_Macrophage.RDS")




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


# tag genes
df_main <- df_main_filtered %>%
    # tag house keeping genes
    mutate(housekeeping_gene_stable = (str_detect(gene, paste(housekeeping_stable_list, collapse = "|")))) %>% 
    mutate(housekeeping_gene = (str_detect(gene, paste(housekeeping_list, collapse = "|")))) %>% 
    mutate(mouse_id = str_extract(cluster_id, ".*(?=_[^_]+$)")) 
head(df_main)

# join dfs to tag markers
df_main_markers <- df_main %>%
  left_join(markers_df, by = c("gene", "cell_type")) %>%
  mutate(is_marker = ifelse(is.na(is_marker), FALSE, is_marker))
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
      housekeeping_gene_stable ~ "Stable Housekeeping Gene",
      housekeeping_gene ~ "Housekeeping Gene",
      TRUE ~ "Other")) %>%

  select(gmean, gene, cell_type, gene_set, gene_type, res_var, cluster_id)
head(df_main_with_markers, n = 300)


# TODO: dublicates are issue
dim(df_main_with_markers)
n_distinct(df_main_with_markers)
# overlap between lineage specific an response genes: 0


#saveRDS(df_main_with_markers, "data/df_with_markers_extended.rds") # for iros on 30th jan 25


# check expression range
#min_expr <- min(df_reporting$gmean, na.rm = TRUE)
#max_expr <- max(df_reporting$gmean, na.rm = TRUE) 
#min_expr
#max_expr
#breakpoints <- c(min(df_reporting$gmean), quantile(df_reporting$gmean, 0.75), max(df_reporting$gmean))

 
# expression range
quantile_summary <- df_main_with_markers %>%
  group_by(cluster_id, cell_type) %>%
  summarise(
    Q1 = quantile(gmean, 0.25, na.rm = TRUE),
    Median = quantile(gmean, 0.5, na.rm = TRUE),
    Q3 = quantile(gmean, 0.75, na.rm = TRUE)
  )
quantile_summary


df_with_levels <- df_main_with_markers %>%
  left_join(quantile_summary, by = c("cluster_id", 'cell_type')) %>%
  mutate(expression_level = case_when(
    gmean < Q1 ~ "low",
    gmean > Q3 ~ "high",
    TRUE ~ "medium"
  ))
head(df_with_levels)

sum(df_with_levels$expression_level == "high")
sum(df_with_levels$expression_level == "low")
sum(df_with_levels$expression_level == "medium")






# filter out the genes in "other" bc they are so many
summary_table <- df_with_levels %>%
  group_by(cluster_id, cell_type, gene_set, gene_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(gene_type != "Other")
summary_table

# Bar Plot
my_colors <- c("Lineage Specific Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
summary_table$gene_set <- factor(summary_table$gene_set, 
                                          levels = c("HVG", "Other", "LVG"))

barplot_lin_vs_innate <- ggplot(summary_table, aes(x = gene_set, y = count, fill = gene_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~paste(cell_type, cluster_id , sep = " - "), scales = "free", labeller = label_parsed) + 
  labs(x = "Gene Set",
       y = "Count of Genes",
       fill = "Gene type") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text.x = element_text(size = 6, face = "bold"))
 
ggsave('plots/3_m/boxplot_lineage_vs_innate.png', plot = barplot_lin_vs_innate)





# stacked bar plot

summary_table2 <- df_with_levels %>%
  group_by(cell_type, cluster_id, gene_set, gene_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total_count * 100)


my_colors <- c("Lineage Specific Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
summary_table2$gene_set <- factor(summary_table2$gene_set, levels = c("HVG", "Other", "LVG"))

# Create the 100% stacked barplot
barplot_relative <- ggplot(summary_table2, aes(x = gene_set, y = proportion, fill = gene_type)) +
  geom_bar(stat = "identity", position = "fill") +  # Using 'fill' for 100% stacking
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)
  )

ggsave('plots/3_m/relative_abundances.png', plot = barplot_relative)




