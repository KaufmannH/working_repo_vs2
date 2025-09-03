# What is the expression level and variability of innate immune response genes, lineage specific genes and house keeping genes? 
# only in 3m, macrophages, monocytes, DCs, only innate response genes as total


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
    select(-organ)
head(markers_df)


df_main_all <- read.csv("data/combined_data.csv")
colnames(df_main_all)

# 
housekeeping_df <- read.csv('data/Housekeeping_TranscriptsMouse.csv', sep = ";") 
housekeeping_list <- housekeeping_df %>%
    pull(Genes)
housekeeping_list

immune_dict <- readRDS("data/Immunedict/ref_data_Macrophage.RDS")
immune_dict





# filter for age and cell types
df_main_filtered <- df_main_all  %>% 
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(age = as.numeric(str_remove(age, "m"))) %>% 
  filter(age == 3) %>%
  filter(tissue == "Lung") %>%
  filter(manual_final %in% c("leukocyte_dc", "leukocyte_macrophage", "leukocyte_monocyte")) %>%
  mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
  filter(inner_gene == TRUE) %>%
  mutate(cell_type = manual_final) %>% 
  mutate(cell_type = case_when(
      str_detect(cell_type, "leukocyte_dc") ~ "DC",
      str_detect(cell_type, "leukocyte_macrophage") ~ "Macrophage",
      str_detect(cell_type, "leukocyte_monocyte") ~ "Monocyte",
      TRUE ~ as.character(cell_type))) %>%
  select(cluster_id, cell_type, gene, res_var, gmean, hvg, lvg, num_cells_per_cluster)
head(df_main_filtered)

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
      housekeeping_gene ~ "Housekeeping Gene", 
      is_marker ~ "Lineage Specific Gene",
      TRUE ~ "Other")) %>%
  select(gmean, gene, cell_type, gene_set, gene_type, res_var, cluster_id)
head(df_main_with_markers, n = 300)


#saveRDS(df_main_with_markers, "data/df_with_markers_extended.rds") # for iros on 30th jan 25


# check how many are both innate and lineage specific
#innate_lineage_count <- df_main_with_markers %>%
  #filter(gene_set == "Innate & Lineage Specific") %>%
  #summarise(count = n())
#innate_lineage_count

#df_reporting <- df_main_with_markers %>%
  #group_by(gene, gene_set) %>%
  #summarise(mean_expression = mean(gmean, na.rm = TRUE), .groups = 'drop')
  #pivot_wider(names_from = gene_set, values_from = mean_expression, values_fill = list(mean_expression = 0)) 
#tail(df_reporting, n = 30)


# check the ones that have same value
# every gene is there once within each cluster



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

sum(df_with_levels$expression_level == "high")
sum(df_with_levels$expression_level == "low")
sum(df_with_levels$expression_level == "medium")

head(df_with_levels)



# get gene set categories
my_colors <- c("LVG" = "#5AB4AC", "Other" = "#556967", "HVG" = "#A4C089")
df_with_levels$gene_set <- factor(df_with_levels$gene_set, 
                                          levels = c("HVG", "Other", "LVG"))

boxplot_expression <- ggplot(df_with_levels, aes(x = gene_set, y = gmean, fill = gene_set)) +
  geom_boxplot(color = "black", size = 0.25) + 
   facet_wrap(~paste(cell_type, cluster_id , sep = " - "), scales = "free", labeller = label_parsed) + 
  labs(y = "Gene Mean Expression (gmean)", fill = "Expression Levels" ) +
  theme_classic() +
  scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),  
        strip.text.x = element_text(size = 7, colour = "black")) 
ggsave('plots/3_m/boxplot_lineage_expression.png', plot = boxplot_expression)



# get expression level categories
my_colors <- c("low" = "#5AB4AC", "medium" = "#556967", "high" = "#A4C089")

df_with_levels$expression_level <- factor(df_with_levels$expression_level, 
                                          levels = c("high", "medium", "low"))
boxplot_levels <- ggplot(df_with_levels, aes(x = expression_level, y = gmean, fill = expression_level)) +
  geom_boxplot(color = "black", size = 0.25) + 
  facet_wrap(~paste(cell_type, cluster_id , sep = " - "), scales = "free", labeller = label_parsed) + 
  labs(y = "Gene Mean Expression (gmean)", fill = "Expression Levels" ) +
  theme_classic() +
  scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),  
        strip.text.x = element_text(size = 7, colour = "black")) 
ggsave('plots/3_m/boxplot_lineage_expression_levels.png', plot = boxplot_levels)



# get category numbers
summary_table1 <- df_with_levels %>%
  group_by(cluster_id, expression_level, gene_set, cell_type) %>%
  summarise(count = n(), .groups = 'drop')
summary_table1

# Bar Plot
my_colors <- c("low" = "#5AB4AC", "medium" = "#556967", "high" = "#A4C089")

barplot_counts <- ggplot(summary_table, aes(x = gene_set, y = count, fill = expression_level)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  #facet_wrap(~ cluster_id, scales = "free_y") + 
  facet_wrap(~paste(cell_type, cluster_id , sep = " - "), scales = "free", labeller = label_parsed) + 
  labs(x = "Gene Set",
       y = "Count of Genes",
       fill = "Expression Level") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text.x = element_text(size = 7, face = "bold"))
 
ggsave('plots/3_m/boxplot_lineage_expression_counts.png', plot = barplot_counts)



# get gene type 
summary_table2 <- df_with_levels %>%
  group_by(cluster_id,cell_type, gene_set, gene_type) %>%
  summarise(count = n(), .groups = 'drop')
summary_table2

# Bar Plot
my_colors <- c("Lineage Specific Gene" = "#A4C089", "Other" = "#556967", "Housekeeping Gene" = "#5AB4AC")
summary_table2$gene_set <- factor(summary_table2$gene_set, 
                                          levels = c("HVG", "Other", "LVG"))

barplot_genetype <- ggplot(summary_table2, aes(x = gene_set, y = count, fill = gene_type)) +
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
 
ggsave('plots/3_m/boxplot_lineage_expression_housekeeping.png', plot = barplot_genetype)






# complex heatmap

heatmap_data <- df_main_with_markers %>%
  select(gene, gmean, cluster_id) %>%
  pivot_wider(names_from = cluster_id, values_from = gmean) 
print(heatmap_data, n = 30)


heatmap_matrix <- heatmap_data %>% 
  select(c(-gene) %>%
  as.matrix()
rownames(heatmap_matrix) <- heatmap_data$gene
heatmap_matrix[is.na(heatmap_matrix)] <- 0
head(heatmap_matrix)

df_annotation <- df_main_with_markers %>%
  select(gene_set, gene)
df_annotation

colors <- colorRampPalette(c("navy", "white", "firebrick"))(100)

heatmap <- Heatmap(
  heatmap_matrix,
  col = colors,
  show_row_names = FALSE,
  show_column_names = TRUE
  bottom_annotation = bottom_annotation)
  #row_annotation = row_annotation,

png('plots/3_m/heatmap_lineage_expression_complex.png', width = 1000, height = 800, res = 300)  
draw(heatmap)
dev.off()



# reporting ggplot
heatmap_plot <- ggplot(df_reporting, aes(x = gene, y = mouse_id, fill = mean_expression)) +
  geom_tile() +
  geom_text(aes(label = gene_set), color = "black", size = 3, vjust = 1) + 
  scale_fill_gradientn(colors = colors, values = scales::rescale(breakpoints), name = "Mean Gene Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank()
  ) 

ggsave('plots/3_m/heatmap_lineage_expression_3m.png', plot = heatmap_plot, width = 40, height = 10, dpi = 300)

# to make it wider, add to ggsave
#width = 8, height = 4, dpi = 150      


