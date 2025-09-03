
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(purrr)


# load tagged df
df_main <- readRDS("data/df_all_cells_lps_tag.rds")
colnames(df_main)
# TEST 6 (in A5 exploring distributions)
# iros code and reporting



# aggregated expression
df_main_filtered <- df_main %>%
    mutate(tissue = str_extract(cluster_id, "^[^_]+")) %>%
    filter(tissue %in% c( "Marrow"))
(df_main_filtered$cluster_id)

long_df_specific <- df_main_filtered %>%
   rename(
    `LPS response early` = lps_stim_early,
    `LPS response late` = lps_stim_late )%>%
  pivot_longer(
    cols = c( "LPS response early", "LPS response late"),
    names_to = "category",
    values_to = "values"
  ) 
head(long_df_specific)


colors <- c('#111D4E', '#4472C4', '#5FB4E5', '#3CBFAE', '#70AD47', '#EE8339', '#8E0C72', "#3B052C")

cleaned_long_df <- long_df_specific %>%
   filter(values == TRUE) 

df_agg <- cleaned_long_df %>%
  group_by(category, manual_final, gene) %>%
  summarise(gmean = mean(gmean, na.rm = TRUE), .groups = "drop") %>%
  mutate(cell_type = recode(manual_final,
    "erythroid_progenitor" = "Erythroid Progenitor",
    "granulocyte_progenitor" = "Granulocyte Progenitor",
    "leukocyte_granulocyte" = "Granulocyte",
    "leukocyte_macrophage" = "Macrophage",
    "leukocyte_monocyte" = "Monocyte",
    "leukocyte_progenitor" = "Leukocyte Progenitor",
    "lymphocyte_B" = "B Lymphocyte",
    "lymphocyte_T" = "T Lymphocyte",
    "mixed" = "Mixed Cells")) %>%
  filter(!cell_type %in% c("Mixed Cells"))
df_agg

df_agg$cell_type <- factor(df_agg$cell_type, levels = c("Erythroid Progenitor","Granulocyte Progenitor", "Granulocyte","Macrophage", "Monocyte",  "Leukocyte Progenitor","B Lymphocyte" , "T Lymphocyte"))
df_agg

output_dir <- "plots/3_m/Test_6/all_clusters"
for (cat in unique(df_agg$category)) {
  df_subset <- df_agg %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cell_type), y = gmean, fill = as.factor(cell_type))) +
    geom_violin(trim = FALSE, size = 0.3) +
    #geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
    #scale_color_manual(values = colors) +
    labs( x = "Cell types", y = "Gene expression") +
    theme_classic() +
    coord_flip() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

  filename <- paste0(output_dir, "/violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}

for (cat in unique(df_agg$category)) {
  df_subset <- df_agg %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cell_type), y = gmean, fill = as.factor(cell_type))) +
    geom_violin(trim = FALSE, size = 0.3) +
    #geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
    #scale_color_manual(values = colors) +
    ylim(NA, quantile(df_subset$gmean, 0.95)) +
    labs( x = "Cell types", y = "Gene expression") +
    theme_classic() +
    coord_flip() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

  filename <- paste0(output_dir, "/zoom_violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}


