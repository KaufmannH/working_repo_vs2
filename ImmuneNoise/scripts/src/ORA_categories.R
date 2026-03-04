# overrepresentation analysis properly
library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(dplyr)

strat_df <- readRDS("ImmuneNoise/pansci/data/strat_df.rds") 
data_source <- "pansci"


long_strat_df <- strat_df |>
    pivot_longer(
      cols      = ends_with(" gene"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    select(gene, gmean, gene_set, res_var, category) 

   all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

long_strat_df <- long_strat_df |> mutate(category = factor(category, levels = all_cats))

col_dict <- c(
  "Housekeeping gene" = "#1f78b4",  
  "Immune response gene"         = "#f5b003", 
  "Control gene"         = "#bfbfbf" ,
   "Other gene"  = "#408709")

unique(long_strat_df$gene_set)

overrepresentaiton_analysis <- function(df, selected_gene_set, selected_category) {

     # total number of genes
     N <- length(unique(df$gene)) 
     # total number of genes in the gene set
     M <- length(df |> filter(gene_set == selected_gene_set) |> pull(gene) |> unique())  
     # number of genes in the bin
     n <- length(df |> filter(category == selected_category) |> pull(gene) |> unique())  
     # overlap (genes in set & in bin)
     i <- length(intersect(df |> filter(gene_set == selected_gene_set) |> pull(gene) |> unique(), 
                           df |> filter(category == selected_category) |> pull(gene) |> unique()))


     # p-value for over-representation 
     p_over <- phyper(i - 1, M, N - M, n, lower.tail = FALSE)

     # expected overlap under null
     expected <- n * (M / N)

     # fold enrichment
     fold_enrichment <- i / expected

     results_df <- list(p_value = p_over,
          expected = expected,
          observed = i,
          fold_enrichment = fold_enrichment)
     return(results_df)
}



result_df <- data.frame()
category_list <- unique(long_strat_df$category)
gene_set_list <- unique(long_strat_df$gene_set)

for (cat in category_list) {
  for (gene_set_name in gene_set_list) {
    result <- overrepresentaiton_analysis(long_strat_df, gene_set_name, cat)
     result$category <- cat
     result$gene_set <- gene_set_name
     result_df <- dplyr::bind_rows(result_df, result)
     }
}
head(result_df)


plot_df <-  result_df %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))

plot <- ggplot(plot_df, aes(x = category, y = fold_enrichment, fill = gene_set)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = col_dict) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(
    aes(label = sig),
    position = position_dodge(width = 0.8),
    vjust = -0.2, size = 3
  ) +
  labs(x = "Category", y = "Fold enrichment", fill = "Gene set") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
     theme(
          strip.background      = element_blank(),
          axis.text.x           = element_text(size = 20, angle = 45, hjust = 1),
          axis.text.y           = element_text(size = 20),
          axis.title.x      = element_text(size = 20),
          axis.title.y      = element_text(size = 20))


  ggsave(paste0('ImmuneNoise/', data_source, "/plots/3_m/Test_16/enrichment_all_clusters.png"), plot, width = 12, height = 5)



