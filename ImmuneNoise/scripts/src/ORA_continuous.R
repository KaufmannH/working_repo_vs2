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




overrepresentaiton_analysis <- function(df, selected_gene_set, selected_category) {

     # total number of genes
     N <- length(unique(df$gene)) 
     # total number of genes in the gene set
     M <- length(df |> filter(`Housekeeping gene` == TRUE) |> pull(gene) |> unique())  
     # number of genes in the bin
     #n <- length(df |> filter(category == selected_category) |> pull(gene) |> unique())  
     n <- length(df |> filter(resvar_bin == selected_category) |> pull(gene) |> unique())  

     # overlap (genes in set & in bin)
     i <- length(intersect(df |> filter(`Housekeeping gene` == TRUE) |> pull(gene) |> unique(), 
                           df |> filter(resvar_bin == selected_category) |> pull(gene) |> unique()))


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




# ENRICHMENT BINS

# Spleen

df <- readRDS('ImmuneNoise/droplet/data/gene_set_df.rds')
data_source <- "droplet"
spleen_df <- df |> filter(tissue == "Spleen")

bin_width <- 0.1 

df_bin <- spleen_df %>%
  select(gene, res_var, `Housekeeping gene`) %>%
  distinct(gene, `Housekeeping gene`, res_var) %>%
  mutate(
    bin_start = floor(res_var / bin_width) * bin_width,
    resvar_bin = sprintf("[%.3f, %.3f)", bin_start, bin_start + bin_width)) %>%
  select(-bin_start)


bins <- sort(unique(df_bin$resvar_bin))

results_df <- map_dfr(bins, function(b) {
  out <- overrepresentaiton_analysis(
    df = df_bin,
    selected_gene_set = NULL,   
    selected_category = b )

  tibble(
    resvar_bin = b,
    p_value = out$p_value,
    expected = out$expected,
    observed = out$observed,
    fold_enrichment = out$fold_enrichment)}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(fold_enrichment))



# get numeric x-position per bin (use bin midpoint)
plot_df <- results_df %>%
  mutate(
    bin_lo = as.numeric(str_match(resvar_bin, "\\[([0-9\\.\\-eE]+),")[,2]),
    bin_hi = as.numeric(str_match(resvar_bin, ",\\s*([0-9\\.\\-eE]+)\\)")[,2]),
    res_var_mid = (bin_lo + bin_hi) / 2 )


plot <- ggplot(plot_df, aes(x = res_var_mid, y = fold_enrichment)) +
  geom_line(colour = '#74C4E7') +
  geom_point(size = 1, color = '#74C4E7') +
   geom_vline(xintercept = c(1, 5), colour = "#1b9e77", linewidth = 0.5, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "#1b9e77", linewidth = 0.5, linetype = "dashed") +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +  
  labs(
    title =   "Over-representation of housekeeping genes",
    x = "res_var (bin midpoint)",
    y = "Fold enrichment") +
  theme_classic(base_size = 20)

ggsave(paste0('ImmuneNoise/', data_source, "/plots/all_ages/spleen_enrichment_hk.png"), plot, width = 12, height = 5)



# Liver

df <- readRDS('ImmuneNoise/droplet/data/gene_set_df.rds')
data_source <- "droplet"
liver_df <- df |> filter(tissue == "Liver")


bin_width <- 0.1 

df_bin <- liver_df %>%
  select(gene, res_var, `Housekeeping gene`) %>%
  distinct(gene, `Housekeeping gene`, res_var) %>%
  mutate(
    bin_start = floor(res_var / bin_width) * bin_width,
    resvar_bin = sprintf("[%.3f, %.3f)", bin_start, bin_start + bin_width)) %>%
  select(-bin_start)


bins <- sort(unique(df_bin$resvar_bin))

results_df <- map_dfr(bins, function(b) {
  out <- overrepresentaiton_analysis(
    df = df_bin,
    selected_gene_set = NULL,   
    selected_category = b )

  tibble(
    resvar_bin = b,
    p_value = out$p_value,
    expected = out$expected,
    observed = out$observed,
    fold_enrichment = out$fold_enrichment)}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(fold_enrichment))



# get numeric x-position per bin (use bin midpoint)
plot_df <- results_df %>%
  mutate(
    bin_lo = as.numeric(str_match(resvar_bin, "\\[([0-9\\.\\-eE]+),")[,2]),
    bin_hi = as.numeric(str_match(resvar_bin, ",\\s*([0-9\\.\\-eE]+)\\)")[,2]),
    res_var_mid = (bin_lo + bin_hi) / 2 )


plot <- ggplot(plot_df, aes(x = res_var_mid, y = fold_enrichment)) +
  geom_line(colour = '#74C4E7') +
  geom_point(size = 1, color = '#74C4E7') +
   geom_vline(xintercept = c(1, 5), colour = "#1b9e77", linewidth = 0.5, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "#1b9e77", linewidth = 0.5, linetype = "dashed") +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +  
  labs(
    title =   "Over-representation of housekeeping genes",
    x = "res_var (bin midpoint)",
    y = "Fold enrichment") +
  theme_classic(base_size = 20)

ggsave(paste0('ImmuneNoise/', data_source, "/plots/all_ages/liver_enrichment_hk.png"), plot, width = 12, height = 5)


