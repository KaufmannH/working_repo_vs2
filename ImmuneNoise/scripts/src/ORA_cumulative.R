

library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(dplyr)




# CUMULATIVE ENRICHMENT

cum_enrichment <- function(df, dataset_name, thresholds,
                               gene_col = "gene",
                               res_var_col = "res_var",
                               gene_set_col) {

  df0 <- df |>
    transmute(
      gene    = .data[[gene_col]],
      res_var = .data[[res_var_col]],
      in_set  = .data[[gene_set_col]]
    ) |>
    distinct(gene, .keep_all = TRUE) |>
    filter(is.finite(res_var))

  # coerce in_set to logical robustly (TRUE/FALSE, 1/0, "TRUE"/"FALSE")
  df0 <- df0 |>
    mutate(in_set = dplyr::case_when(
      is.logical(in_set) ~ in_set,
      is.numeric(in_set) ~ in_set != 0,
      is.character(in_set) ~ tolower(in_set) %in% c("true", "t", "1", "yes", "y"),
      TRUE ~ as.logical(in_set)
    ))

  N <- n_distinct(df0$gene)
  M <- df0 |>
    filter(in_set == TRUE) |>
    summarise(M = n_distinct(gene)) |>
    pull(M)

  out <- purrr::map_dfr(thresholds, function(t) {

    sub <- df0 |> filter(res_var < t)

    n <- n_distinct(sub$gene)
    k <- sub |>
      filter(in_set == TRUE) |>
      summarise(k = n_distinct(gene)) |>
      pull(k)

    expected <- n * (M / N)
    fold_enrichment <- if (n == 0 || M == 0) NA_real_ else (k / n) / (M / N)
    p_value <- if (n == 0) NA_real_ else stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)

    tibble(
      dataset = dataset_name,
      threshold = t,
      n_left = n,
      observed = k,
      expected = expected,
      fold_enrichment = fold_enrichment,
      p_value = p_value
    )
  }) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  out
}


# 1. liver

# load data 

df_facs_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/combined_data.csv")
df_droplet_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv") 
df_pansci_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/pansci/data/prepped_for_strat/combined_data.csv") 

  # select only hepatocytes (delete in future)
df_facs <- df_facs_raw |> mutate(technology = 'facs') |> filter(cell_type == 'hepatocyte')
df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') |> filter(cell_type == 'hepatocyte')
df_easysci <- df_pansci_raw |> mutate(technology = 'easysci') # already selected hepatocytes in previous steps



# tag gene sets

# general house keeping genes
housekeeping_df <- read.csv('ImmuneNoise/reference_gene_sets/Housekeeping_TranscriptsMouse.csv', sep = ";") 
housekeeping_list <- housekeeping_df |>
    pull(Genes)

df_facs_tagged <- df_facs |>  mutate(Housekeeping_gene = gene %in% housekeeping_list)
df_droplet_tagged <- df_droplet |>  mutate(Housekeeping_gene = gene %in% housekeeping_list)
df_easysci_tagged <- df_easysci |>  mutate(Housekeeping_gene = gene %in% housekeeping_list)


datasets <- list(FACS = df_facs_tagged, 
                Droplet =  df_droplet_tagged, 
                EasySci = df_easysci_tagged)
head(datasets$FACS)



# run
selected_gene_set <- "Housekeeping_gene"
bin_width <- 0.1

all_res <- purrr::imap_dfr(datasets, function(df, nm) {
  df |> transmute(dataset = nm, res_var = .data[["res_var"]])
}) |>
  filter(is.finite(res_var))

thr_min <- floor(min(all_res$res_var) / bin_width) * bin_width + bin_width
thr_max <- floor(max(all_res$res_var) / bin_width) * bin_width + bin_width
thresholds <- seq(thr_min, thr_max, by = bin_width)

cum_df_all <- purrr::imap_dfr(datasets, function(df, nm) {
  cum_enrichment(
    df = df,
    dataset_name = nm,
    thresholds = thresholds,
    gene_set_col = selected_gene_set
  )
})



plot <- ggplot(cum_df_all, aes(x = threshold, y = fold_enrichment, colour = dataset, group = dataset)) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = c(1, 5), colour = "black", linewidth = 0.5, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = c(0.5, 1, 3, 5), y = Inf,   label = c("0.5", '1' ,"3", '5'),  hjust = -0.5, vjust = 1, colour = "grey40", size = 4) +
  scale_x_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_log()) +
  scale_colour_manual(values = c(FACS = "#74C4E7", Droplet = "#F7D639", EasySci = "#94BC47")) +
  theme_classic(base_size = 16) +
  labs(
    title = 'Cumulative enrichment of housekeeping genes in hepatocytes',
    x = "Residual variance",
    y = "Fold enrichment") +
  theme( plot.title = element_text(size = 14))


ggsave(('ImmuneNoise/comparison/enrichment/liver_cum_enrichment_hk.png'), plot, width = 12, height = 5)





# 2. spleen

# load data 




# run
selected_gene_set <- "Housekeeping_gene"
bin_width <- 0.1

all_res <- purrr::imap_dfr(datasets, function(df, nm) {
  df |> transmute(dataset = nm, res_var = .data[["res_var"]])
}) |>
  filter(is.finite(res_var))

thr_min <- floor(min(all_res$res_var) / bin_width) * bin_width + bin_width
thr_max <- floor(max(all_res$res_var) / bin_width) * bin_width + bin_width
thresholds <- seq(thr_min, thr_max, by = bin_width)

cum_df_all <- purrr::imap_dfr(datasets, function(df, nm) {
  cum_enrichment(
    df = df,
    dataset_name = nm,
    thresholds = thresholds,
    gene_set_col = selected_gene_set
  )
})


plot <- ggplot(cum_df_all, aes(x = threshold, y = fold_enrichment, colour = dataset, group = dataset)) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = c(1, 5), colour = "black", linewidth = 0.5, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = c(0.5, 1, 3, 5), y = Inf,   label = c("0.5", '1' ,"3", '5'),  hjust = -0.5, vjust = 1, colour = "grey40", size = 4) +
  scale_x_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_log()) +
  scale_colour_manual(values = c(FACS = "#74C4E7", Droplet = "#F7D639", `10X` = "#E77E01")) +
  theme_classic(base_size = 16) +
  labs(
    title = 'Cumulative enrichment of housekeeping genes in spleen',
    x = "Residual variance",
    y = "Fold enrichment") +
  theme( plot.title = element_text(size = 14))


ggsave(('ImmuneNoise/comparison/enrichment/spleen_cum_enrichment_hk.png'), plot, width = 12, height = 5)







# 0ld

head(liver_df)
bin_width <- 0.1

# 1) One row per gene, keep HK flag + res_var, make numeric thresholds (bin upper edges)
df0 <- liver_df |>
  select(gene, res_var, `Housekeeping gene`) |>
  distinct(gene, .keep_all = TRUE) |>
  filter(is.finite(res_var)) |>
  mutate(
    bin_start = floor(res_var / bin_width) * bin_width,
    thr = bin_start + bin_width   # this is the threshold t (upper edge)
  )

# thresholds to evaluate (sorted unique upper edges)
thresholds <- sort(unique(df0$thr))

head(df0)


# 2) Background totals
N <- n_distinct(df0$gene)                     
M <- df0 |>
  dplyr::filter(`Housekeeping gene` == TRUE) |>
  dplyr::summarise(M = dplyr::n_distinct(gene)) |>
  dplyr::pull(M)

# 3) Cumulative enrichment for "everything left of threshold"
cum_df <- purrr::map_dfr(thresholds, function(t) {

  sub <- df0 |> dplyr::filter(res_var < t)

  n <- dplyr::n_distinct(sub$gene)

  k <- sub |>
    dplyr::filter(`Housekeeping gene` == TRUE) |>
    dplyr::pull(gene) |>
    dplyr::n_distinct()

  expected <- n * (M / N)
  fold_enrichment <- if (n == 0 || M == 0) NA_real_ else (k / n) / (M / N)

  p_value <- if (n == 0) NA_real_ else stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)

  tibble::tibble(
    threshold = t,
    n_left = n,
    observed = k,
    expected = expected,
    fold_enrichment = fold_enrichment,
    p_value = p_value
  )
}) |>
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

# 4) Plot (fold enrichment over threshold)
plot <- ggplot(cum_df, aes(x = threshold, y = fold_enrichment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_classic(base_size = 16) +
  labs(
    x = "Residual variance threshold (t): genes with res_var < t",
    y = "Fold enrichment of housekeeping genes"
  )


ggsave(paste0('ImmuneNoise/', data_source, "/plots/all_ages/liver_cum_enrichment_hk.png"), plot, width = 12, height = 5)












# CUMULATIVE RESVAR

# 1. Liver

# load data 

df_facs_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/combined_data.csv")
df_droplet_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv") 
df_pansci_raw <- read_csv("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/pansci_preprocessing/data/prepped_for_strat/combined_data.csv") 

  # select only hepatocytes (delete in future)
df_facs <- df_facs_raw |> mutate(technology = 'facs') |> filter(cell_type == 'hepatocyte')
df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') |> filter(cell_type == 'hepatocyte')
df_easysci <- df_pansci_raw |> mutate(technology = 'easysci') # already selected hepatocytes in previous steps


datasets <- list(
  FACS    = df_facs,
  Droplet = df_droplet,
  EasySci = df_easysci)


cum_df_all <- imap_dfr(datasets, ~{
  d <- .x |>
    select(gene, cluster_id, res_var) |>
    filter(is.finite(res_var), res_var > 0) |>
    arrange(res_var)

  n_total <- nrow(d)

  d |>
    mutate(
      technology = .y,
      cum_n      = row_number(),
      cum_frac   = if (n_total > 0) cum_n / n_total else NA_real_
    )
})


plot <- ggplot(cum_df_all, aes(x = res_var, y = cum_frac, colour = technology, group = technology)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  geom_line(linewidth = 0.7) +
  coord_cartesian(xlim = c(0.00000001, 100))  +
  geom_vline(xintercept = c(1, 5), colour = "black", linewidth = 0.3, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "black", linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = c(0.5, 1, 3, 5), y = Inf,   label = c("0.5", '1' ,"3", '5'),  hjust = -0.1, vjust = 1, colour = "grey40", size = 3.5) +
  scale_x_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_log()) +
  scale_colour_manual(values = c(FACS = "#74C4E7", Droplet = "#F7D639", EasySci = "#94BC47")) +
  theme_classic(base_size = 20) +
  labs(
    title = "Cumulative residual variance - Hepatocytes",
    x = "Residual variance",
    y = "Cumulative fraction of genes")


ggsave('ImmuneNoise/comparison/enrichment/liver_cum_resvar.png', plot, width = 12, height = 5)







# 1. Spleen


# load data 

df_facs_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/combined_data.csv")
df_droplet_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv") 
df_10X_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/spleen_10X/data/combined_data.csv")

  # select only spleen 
df_facs <- df_facs_raw |> mutate(technology = 'facs') |> filter(tissue == 'Spleen')
df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') |> filter(tissue == 'Spleen')
df_10X <- df_10X_raw |> mutate(technology = '10X') 


datasets <- list(
  FACS    = df_facs,
  Droplet = df_droplet,
  `10X`   = df_10X)


cum_df_all <- imap_dfr(datasets, ~{
  d <- .x |>
    select(gene, cluster_id, res_var) |>
    filter(is.finite(res_var), res_var > 0) |>
    arrange(res_var)

  n_total <- nrow(d)

  d |>
    mutate(
      technology = .y,
      cum_n      = row_number(),
      cum_frac   = if (n_total > 0) cum_n / n_total else NA_real_
    )
})



plot <- ggplot(cum_df_all, aes(x = res_var, y = cum_frac, colour = technology, group = technology)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = c(1, 5), colour = "black", linewidth = 0.3, linetype = "solid") +
  geom_vline(xintercept = c(0.5, 3), colour = "black", linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = c(0.5, 1, 3, 5), y = 0.95, label = c("0.5", '1' ,"3", '5'),  hjust = c(1.2, -0.3, 1.2, -0.3), colour = "grey40", size = 4) +
  scale_x_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_log()) +
  coord_cartesian(xlim = c(0.00000001, 100)) +
  scale_colour_manual(values = c(FACS = "#74C4E7", Droplet = "#F7D639",   `10X`  = "#E77E01")) +
  theme_classic(base_size = 16) +
  labs(
    title = "Cumulative residual variance - Spleen",
    x = "Residual variance",
    y = "Cumulative fraction of genes")


ggsave('ImmuneNoise/comparison/enrichment/spleen_cum_resvar.png', plot, width = 12, height = 5)





# trying the felxible CUTOFF based on the cumulative dist
# LIVER

df_facs_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/combined_data.csv")
df_droplet_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv") 
df_pansci_raw <- read_csv("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/pansci_preprocessing/data/prepped_for_strat/combined_data.csv") 

  # select only hepatocytes (delete in future)
df_facs <- df_facs_raw |> mutate(technology = 'facs') |> filter(cell_type == 'hepatocyte')
df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') |> filter(cell_type == 'hepatocyte')
df_easysci <- df_pansci_raw |> mutate(technology = 'easysci') # already selected hepatocytes in previous steps


# select threshold for the gmean cutoff 
gmean_threshold <- 0.2 

datasets <- list(
  FACS    = df_facs,
  Droplet = df_droplet,
  EasySci = df_easysci) |>
  map(~ filter(.x, is.finite(gmean), gmean >= gmean_threshold))

min(datasets$FACS$gmean)

# what genes still have super low var
low_resvar_genes <- imap_dfr(datasets, ~ {
  .x |>
    filter(res_var < 1e-5) |>
    select(gene, cluster_id, res_var, gmean) |>
    mutate(technology = .y)
})

low_resvar_genes



cum_df_all <- imap_dfr(datasets, ~{
  d <- .x |>
    select(gene, cluster_id, res_var) |>
    filter(is.finite(res_var), res_var > 0) |>
    arrange(res_var)

  n_total <- nrow(d)

  d |>
    mutate(
      technology = .y,
      cum_n      = row_number(),
      cum_frac   = if (n_total > 0) cum_n / n_total else NA_real_
    )
})


find_density_cutoff <- function(x, frac = 0.3, n = 512) {
  d <- density(log10(x[x > 0 & is.finite(x)]), n = n)
  
  peak_idx <- which.max(d$y)
  target <- frac * d$y[peak_idx]
  
  # upper cutoff: search right of peak
  upper_region <- (peak_idx + 1):length(d$y)
  upper_idx <- upper_region[which.min(abs(d$y[upper_region] - target))]
  
  # lower cutoff: search left of peak
  lower_region <- 1:(peak_idx - 1)
  lower_idx <- lower_region[which.min(abs(d$y[lower_region] - target))]
  
  list(
    res_var_upper = 10^d$x[upper_idx],
    res_var_lower = 10^d$x[lower_idx],
    density_at_upper = d$y[upper_idx],
    density_at_lower = d$y[lower_idx],
    peak_density = d$y[peak_idx],
    peak_res_var = 10^d$x[peak_idx]
  )
}


frac <- 0.3

cutoffs <- imap_dfr(datasets, ~ {
  c <- find_density_cutoff(.x$res_var, frac = frac)
  tibble(
    technology = .y,
    res_var_lower = c$res_var_lower,
    res_var_upper = c$res_var_upper,
    peak_res_var = c$peak_res_var
  )
})

cutoffs


# reshape to long for plotting both cutoffs at once
cutoffs_long <- cutoffs |>
  pivot_longer(cols = c(res_var_lower, res_var_upper),
               names_to = "side", values_to = "res_var_cutoff")


density_df <- imap_dfr(datasets, ~ mutate(.x, technology = .y))


tech_cols_fill <- c(FACS = "#74C4E7", Droplet = "#F7D639", EasySci = "#94BC47")

plot <- ggplot(density_df, aes(x = res_var, fill = technology, color = technology)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +  
  geom_vline(data = cutoffs_long,
             aes(xintercept = res_var_cutoff, colour = technology),
             linetype = "dotted", linewidth = 0.8) +
  geom_text(data = cutoffs_long,
            aes(x = res_var_cutoff, y = Inf, colour = technology,
                label = round(res_var_cutoff, 2)),
            vjust = 1.5, hjust = -0.1, size = 3.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = tech_cols_fill) +
  scale_colour_manual(values = tech_cols_fill) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = c(10^-3, NA)) +
  guides(colour = "none") +  
  labs(
    x = "Residual variance",
    y = "Density",
    fill = "Technology")

ggsave('ImmuneNoise/comparison/enrichment/liver_density.png', plot, width = 12, height = 5)





# SPLEEN

df_facs_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/combined_data.csv")
df_droplet_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv") 
df_10X_raw <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/spleen_10X/data/combined_data.csv")

  # select only spleen 
df_facs <- df_facs_raw |> mutate(technology = 'facs') |> filter(tissue == 'Spleen')
df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') |> filter(tissue == 'Spleen')
df_10X <- df_10X_raw |> mutate(technology = '10X') 


# select threshold for the gmean cutoff
gmean_threshold <- 0.2

datasets <- list(
  FACS    = df_facs,
  Droplet = df_droplet,
  `10X`   = df_10X) |>
  map(~ filter(.x, is.finite(gmean), gmean >= gmean_threshold))

(datasets$FACS$cluster_id)

# what genes still have super low var
low_resvar_genes <- imap_dfr(datasets, ~ {
  .x |>
    filter(res_var < 1e-5) |>
    select(gene, cluster_id, res_var, gmean) |>
    mutate(technology = .y)
})

low_resvar_genes


cum_df_all <- imap_dfr(datasets, ~{
  d <- .x |>
    select(gene, cluster_id, res_var) |>
    filter(is.finite(res_var), res_var > 0) |>
    arrange(res_var)

  n_total <- nrow(d)

  d |>
    mutate(
      technology = .y,
      cum_n      = row_number(),
      cum_frac   = if (n_total > 0) cum_n / n_total else NA_real_
    )
})


find_density_cutoff <- function(x, frac = 0.3, n = 512) {
  d <- density(log10(x[x > 0 & is.finite(x)]), n = n)

  peak_idx <- which.max(d$y)
  target <- frac * d$y[peak_idx]

  upper_region <- (peak_idx + 1):length(d$y)
  upper_idx <- upper_region[which.min(abs(d$y[upper_region] - target))]

  lower_region <- 1:(peak_idx - 1)
  lower_idx <- lower_region[which.min(abs(d$y[lower_region] - target))]

  list(
    res_var_upper = 10^d$x[upper_idx],
    res_var_lower = 10^d$x[lower_idx],
    density_at_upper = d$y[upper_idx],
    density_at_lower = d$y[lower_idx],
    peak_density = d$y[peak_idx],
    peak_res_var = 10^d$x[peak_idx]
  )
}


frac <- 0.1  # percentage of the peak that it drops to

cutoffs <- imap_dfr(datasets, ~ {
  c <- find_density_cutoff(.x$res_var, frac = frac)
  tibble(
    technology = .y,
    res_var_lower = c$res_var_lower,
    res_var_upper = c$res_var_upper,
    peak_res_var = c$peak_res_var
  )
})

cutoffs


cutoffs_long <- cutoffs |>
  pivot_longer(cols = c(res_var_lower, res_var_upper),
               names_to = "side", values_to = "res_var_cutoff")


density_df <- imap_dfr(datasets, ~ mutate(.x, technology = .y))


tech_cols_fill <- c(FACS = "#74C4E7", Droplet = "#F7D639", `10X` = "#E77E01")

plot <- ggplot(density_df, aes(x = res_var, fill = technology, color = technology)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
  geom_vline(data = cutoffs_long,
             aes(xintercept = res_var_cutoff, colour = technology),
             linetype = "dotted", linewidth = 0.8) +
  geom_text(data = cutoffs_long,
            aes(x = res_var_cutoff, y = Inf, colour = technology,
                label = round(res_var_cutoff, 2)),
            vjust = 1.5, hjust = -0.1, size = 3.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = tech_cols_fill) +
  scale_colour_manual(values = tech_cols_fill) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = c(10^-3, NA)) +
  guides(colour = "none") +
  labs(
    x = "Residual variance",
    y = "Density",
    fill = "Technology")

ggsave('ImmuneNoise/comparison/enrichment/spleen_density.png', plot, width = 12, height = 5)






# faceted by sex (bc the sct ran per sex group)

# select threshold for the gmean cutoff
gmean_threshold <- 0.2

datasets <- list(
  FACS    = df_facs,
  Droplet = df_droplet,
  `10X`   = df_10X) |>
  map(~ .x |>
        filter(is.finite(gmean), gmean >= gmean_threshold) |>
        mutate(sex = if_else(str_detect(cluster_id, "male") & !str_detect(cluster_id, "female"),
                             "male", "female")))


# what genes still have super low var
low_resvar_genes <- imap_dfr(datasets, ~ {
  .x |>
    filter(res_var < 1e-5) |>
    select(gene, cluster_id, res_var, gmean, sex) |>
    mutate(technology = .y)
})

low_resvar_genes


find_density_cutoff <- function(x, frac = 0.3, n = 512) {
  d <- density(log10(x[x > 0 & is.finite(x)]), n = n)

  peak_idx <- which.max(d$y)
  target <- frac * d$y[peak_idx]

  upper_region <- (peak_idx + 1):length(d$y)
  upper_idx <- upper_region[which.min(abs(d$y[upper_region] - target))]

  lower_region <- 1:(peak_idx - 1)
  lower_idx <- lower_region[which.min(abs(d$y[lower_region] - target))]

  list(
    res_var_upper = 10^d$x[upper_idx],
    res_var_lower = 10^d$x[lower_idx],
    density_at_upper = d$y[upper_idx],
    density_at_lower = d$y[lower_idx],
    peak_density = d$y[peak_idx],
    peak_res_var = 10^d$x[peak_idx]
  )
}


frac <- 0.1


# compute cutoffs per technology AND sex
cutoffs <- imap_dfr(datasets, ~ {
  .x |>
    group_by(sex) |>
    group_modify(~ {
      c <- find_density_cutoff(.x$res_var, frac = frac)
      tibble(
        res_var_lower = c$res_var_lower,
        res_var_upper = c$res_var_upper,
        peak_res_var = c$peak_res_var
      )
    }) |>
    ungroup() |>
    mutate(technology = .y)
})

cutoffs


cutoffs_long <- cutoffs |>
  pivot_longer(cols = c(res_var_lower, res_var_upper),
               names_to = "side", values_to = "res_var_cutoff")


density_df <- imap_dfr(datasets, ~ mutate(.x, technology = .y))


tech_cols_fill <- c(FACS = "#74C4E7", Droplet = "#F7D639", `10X` = "#E77E01")


make_density_plot <- function(density_data, cutoff_data, sex_label) {
  ggplot(density_data, aes(x = res_var, fill = technology, color = technology)) +
    geom_density(alpha = 0.5, linewidth = 0.5) +
    scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
    geom_vline(data = cutoff_data,
               aes(xintercept = res_var_cutoff, colour = technology),
               linetype = "dotted", linewidth = 0.8) +
    geom_text(data = cutoff_data,
              aes(x = res_var_cutoff, y = Inf, colour = technology,
                  label = round(res_var_cutoff, 2)),
              vjust = 1.5, hjust = -0.1, size = 3.5,
              inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = tech_cols_fill) +
    scale_colour_manual(values = tech_cols_fill) +
    theme_classic(base_size = 20) +
    coord_cartesian(xlim = c(10^-3, NA)) +
    guides(colour = "none") +
    labs(
      x = "Residual variance",
      y = "Density",
      fill = "Technology",
      title = sex_label)
}


p_male <- make_density_plot(
  density_df |> filter(sex == "male"),
  cutoffs_long |> filter(sex == "male"),
  "Male")

p_female <- make_density_plot(
  density_df |> filter(sex == "female"),
  cutoffs_long |> filter(sex == "female"),
  "Female")


combined_plot <- p_male / p_female

ggsave('ImmuneNoise/comparison/enrichment/spleen_density_by_sex.png',
       combined_plot, width = 12, height = 10)



# the peaks look asymmetric. diagnose. 

d <- density(log10(datasets$FACS$res_var[datasets$FACS$res_var > 0]))
plot(d)
abline(v = log10(c(cutoffs$res_var_lower[1], cutoffs$res_var_upper[1])), col = "red")
abline(h = 0.3 * max(d$y), col = "blue", lty = 2)