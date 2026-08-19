



library(Seurat)
library(tidyverse)
library(arrow)
library(patchwork)



# 1. Mean and resvar distribution
# ---------------------------------------------------


# Paths
# ---------------------------------------------------

processed_dir <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/bundles'
seurat_dir    <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/seurat_objects'
bootstrap_dir <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/bootstrap'
plot_dir      <- 'droplet/plots/regnoise_vs2'


# Function to do it for all liver units
# ---------------------------------------------------

process_unit <- function(results_path) {

  results_obj <- readRDS(results_path)

  # unit info
  unit <- results_obj |>
    keep(\(x) is.atomic(x) && length(x) == 1) |>
    as_tibble() |>
    rename(n_cells_unit = n_cells)

  id <- unit$unit_id  

  so_path <- list.files(seurat_dir,    pattern = paste0(id, "\\.rds$"),     full.names = TRUE)
  bs_path <- list.files(bootstrap_dir, pattern = paste0(id, "\\.parquet$"), full.names = TRUE)
  so    <- readRDS(so_path)
  bs_df <- read_parquet(bs_path) |> as_tibble()


  # cluster-based info
  clusters <- map_dfr(results_obj$annotation, as_tibble) |>
    left_join(enframe(unlist(results_obj$hvg$per_cluster), "cluster", "n_hvg"), by = "cluster") |>
    rename(n_cells_cluster = n_cells) |>
    rename(cluster_id = cluster) |>
    mutate(cluster = paste(gsub(" ", "_", suggested_label), cluster_id, sep = "_"))

  # bootstrapping
  per_gene <- bs_df |>
    rename(cluster_id = cluster) |>
    left_join(clusters, by = "cluster_id") |>
    bind_cols(unit)

  # gmean
  expression_df <- AverageExpression(so, assays = "SCT", slot = "data",
                                     group.by = "seurat_clusters")$SCT |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "cluster_id", values_to = "mean_expr") |>
    dplyr::mutate(cluster_id = gsub("^g", "", cluster_id))

  assembled_df <- per_gene |>
    left_join(expression_df, by = c("gene", "cluster_id")) |>
    filter(mean_expr > 0)

  # plots
  cluster_cols <- setNames(
    scales::hue_pal()(n_distinct(assembled_df$cluster)),
    sort(unique(assembled_df$cluster)))

  

  plot_resvar <- ggplot(assembled_df, aes(x = res_var, fill = cluster)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 500, colour = NA) +
    scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
    scale_fill_manual(values = cluster_cols) +
    theme_classic(base_size = 20) +
    labs(x = "Residual variance (log10)", y = "Count", fill = "Cluster",
         title = paste(unit$unit_id))

  plot_mean <- ggplot(assembled_df, aes(x = mean_expr, fill = cluster)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 500, colour = NA) +
    scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
    scale_fill_manual(values = cluster_cols) +
    theme_classic(base_size = 20) +
    labs(x = "Mean expression (log10)", y = "Count", fill = "Cluster",
         title = paste(unit$unit_id))

  plot <- plot_mean | plot_resvar
  ggsave(file.path(plot_dir, paste0('resvar_mean_histogram_', id, '.png')),
         plot, width = 20, height = 5)

  invisible(id)
}


# Call function
# ---------------------------------------------------
results_files <- list.files(processed_dir, pattern = '.rds', full.names = TRUE)

walk(results_files, possibly(process_unit, otherwise = NULL))








# 2. Define HVG and LVG
# ---------------------------------------------------


# Paths
# ---------------------------------------------------
processed_dir <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/bundles'
seurat_dir    <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/seurat_objects'
bootstrap_dir <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_20_liver_droplet/bootstrap'
plot_dir      <- 'droplet/plots/regnoise_vs2'

hvg_cutoff <- 0.95
lvg_cutoff <- 0.05
var_cols   <- c(LVG = "#2c7bb6", Intermediate = "grey70", HVG = "#8b0000")




#  Function to do it for all liver units
# ---------------------------------------------------
process_variability <- function(results_path) {

  results_obj <- readRDS(results_path)

  # unit info
  unit <- results_obj |>
    keep(\(x) is.atomic(x) && length(x) == 1) |>
    as_tibble() |>
    rename(n_cells_unit = n_cells)

  id <- unit$unit_id

  so_path <- list.files(seurat_dir,    pattern = paste0(id, "\\.rds$"),     full.names = TRUE)
  bs_path <- list.files(bootstrap_dir, pattern = paste0(id, "\\.parquet$"), full.names = TRUE)

  if (length(so_path) != 1 || length(bs_path) != 1) {
    warning("Skipping ", id, ": expected 1 seurat + 1 bootstrap file, found ",
            length(so_path), " / ", length(bs_path))
    return(invisible(NULL))
  }

  so    <- readRDS(so_path)
  bs_df <- read_parquet(bs_path) |> as_tibble()

  # cluster-based info
  clusters <- map_dfr(results_obj$annotation, as_tibble) |>
    left_join(enframe(unlist(results_obj$hvg$per_cluster), "cluster", "n_hvg"), by = "cluster") |>
    rename(n_cells_cluster = n_cells) |>
    rename(cluster_id = cluster) |>
    mutate(cluster = paste(gsub(" ", "_", suggested_label), cluster_id, sep = "_"))

  # bootstrapping
  per_gene <- bs_df |>
    rename(cluster_id = cluster) |>
    left_join(clusters, by = "cluster_id") |>
    bind_cols(unit)

  # gmean
  expression_df <- AverageExpression(so, assays = "SCT", slot = "data",
                                     group.by = "seurat_clusters")$SCT |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "cluster_id", values_to = "mean_expr") |>
    dplyr::mutate(cluster_id = gsub("^g", "", cluster_id))

  assembled_df <- per_gene |>
    left_join(expression_df, by = c("gene", "cluster_id")) |>
    filter(mean_expr > 0)

  # ---- 2a. mean / resvar histograms --------------------------------
  cluster_cols <- setNames(
    scales::hue_pal()(n_distinct(assembled_df$cluster)),
    sort(unique(assembled_df$cluster)))

  plot_resvar_hist <- ggplot(assembled_df, aes(x = res_var, fill = cluster)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 500, colour = NA) +
    scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
    scale_fill_manual(values = cluster_cols) +
    theme_classic(base_size = 20) +
    labs(x = "Residual variance (log10)", y = "Count", fill = "Cluster", title = id)

  plot_mean <- ggplot(assembled_df, aes(x = mean_expr, fill = cluster)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 500, colour = NA) +
    scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
    scale_fill_manual(values = cluster_cols) +
    theme_classic(base_size = 20) +
    labs(x = "Mean expression (log10)", y = "Count", fill = "Cluster", title = id)

  ggsave(file.path(plot_dir, paste0('resvar_mean_histogram_', id, '.png')),
         plot_mean | plot_resvar_hist, width = 20, height = 5)

  # ---- 2b. Define HVG and LVG --------------------------------------
  assembled_tagged_df <- assembled_df |>
    mutate(
      hi_cut = quantile(res_var, hvg_cutoff, na.rm = TRUE),
      lo_cut = quantile(res_var, lvg_cutoff, na.rm = TRUE),
      hvg    = q010 >= hi_cut,
      lvg    = q090 <= lo_cut,
      variability_class = factor(
        case_when(hvg ~ "HVG", lvg ~ "LVG", TRUE ~ "Intermediate"),
        levels = c("LVG", "Intermediate", "HVG"))) |>
    ungroup() |>
    relocate(all_of(c('gene', 'mean_expr', 'res_var', 'cluster', "n_cells_cluster",
                      'unit_id', 'variability_class', 'hvg', 'lvg', "ci_hi", "ci_lo")))

  # borderline diagnostic (printed per unit)
  message("\n=== ", id, " ===")
  print(table(assembled_tagged_df$variability_class))
  assembled_tagged_df |>
    mutate(
      call = case_when(
        hvg ~ "HVG (stable)",
        res_var >= hi_cut & !hvg ~ "HVG (not stable)",
        lvg ~ "LVG (stable)",
        res_var <= lo_cut & !lvg ~ "LVG (not stable)",
        TRUE ~ "Intermediate")) |>
    count(call) |>
    print()

  # ---- 2c. HVG / LVG cutoff plot -----------------------------------
  cuts <- assembled_tagged_df |>
    group_by(cluster) |>
    summarise(hi_cut = first(hi_cut), lo_cut = first(lo_cut), .groups = "drop") |>
    mutate(hvg_label = paste0(hvg_cutoff * 100, "%"),
           lvg_label = paste0(lvg_cutoff * 100, "%"))

  plot_cutoffs <- ggplot(assembled_tagged_df, aes(x = res_var, fill = variability_class)) +
    geom_rect(data = cuts, aes(xmin = hi_cut, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "#f4a6a6", alpha = 0.2, inherit.aes = FALSE) +
    geom_rect(data = cuts, aes(xmin = 1e-3, xmax = lo_cut, ymin = -Inf, ymax = Inf),
              fill = "#9ecae1", alpha = 0.2, inherit.aes = FALSE) +
    geom_histogram(alpha = 0.8, position = "identity", bins = 100, colour = NA) +
    geom_vline(data = cuts, aes(xintercept = hi_cut),
               colour = "#8b0000", linewidth = 0.01, linetype = "solid") +
    geom_vline(data = cuts, aes(xintercept = lo_cut),
               colour = "#2c7bb6", linewidth = 0.01, linetype = "solid") +
    geom_text(data = cuts, aes(x = hi_cut, y = Inf, label = hvg_label),
              colour = "#8b0000", angle = 0, vjust = 1.2, hjust = -0.5,
              size = 3, inherit.aes = FALSE) +
    geom_text(data = cuts, aes(x = lo_cut, y = Inf, label = lvg_label),
              colour = "#2c7bb6", angle = 0, vjust = 1.2, hjust = 1.1,
              size = 3, inherit.aes = FALSE) +
    scale_x_log10(breaks = scales::breaks_log(n = 5), labels = scales::label_log()) +
    scale_fill_manual(values = var_cols) +
    facet_wrap(~ cluster, ncol = 4) +
    coord_cartesian(xlim = c(1e-3, 10), ylim = c(0, 1200)) +
    theme_classic(base_size = 20) +
    labs(x = "Residual variance", y = "Count", fill = NULL, title = id)

  ggsave(file.path(plot_dir, paste0('hvg_lvg_cutoffs_', id, '.png')),
         plot_cutoffs, width = 17, height = 10)

  invisible(assembled_tagged_df)
}


# Call function
# ---------------------------------------------------

results_files <- list.files(processed_dir, pattern = '\\.rds$', full.names = TRUE)
tagged_list <- map(results_files, possibly(process_variability, otherwise = NULL))



# Save output
# ---------------------------------------------------

assembled_tagged_all <- tagged_list |>   list_rbind()          
write_tsv(assembled_tagged_all, 'droplet/data/assembled_df_annotated_liver_droplet.tsv')





# 4. Check number of cells and HVGs | LVGs per cluster
# ---------------------------------------------------

assembled_tagged_all <- read_tsv('droplet/data/assembled_df_annotated_liver_droplet.tsv')

young_df <- assembled_tagged_all |>
  filter(unit_id %in% c('Liver_3m_female', 'Liver_3m_male')) 

dim(young_df)




dark_col  <- "#3a6f4c"                      
light_col <- colorspace::lighten(dark_col, 0.4) 

cond_cols <- c(Liver_3m_female = light_col, Liver_3m_male = dark_col)
table(plot_df)

plot_df <- young_df |>
  distinct(cluster, unit_id, n_cells_cluster) |>
  mutate( unit_id = factor(unit_id, levels = c("Liver_3m_female", "Liver_3m_male")))

plot <- ggplot(plot_df, aes(x = unit_id, y = n_cells_cluster, fill = unit_id)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_cells_cluster), vjust = -0.3, size = 4) +
  scale_fill_manual(values = cond_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "number of cells", fill = NULL)

ggsave("droplet/plots/regnoise_vs2/cluster_ncells_liver_3m.png", plot, width = 20, height = 10)





assembled_tagged_all <- read_tsv('droplet/data/assembled_df_annotated_liver_droplet.tsv')

# one row per cluster × unit, with age group + sex pulled out of unit_id
plot_base <- assembled_tagged_all |>
  distinct(cluster, unit_id, n_cells_cluster) |>
  mutate(
    age_group = str_extract(unit_id, "\\d+m"),          # "Liver_30m_female" -> "30m"
    sex       = str_extract(unit_id, "female|male"),
    sex       = factor(sex, levels = c("female", "male")))

# colours keyed on sex so they generalise across every age group
dark_col  <- "#3a6f4c"
light_col <- colorspace::lighten(dark_col, 0.4)
sex_cols  <- c(female = light_col, male = dark_col)

plot_dir <- "droplet/plots/regnoise_vs2"

plot_age_group <- function(ag) {

  df <- plot_base |>
    filter(age_group == ag) |>
    droplevels()

  # facets = cell types present in THIS age group only (absent ones just don't appear)
  p <- ggplot(df, aes(x = sex, y = n_cells_cluster, fill = sex)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_cells_cluster), vjust = -0.3, size = 4) +
    scale_fill_manual(values = sex_cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "number of cells", fill = NULL,
         title = paste("Liver", ag))

  ggsave(file.path(plot_dir, paste0("cluster_ncells_liver_", ag, ".png")),
         p, width = 20, height = 10)

  invisible(ag)
}

# 4. Loop over every age group
# ---------------------------------------------------
age_groups <- plot_base |> distinct(age_group) |> pull(age_group)

walk(age_groups, plot_age_group)






