

library(Seurat)
library(tidyverse)
library(arrow)
library(patchwork)
library(Matrix)



base_dir      <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_08_12_spleen'
processed_dir <- file.path(base_dir, 'bundles')
seurat_dir    <- file.path(base_dir, 'seurat_objects')
bootstrap_dir <- file.path(base_dir, 'bootstrap')
plot_dir      <- 'droplet/plots/regnoise_vs2'   # <- spleen dir has no "droplet"; kept to match your single-spleen script

# Cutoffs
# ---------------------------------------------------
# NOTE: carried over from your single-dataset spleen script (0.90 / 0.10),
#       NOT the 0.95 / 0.05 used in the liver multi-unit script.
hvg_cutoff <- 0.90
lvg_cutoff <- 0.10
var_cols   <- c(LVG = "#2c7bb6", Intermediate = "grey70", HVG = "#8b0000")


# Function to do it for all spleen units
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

  # ---- mean / resvar histograms ------------------------------------
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

  # ---- Define HVG and LVG ------------------------------------------
  # NOTE: hi_cut / lo_cut are computed globally over the whole unit
  #       (no group_by(cluster)) -- same as both of your originals.
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

  # ---- HVG / LVG cutoff plot ---------------------------------------
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
tagged_list   <- map(results_files, possibly(process_variability, otherwise = NULL))


# Save output
# ---------------------------------------------------
assembled_tagged_all <- tagged_list |> list_rbind()
write_tsv(assembled_tagged_all, 'droplet/data/assembled_df_annotated_spleen_full.tsv')


