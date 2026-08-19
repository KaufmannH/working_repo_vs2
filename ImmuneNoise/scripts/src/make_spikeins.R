

library(Seurat)
library(tidyverse)
library(Matrix)
library(arrow)

so <- readRDS('/home/hkaufm49/analyses/Regnoise_vs2/data/Spleen-3m-female.rds')
#export
saveRDS(seu_spiked, '/home/hkaufm49/analyses/Regnoise_vs2/data/rawdata/Spleen-3m-female.rds')

head(so@meta.data)
cells_cluster <- so@meta.data |>
    group_by(cell_ontology_class) |>
    count()
cells_cluster


# 1. Settings 
# -----------------------------------------------------------

specs <- list(
  list(name = "spike1", group = "macrophage dendritic cell progenitor",
       n_cells = 10,  level = "cluster", dist = "poisson", truncate = TRUE),
  list(name = "spike2", group = "T cell",
       n_cells = 10,  level = "cluster", dist = "nb",      truncate = TRUE),
  list(name = "spike3", group = "T cell",
       n_cells = 500, level = "cluster", dist = "poisson", truncate = TRUE)
)

seu_spiked <- add_spike_genes(so, specs)     # assay inferred from DefaultAssay(seu)

# checks
Assays(seu_spiked)                                   # "originalexp"
DefaultAssay(seu_spiked)                             # "originalexp"
ncol(seu_spiked[["originalexp"]]@meta.features)      # should match the original (e.g. 55)
identical(seu_spiked$orig.ident, so$orig.ident)     # TRUE
Misc(seu_spiked, "spike_info")

si <- as.data.frame(Misc(seu_spiked, "spike_info"))
si
sp <- si$gene       



# 2. Functions 
# -----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- zero-truncated draws (every value >= 1) ----------------------------
rztpois <- function(n, lambda) {
  u <- runif(n, dpois(0, lambda), 1)
  qpois(u, lambda)
}

rztnbinom <- function(n, mu, size) {
  u <- runif(n, dnbinom(0, mu = mu, size = size), 1)
  qnbinom(u, mu = mu, size = size)
}

# ---- reference expression level for a cluster ---------------------------
cluster_ref_level <- function(cm, cells, stat = c("median", "mean"), min_cells = 5) {
  stat <- match.arg(stat)
  sub  <- cm[, cells, drop = FALSE]
  det  <- Matrix::rowSums(sub > 0) >= min_cells
  gm   <- Matrix::rowMeans(sub[det, , drop = FALSE])
  if (stat == "mean") mean(gm) else median(gm)
}

# ---- main ---------------------------------------------------------------
# specs: list of lists, one entry per spike gene. Fields:
#   name      : "spike1"
#   group     : one or more values of group_col to target
#   n_cells   : number of positive cells inside that group
#   level     : "cluster" (match the cluster's typical gene) or a numeric lambda
#   scaling   : "per_positive_cell" (default) | "per_cluster"
#   dist      : "poisson" (default) | "nb" | "constant"
#   truncate  : TRUE (default) = guarantee exactly n_cells detected cells
#   bg_level  : optional lambda for cells outside the target group (default 0)
#
# assay = NULL -> uses (and preserves) the object's current DefaultAssay,
#                 e.g. "originalexp" for anndata / SingleCellExperiment conversions.
#
# The assay is swapped IN PLACE: meta.data, orig.ident, reductions, idents and
# per-gene meta.features are all preserved. Only the counts matrix grows by the
# spike rows.
add_spike_genes <- function(obj, specs,
                            group_col   = "cell_ontology_class",
                            assay       = NULL,
                            ref_stat    = "median",
                            depth_scale = TRUE,
                            nb_size     = 1,
                            keep_idents = TRUE,
                            seed        = 42) {

  set.seed(seed)

  assay <- assay %||% DefaultAssay(obj)          # infer, don't assume "RNA"
  if (!assay %in% Assays(obj))
    stop("assay '", assay, "' not found. Available: ",
         paste(Assays(obj), collapse = ", "))

  cm  <- GetAssayData(obj, assay = assay, layer = "counts")  # slot= on Seurat v4
  md  <- obj@meta.data
  lib <- Matrix::colSums(cm)

  nms <- vapply(specs, `[[`, "", "name")
  if (any(nms %in% rownames(cm)))
    stop("spike name already present in object: ",
         paste(intersect(nms, rownames(cm)), collapse = ", "))

  spike_mat <- Matrix(0, nrow = length(specs), ncol = ncol(cm), sparse = TRUE,
                      dimnames = list(nms, colnames(cm)))
  info <- list()

  for (s in specs) {

    grp_cells <- rownames(md)[!is.na(md[[group_col]]) & md[[group_col]] %in% s$group]
    if (length(grp_cells) < s$n_cells)
      stop(s$name, ": only ", length(grp_cells), " cells in '",
           paste(s$group, collapse = "/"), "' but n_cells = ", s$n_cells)

    pos <- sample(grp_cells, s$n_cells)

    lam <- if (is.numeric(s$level)) s$level
           else cluster_ref_level(cm, grp_cells, stat = ref_stat)

    if ((s$scaling %||% "per_positive_cell") == "per_cluster")
      lam <- lam * length(grp_cells) / s$n_cells

    lam_i <- rep(lam, length(pos))
    if (depth_scale) lam_i <- lam_i * (lib[pos] / mean(lib[grp_cells]))

    dist  <- s$dist     %||% "poisson"
    trunc <- s$truncate %||% TRUE

    counts <- switch(dist,
      poisson  = if (trunc) rztpois(length(pos), lam_i)
                 else       rpois(length(pos), lam_i),
      nb       = if (trunc) rztnbinom(length(pos), lam_i, nb_size)
                 else       rnbinom(length(pos), mu = lam_i, size = nb_size),
      constant = if (trunc) pmax(1, round(lam_i)) else round(lam_i),
      stop("unknown dist: ", dist))

    spike_mat[s$name, pos] <- counts

    if (!is.null(s$bg_level) && s$bg_level > 0) {
      other  <- setdiff(colnames(cm), grp_cells)
      bg_lam <- s$bg_level * if (depth_scale) lib[other] / mean(lib) else 1
      spike_mat[s$name, other] <- rpois(length(other), bg_lam)
    }

    info[[s$name]] <- data.frame(
      gene            = s$name,
      group           = paste(s$group, collapse = "/"),
      n_cluster       = length(grp_cells),
      n_pos_requested = s$n_cells,
      n_pos_realised  = sum(counts > 0),
      lambda          = lam,
      dist            = dist,
      truncated       = trunc,
      mean_in_pos     = if (any(counts > 0)) mean(counts[counts > 0]) else 0,
      mean_in_cluster = sum(counts) / length(grp_cells),
      stringsAsFactors = FALSE)
  }

  # ---- swap the assay in place, keeping everything else intact ----------
  new_counts <- rbind(cm, spike_mat)
  old_assay  <- obj[[assay]]
  v5 <- inherits(old_assay, "Assay5")

  new_obj <- obj                                       # keep meta.data, reductions, idents
  new_obj[[assay]] <- if (v5) CreateAssay5Object(counts = new_counts)
                      else     CreateAssayObject(counts = new_counts)
  DefaultAssay(new_obj) <- assay

  # carry over per-gene meta.features (a fresh assay resets it to 0 cols).
  # spike rows aren't in old_mf -> looked up as NA, the correct default.
  old_mf <- old_assay@meta.features
  if (ncol(old_mf) > 0) {
    new_mf <- new_obj[[assay]]@meta.features                     # full gene set, 0 cols
    new_mf[colnames(old_mf)] <- old_mf[rownames(new_mf), , drop = FALSE]
    new_obj[[assay]]@meta.features <- new_mf
  }

  if (keep_idents) Idents(new_obj) <- Idents(obj)[colnames(new_obj)]

  Misc(new_obj, "spike_info") <- do.call(rbind, info)   # ground truth
  new_obj
}




# check results of spike-ins



# ---- read (unchanged) -------------------------------------------------------
results_obj <- readRDS('/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_14_part3/bundles/Spleen_3m_female.rds')
so <- readRDS("/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_14_part3/seurat_objects/Spleen_3m_female.rds")
bs_path <- '/home/hkaufm49/analyses/Regnoise_vs2/data/processed/results_2026_07_14_part3/bootstrap/Spleen_3m_female.parquet'
bs_df <- read_parquet(bs_path) |> as_tibble()

clusters <- map_dfr(results_obj$annotation, as_tibble) |>
  left_join(enframe(unlist(results_obj$hvg$per_cluster), "cluster", "n_hvg"), by = "cluster")
clusters

run <- results_obj |>
  keep(\(x) is.atomic(x) && length(x) == 1) |>
  as_tibble()

per_gene <- bs_df |>
  mutate(cluster = as.character(cluster)) |>
  left_join(clusters, by = "cluster") |>
  bind_cols(run) |>
  relocate(starts_with("q"), .after = last_col())

expression_df <- AverageExpression(so, assays = "SCT", slot = "data",
                                   group.by = "seurat_clusters")$SCT |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "mean_expr") |>
  dplyr::mutate(cluster = gsub("^g", "", cluster))

assembled_df <- per_gene |>
  left_join(expression_df, by = c("gene", "cluster")) |>
  mutate(
    cluster_label = paste(cluster, gsub("\\s+", "_", trimws(suggested_label)), sep = "_"),
    cluster_label = forcats::fct_reorder(cluster_label, as.integer(cluster)))


# ---- flag spikes ------------------------------------------------------------
plot_df <- assembled_df |>
  filter(mean_expr > 0) |>
  mutate(is_spike = grepl("^spike[0-9]+$", gene))

spikes <- plot_df |>
  filter(is_spike) |>
  mutate(spike_cluster = paste(gene, cluster, sep = "_"))   # one label per line

# one colour per spike x cluster line; everything else is grey
line_ids   <- sort(unique(spikes$spike_cluster))
spike_cols <- setNames(scales::hue_pal()(length(line_ids)), line_ids)

q97 <- quantile(plot_df$res_var, 0.97, na.rm = TRUE)


# ---- resvar -----------------------------------------------------------------
plot_resvar <- ggplot() +
  # background: every gene, grey
  geom_histogram(data = plot_df, aes(x = res_var),
                 fill = "grey80", colour = NA, bins = 500) +
  geom_vline(xintercept = q97, colour = "black", linewidth = 0.3, linetype = "dashed") +
  # spikes: coloured vertical lines on top
  geom_vline(data = spikes, aes(xintercept = res_var, colour = spike_cluster),
             linewidth = 0.8) +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
  scale_colour_manual(values = spike_cols) +
  theme_classic(base_size = 20) +
  labs(x = "Residual variance (log10)", y = "Count",
       colour = "Spike", title = results_obj$unit_id)

# ---- mean expression --------------------------------------------------------
plot_mean <- ggplot() +
  geom_histogram(data = plot_df, aes(x = mean_expr),
                 fill = "grey80", colour = NA, bins = 500) +
  geom_vline(data = spikes, aes(xintercept = mean_expr, colour = spike_cluster),
             linewidth = 0.8) +
  scale_x_log10(breaks = scales::breaks_log(n = 4), labels = scales::label_log()) +
  scale_colour_manual(values = spike_cols) +
  theme_classic(base_size = 20) +
  labs(x = "Mean expression (log10)", y = "Count",
       colour = "Spike", title = results_obj$unit_id)



plot <- plot_mean | plot_resvar
ggsave('droplet/plots/regnoise_vs2/spleen_resvar_dist_spikes.png', plot, width = 20, height = 5)



# check 
# in cluster 3 and 4, how many cells express the spike gene in what clusters


raw <- GetAssayData(so, assay = "originalexp", layer = "counts")["spike3", ]
cl  <- so$seurat_clusters

testing <- tibble(cell = names(raw), count = as.numeric(raw), cluster = cl) |>
  group_by(cluster) |>
  summarise(
    n_cells      = n(),
    n_pos        = sum(count > 0),
    frac_pos     = n_pos / n_cells,
    mean_all     = mean(count),
    mean_in_pos  = if (n_pos > 0) mean(count[count > 0]) else 0,
    max_count    = max(count),
    total_counts = sum(count)) |>
  arrange(desc(frac_pos))


testing