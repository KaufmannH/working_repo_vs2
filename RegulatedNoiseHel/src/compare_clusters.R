library(tidyverse)
library(ggplot2)
# plot num of cells per seurat clusters per age sex group 


in_dir  <- "/home/hkaufm49/analyses/working_repo/RegulatedNoiseHel/data/processed/reproduced/seurat_objects"
out_dir <- "/home/hkaufm49/analyses/working_repo/RegulatedNoiseHel/data/processed/reproduced/seurat_objects"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Set these to your meta.data column names
cell_ontology_col <- "cell_ontology_class"     # e.g. "cell_ontology", "celltype", "CellOntology", ...
cluster_col       <- "seurat_clusters"   # or set to your own meta col; fallback below uses Idents()

# ---- packages ----
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# ---- helpers ----
read_seurat_safe <- function(path) {
  tryCatch(readRDS(path), error = function(e) NULL)
}

get_meta_vec <- function(so, col, fallback = NULL) {
  md <- so[[]]
  if (!is.null(col) && col %in% colnames(md)) return(md[[col]])
  if (!is.null(fallback)) return(fallback(so))
  rep(NA_character_, ncol(so))
}

# ---- main ----
files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)

res <- map(files, function(f) {
  so <- read_seurat_safe(f)
  if (is.null(so)) return(NULL)

  obj_name <- tools::file_path_sans_ext(basename(f))

  ontology <- get_meta_vec(so, cell_ontology_col)
  cluster  <- get_meta_vec(
    so,
    cluster_col,
    fallback = function(x) as.character(Idents(x))
  )

  tibble(
    object   = obj_name,
    ontology = as.character(ontology),
    cluster  = as.character(cluster)
  )
}) |> compact()

if (length(res) == 0) stop("No readable .rds Seurat objects found.")

df <- bind_rows(res)

# 1) total cells per object
tbl_total <- df |>
  count(object, name = "cells_total") |>
  arrange(desc(cells_total), object)

# 2) per cell ontology call per object
tbl_ontology <- df |>
  mutate(ontology = if_else(is.na(ontology) | ontology == "", "NA", ontology)) |>
  count(object, ontology, name = "cells") |>
  arrange(object, desc(cells), ontology)

# 3) per cluster per object
tbl_cluster <- df |>
  mutate(cluster = if_else(is.na(cluster) | cluster == "", "NA", cluster)) |>
  count(object, cluster, name = "cells") |>
  arrange(object, desc(cells), cluster)

# 4) ontology x cluster per object
tbl_onto_x_cluster <- df |>
  mutate(
    ontology = if_else(is.na(ontology) | ontology == "", "NA", ontology),
    cluster  = if_else(is.na(cluster)  | cluster  == "", "NA", cluster)
  ) |>
  count(object, ontology, cluster, name = "cells") |>
  arrange(object, desc(cells), ontology, cluster)

# ---- save ----
write.csv(tbl_total,           file.path(out_dir, "cells_total_per_object.csv"), row.names = FALSE)
write.csv(tbl_ontology,        file.path(out_dir, "cells_per_ontology_per_object.csv"), row.names = FALSE)
write.csv(tbl_cluster,         file.path(out_dir, "cells_per_cluster_per_object.csv"), row.names = FALSE)
write.csv(tbl_onto_x_cluster,  file.path(out_dir, "cells_per_ontology_x_cluster_per_object.csv"), row.names = FALSE)

one_tsv <- bind_rows(
  total             = tbl_total,
  ontology          = tbl_ontology,
  cluster           = tbl_cluster,
  ontology_x_cluster= tbl_onto_x_cluster,
  .id = "table"
)

write_tsv(one_tsv, file.path(out_dir, "cell_count_tables.tsv"))

head(one_tsv)


p <- tbl_cluster %>%
  mutate(object = factor(object)) %>%
  ggplot(aes(x = object, y = cells)) +
  geom_boxplot(outlier_alpha = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 2.2) +
  labs(
    x = "Object (age group)",
    y = "Cells per cluster",
    title = "Cells per cluster by object"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave('/home/hkaufm49/analyses/working_repo/RegulatedNoiseHel/data/processed/reproduced/figures/cells_per_cluster.png', p)