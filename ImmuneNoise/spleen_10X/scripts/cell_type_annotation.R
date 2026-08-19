# annotating the additional spleen 10X dataset: lymphoid vs myeloid lineages

library(dplyr)
library(readr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(RColorBrewer)



# ── 0. Setup ─────────────────────────────────────────────────────────────────
plot_dir <- "ImmuneNoise/spleen_10X/plots/cell_type_annotation" 
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

so <- readRDS("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_additional_vs1/seurat_objects/Spleen_2m_male_processed.rds")

cat("Object loaded:", ncol(so), "cells,", nrow(so), "genes\n")

DefaultAssay(so) <- "SCT"

# ── 1. Lineage marker signatures ────────────────────────────────────────────
# Broad, well-established markers for mouse spleen lineages.

lineage_markers <- list(
  Lymphoid = c(
    "Cd3d", "Cd3e", "Cd3g",           # T cells
    "Cd4", "Cd8a", "Cd8b1",           # T cell subtypes
    "Ncr1", "Klrb1c", "Nkg7",         # NK
    "Cd19", "Ms4a1", "Cd79a", "Cd79b",# B cells
    "Jchain", "Sdc1"                   # Plasma
  ),
  Myeloid = c(
    "Lyz2", "Csf1r", "Ly6c2",         # Monocytes
    "Adgre1", "Cd68", "C1qa", "C1qb", # Macrophages
    "Itgax", "H2-Aa", "H2-Ab1",       # DCs / MHC-II
    "Xcr1", "Siglech", "Bst2",        # DC subtypes
    "S100a8", "S100a9", "Ly6g"         # Neutrophils
  ),
  Other = c(
    "Hba-a1", "Hba-a2", "Hbb-bs",     # Erythroid
    "Pf4", "Ppbp",                     # Platelets
    "Pecam1", "Cdh5",                  # Endothelial
    "Col1a1", "Dcn"                    # Fibroblast
  )
)

# Filter to genes present in the object
lineage_markers_filt <- lapply(lineage_markers, function(g) g[g %in% rownames(so)])
for (nm in names(lineage_markers_filt)) {
  cat(nm, ":", length(lineage_markers_filt[[nm]]), "/",
      length(lineage_markers[[nm]]), "markers found\n")
}

# ── 2. Score and assign ─────────────────────────────────────────────────────
for (nm in names(lineage_markers_filt)) {
  so <- AddModuleScore(
    so,
    features = list(lineage_markers_filt[[nm]]),
    name     = paste0("score_", nm),
    seed     = 42
  )
}

# AddModuleScore appends "1" — grab the score columns
score_df <- data.frame(
  Lymphoid = so@meta.data[["score_Lymphoid1"]],
  Myeloid  = so@meta.data[["score_Myeloid1"]],
  Other    = so@meta.data[["score_Other1"]]
)

# Direct assignment — no mapping step needed
so$lineage <- colnames(score_df)[apply(score_df, 1, which.max)]
so$lineage <- factor(so$lineage, levels = c("Lymphoid", "Myeloid", "Other"))

cat("\n── Assignment summary ──\n")
print(table(so$lineage))

# ── 3. Colors ────────────────────────────────────────────────────────────────
lineage_cols <- c(Lymphoid = "#4DBBD5", Myeloid = "#E64B35", Other = "#91D1C2")

# ── 4. PLOT 1: UMAP ─────────────────────────────────────────────────────────
cat("\nGenerating UMAP...\n")

p_umap <- DimPlot(
  so,
  reduction  = "umap",
  group.by   = "lineage",
  cols       = lineage_cols,
  pt.size    = 0.3,
  label      = TRUE,
  label.size = 5,
  repel      = TRUE
) +
  ggtitle("Mouse Spleen — Broad Cell Type Annotation") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(plot_dir, "umap_celltypes.png"),
       plot = p_umap, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "umap_celltypes.pdf"),
       plot = p_umap, width = 10, height = 8)
cat("  Saved: umap_celltypes\n")

# ── 5. PLOT 2: DotPlot — top markers per lineage ────────────────────────────
cat("Generating DotPlot...\n")

dot_markers <- list(
  Lymphoid = c("Cd3e", "Cd4", "Cd8a", "Ncr1", "Cd19", "Ms4a1", "Cd79a", "Jchain"),
  Myeloid  = c("Lyz2", "Csf1r", "Adgre1", "Cd68", "Itgax", "Siglech", "S100a8", "S100a9"),
  Other    = c("Hba-a1", "Hbb-bs", "Pf4", "Pecam1", "Col1a1")
)
dot_genes <- unique(unlist(dot_markers))
dot_genes <- dot_genes[dot_genes %in% rownames(so)]

Idents(so) <- "lineage"

p_dot <- DotPlot(so, features = dot_genes, cols = c("lightgrey", "#CB181D"),
                 dot.scale = 6) +
  RotatedAxis() +
  ggtitle("Lineage Marker Expression") +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x = element_text(size = 9, face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11)
  )


ggsave(file.path(plot_dir, "dotplot_markers.pdf"),
       plot = p_dot, width = 14, height = 5)
cat("  Saved: dotplot_markers\n")

# ── 6. PLOT 3: FeaturePlot — one key marker per lineage ─────────────────────
cat("Generating FeaturePlots...\n")

feat_genes <- c("Cd3e", "Ms4a1", "Lyz2", "Itgax", "S100a8", "Hba-a1")
feat_genes <- feat_genes[feat_genes %in% rownames(so)]

p_feat <- FeaturePlot(
  so,
  features  = feat_genes,
  reduction = "umap",
  pt.size   = 0.15,
  order     = TRUE,
  ncol      = 3,
  cols      = c("lightgrey", "#08519C")
) &
  theme(
    plot.title = element_text(size = 11, face = "italic"),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

n_rows <- ceiling(length(feat_genes) / 3)
ggsave(file.path(plot_dir, "featureplot_key_markers.png"),
       plot = p_feat, width = 14, height = n_rows * 4.5, dpi = 300)
ggsave(file.path(plot_dir, "featureplot_key_markers.pdf"),
       plot = p_feat, width = 14, height = n_rows * 4.5)
cat("  Saved: featureplot_key_markers\n")

# ── 7. Save annotated object ────────────────────────────────────────────────
out_rds <- "ImmuneNoise/spleen_10X/data/Spleen_2m_male_annotated.rds"
saveRDS(so, out_rds)
cat("\nAnnotated object saved to:", out_rds, "\n")

cat("\n Done. All plots in:", plot_dir, "\n")



# export cell type annottion table
annotation_table <- so@meta.data |>
  mutate(cluster_id = paste0(age, '_', sex, '_', seurat_clusters)) |>
  count(cluster_id, lineage) |>
  slice_max(n, n = 1, by = cluster_id) |>
  select(cluster_id, lineage) |>
  arrange(cluster_id)
rownames(annotation_table) <- NULL
head(annotation_table)
write_tsv(annotation_table, 'ImmuneNoise/spleen_10X/data/cell_type_annotation_table.tsv')


