

library(dplyr)
library(ggplot2)
library(scales)

paths <- list(
    og   = "/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_og_3m",
    clip = "/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_clip_3m"
)

load_resvar <- function(base_path, source_label) {
    files_rv <- list.files(paste0(base_path, "/res_var_tables"), pattern = "\\.rds$", full.names = TRUE)
    ids_rv <- sub("^.*_(\\d+m_(male|female)).*$", "\\1", basename(files_rv))
    bind_rows(Map(\(f, id) dplyr::mutate(readRDS(f), age_sex = id), files_rv, ids_rv)) |>
        mutate(source = source_label)
}

combined_df <- bind_rows(Map(load_resvar, paths, names(paths)))

# split into og vs clip
og_rv   <- combined_df |> filter(source == "og")   |> dpyr::select(gene, cluster, age_sex, res_var_og   = res_var, gmean_og   = gmean)
clip_rv <- combined_df |> filter(source == "clip") |> dpyr::select(gene, cluster, age_sex, res_var_clip = res_var, gmean_clip = gmean)

res_var_threshold <- 0.00005


# genes with very low res_var in og
low_in_og <- og_rv |> filter(res_var_og < res_var_threshold)

# match them in clip
comparison <- low_in_og |>
    inner_join(clip_rv, by = c("gene", "cluster", "age_sex")) |>
    mutate(delta = res_var_clip - res_var_og)


# summary stats
nrow(comparison)
summary(comparison$res_var_clip)

summary <- comparison |>
    summarise(
        n_total          = n(),
        still_below_0.0001 = sum(res_var_clip < 0.0001),
        below_0.001        = sum(res_var_clip < 0.001),
        below_1          = sum(res_var_clip < 1),
        above_1          = sum(res_var_clip >= 1) )


# what resvar do the ones have that were super low before? 

p <- ggplot(comparison, aes(x = res_var_og, y = res_var_clip, color = cluster)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "grey") +
    scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
    labs(
        x = "res_var (og)",
        y = "res_var (no clipping)",
        color = "Cluster",
        title = "Change of residual variance\n after not clipping the pearson residuals" ) +
    theme_classic(base_size = 20) +
    theme(
        axis.title   = element_text(size = 20),
        axis.text    = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 20)
    )

ggsave("/home/hkaufm49/analyses/working_repo/ImmuneNoise/comparison/plots/3_m/Test_14/lowresvar_og_vs_clip.png",   p, width = 12, height = 10)


# what are the genes that changed? 

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# 1. genes that crossed above 1 in clip (per cluster)
moved_up <- comparison |>
    filter(res_var_clip >= 1)

table(moved_up$cluster)  # see counts per cluster

# 2. background: all expressed genes (gmean > threshold) from og
gmean_threshold <- 0.01
background_genes <- og_rv |>
    filter(gmean_og > gmean_threshold) |>
    pull(gene) |>
    unique()

length(background_genes)

# 3. map symbols -> Entrez once for everything (clusterProfiler GO works on either,
#    but Entrez avoids ambiguity for some symbols)
all_symbols <- unique(c(moved_up$gene, background_genes))
gene_map <- bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

background_entrez <- gene_map |>
    filter(SYMBOL %in% background_genes) |>
    pull(ENTREZID) |>
    unique()

# 4. enrichment per cluster
clusters <- unique(moved_up$cluster)

enrich_one_cluster <- function(clust) {
    sym <- moved_up |> filter(cluster == clust) |> pull(gene) |> unique()
    entrez <- gene_map |> filter(SYMBOL %in% sym) |> pull(ENTREZID) |> unique()

    if (length(entrez) < 5) {
        message("Cluster ", clust, ": only ", length(entrez), " mapped genes, skipping")
        return(NULL)
    }

    enrichGO(
        gene          = entrez,
        universe      = background_entrez,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",         # BP / MF / CC / ALL
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )
}

enrich_results <- setNames(lapply(clusters, enrich_one_cluster), clusters)

# 5. compact table of top hits per cluster
top_per_cluster <- bind_rows(
    lapply(names(enrich_results), function(cl) {
        res <- enrich_results[[cl]]
        if (is.null(res) || nrow(res@result) == 0) return(NULL)
        res@result |>
            arrange(p.adjust) |>
            head(10) |>
            mutate(cluster = cl)
    })
)

top_per_cluster |> dpyr::select(cluster, ID, Description, p.adjust, Count, geneID)

# 6. quick visualization — top terms per cluster
p <- top_per_cluster |>
    group_by(cluster) |>
    slice_min(p.adjust, n = 5) |>
    ungroup() |>
    mutate(Description = reorder(Description, -p.adjust)) |>
    ggplot(aes(x = -log10(p.adjust), y = Description, fill = cluster)) +
    geom_col() +
    facet_wrap(~ cluster, scales = "free_y") +
    labs(x = "-log10(adj p)", y = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

ggsave("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/go_enrichment_moved_genes.png",
       p, width = 16, height = 12)


library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(data.table)

# 1. genes that crossed above 1 in clip (per cluster)
moved_up <- comparison |>
    filter(res_var_clip >= 1)

table(moved_up$cluster)  # see counts per cluster

# 2. background: all expressed genes (gmean > threshold) from og
gmean_threshold <- 0.01
background_genes <- og_rv |>
    filter(gmean_og > gmean_threshold) |>
    pull(gene) |>
    unique()

length(background_genes)

# 3. map symbols -> Entrez once for everything (clusterProfiler GO works on either,
#    but Entrez avoids ambiguity for some symbols)
all_symbols <- unique(c(moved_up$gene, background_genes))
gene_map <- bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

background_entrez <- gene_map |>
    filter(SYMBOL %in% background_genes) |>
    pull(ENTREZID) |>
    unique()

# 4. enrichment per cluster
clusters <- unique(moved_up$cluster)

enrich_one_cluster <- function(clust) {
    sym <- moved_up |> filter(cluster == clust) |> pull(gene) |> unique()
    entrez <- gene_map |> filter(SYMBOL %in% sym) |> pull(ENTREZID) |> unique()

    if (length(entrez) < 5) {
        message("Cluster ", clust, ": only ", length(entrez), " mapped genes, skipping")
        return(NULL)
    }

    enrichGO(
        gene          = entrez,
        universe      = background_entrez,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",         # BP / MF / CC / ALL
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )
}

enrich_results <- setNames(lapply(clusters, enrich_one_cluster), clusters)

# 5. compact table of top hits per cluster
top_per_cluster <- bind_rows(
    lapply(names(enrich_results), function(cl) {
        res <- enrich_results[[cl]]
        if (is.null(res) || nrow(res@result) == 0) return(NULL)
        res@result |>
            arrange(p.adjust) |>
            head(10) |>
            mutate(cluster = cl)
    })
)

top_per_cluster |> select(cluster, ID, Description, p.adjust, Count, geneID)

# 6. quick visualization — top terms per cluster
p <- top_per_cluster |>
    group_by(cluster) |>
    slice_min(p.adjust, n = 5) |>
    ungroup() |>
    mutate(Description = reorder(Description, -p.adjust)) |>
    ggplot(aes(x = -log10(p.adjust), y = Description, fill = cluster)) +
    geom_col() +
    facet_wrap(~ cluster, scales = "free_y") +
    labs(x = "-log10(adj p)", y = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

ggsave("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/go_enrichment_moved_genes.png",
       p, width = 16, height = 12)

