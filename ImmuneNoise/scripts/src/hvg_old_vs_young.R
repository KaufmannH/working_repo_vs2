


library(readr)
library(dplyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(msigdbr)
library(clusterProfiler)
library(patchwork)


# TF mouse list (from le_network decoupler)
tf_list <- read_lines('droplet/data/mouse_tfs_decoupler.txt')
female_3m <- read_tsv('droplet/data/assembled_df_annotated.tsv')
spleen_full <- read_tsv('droplet/data/assembled_df_annotated_spleen_full.tsv')



# plot p values across old age groups 
# --------------------------------------------
reactome <- msigdbr(species = "Mus musculus", collection = "C2",
                    subcollection = "CP:REACTOME") |>
  select(gs_name, gene_symbol)

ages_old <- c("18m", "21m", "24m", "30m")

spleen_age <- spleen_full |>
  mutate(age = stringr::str_extract(unit_id, "(?<=_)[0-9]+m(?=_)"))

# per-age universe sizes -> common target N (smallest across ages)
universe_sizes <- sapply(ages_old, function(a) {
  spleen_age |>
    filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
    filter(age %in% c("3m", a)) |>
    filter(mean_expr > 0) |>
    pull(gene) |> unique() |> length()
})
target_universe_n <- min(universe_sizes)

classify_fn <- function(x) {
  x <- tolower(x)
  case_when(
    str_detect(x, "immune|interferon|cytokine|inflamm|interleukin|complement|antigen|tnf|toll|innate|adaptive|mhc|chemokine|neutrophil|phagocy|leukocyte|antimicrobial|leishmania|immunoregulat|infection") ~ "Inflammatory response",

    str_detect(x, "cell cycle|mitotic|mitosis|m phase|meiosis|checkpoint|cyclin|chromosome|chromatid|kinetochore|centromere|spindle|anaphase|pre replicative|replication|repair|dna damage|excision|mismatch|double.strand|homologous recombination|p53|translesion|lesion|abasic|ap site|strand|telomere|synthesis of dna|\\batr\\b|\\batm\\b") ~ "Genome maintenance & proliferation",

    str_detect(x, "unfolded protein|upr|hsf1|chaperone|ire1|heat shock|proteostasis|proteasome|autophagy|stress|programmed cell death|apoptosis|senescence") ~ "Proteostasis / stress",

    str_detect(x, "gpcr|g alpha|g protein|rhodopsin|ligand binding|receptor tyrosine|ntrk|rho gtpase|signaling by|signalling|kinase") ~ "Signaling",

    TRUE ~ "Other")
}


# ---- per-age ORA ----
ora_for_age <- function(this_old) {
  set.seed(888)

  matched <- spleen_age |>
    filter(age %in% c("3m", this_old)) |>
    filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
    filter(hvg) |>
    transmute(gene, mean_expr,
              age_group = if_else(age == "3m", "young", "old")) |>
    filter(mean_expr > 0) |>
    mutate(raw_bin  = ggplot2::cut_interval(mean_expr, n = 100),
           expr_bin = if_else(mean_expr > 2, "high_tail", as.character(raw_bin)))

  targets <- matched |>
    count(expr_bin, age_group) |>
    pivot_wider(names_from = age_group, values_from = n, values_fill = 0) |>
    mutate(n_match = pmin(young, old))

  sampled <- matched |>
    inner_join(select(targets, expr_bin, n_match), by = "expr_bin") |>
    filter(n_match > 0) |>
    group_by(expr_bin, age_group) |>
    group_modify(~ slice_sample(.x, n = .x$n_match[1])) |>
    ungroup()

  young_genes <- sampled |> filter(age_group == "young") |> pull(gene) |> unique()
  old_genes   <- sampled |> filter(age_group == "old")   |> pull(gene) |> unique()

  classes <- tibble(gene = union(young_genes, old_genes)) |>
    mutate(hvg_class = case_when(
      gene %in% young_genes & gene %in% old_genes ~ "shared",
      gene %in% young_genes                       ~ "young_only",
      TRUE                                        ~ "old_only"))

  # per-age universe, downsampled to common N (foreground always retained)
  universe_pool <- spleen_age |>
    filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
    filter(age %in% c("3m", this_old)) |>
    filter(mean_expr > 0) |>
    pull(gene) |> unique()

  extra_pool     <- setdiff(universe_pool, classes$gene)
  n_extra        <- max(0, target_universe_n - length(classes$gene))
  universe_sized <- union(classes$gene,
                          sample(extra_pool, min(n_extra, length(extra_pool))))

  run_one <- function(cls) {
    genes <- classes |> filter(hvg_class == cls) |> pull(gene) |> unique()
    if (length(genes) < 5) return(NULL)
    res <- enricher(gene = genes, universe = universe_sized,
                    TERM2GENE = reactome,
                    pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = 15, maxGSSize = 500)
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    as_tibble(res@result) |> mutate(hvg_class = cls, .before = 1)
  }

  map_dfr(c("shared", "young_only", "old_only"), run_one) |>
    mutate(old_age = this_old, .before = 1)
}

ora_trajectory <- map_dfr(ages_old, ora_for_age)

#fold enrichment 
plot_df <- ora_trajectory |>
  separate(GeneRatio, c("gr_k", "gr_n"), sep = "/", convert = TRUE) |>
  separate(BgRatio,   c("bg_k", "bg_n"), sep = "/", convert = TRUE) |>
  mutate(fold_enrichment = (gr_k / gr_n) / (bg_k / bg_n),
         old_age = factor(old_age, levels = ages_old),
         ID = ID |>
           gsub("^REACTOME_", "", x = _) |>
           gsub("_", " ", x = _) |>
           stringr::str_to_sentence(),
         line_id = paste(hvg_class, ID, sep = " | "),
         functional_group = classify_fn(ID))


sig_lines <- plot_df |> filter(p.adjust < 0.005) |> distinct(line_id) |> pull(line_id)
plot_df   <- plot_df |> filter(line_id %in% sig_lines)

plot_df <- plot_df |>
  mutate(functional_group = factor(functional_group,
    levels = c("Inflammatory response",
               "Genome maintenance & proliferation",
               "Proteostasis / stress",
               "Signaling",
               "Other")))

group_cols <- c("Inflammatory response"               = "#E08214", 
                "Genome maintenance & proliferation"  = "#4499C4",
                "Proteostasis / stress"               = "#70AD47", 
                "Signaling"                           = "#9B0049",
                "Other"                               = "grey70")



# ---- panels ----
make_panel <- function(cls, title) {
  df  <- plot_df |> filter(hvg_class == cls)
  avg <- df |>
    group_by(functional_group, old_age) |>
    summarise(mean_fe = mean(fold_enrichment), .groups = "drop")

  ggplot(df, aes(x = old_age, y = fold_enrichment)) +
    geom_line(aes(group = line_id, colour = functional_group),
              linewidth = 0.5, alpha = 0.3) +
    geom_point(aes(colour = functional_group), size = 1.6, alpha = 0.3) +
    geom_line(data = avg,
              aes(y = mean_fe, group = functional_group, colour = functional_group),
              linewidth = 1.8) +
    geom_point(data = avg,
               aes(y = mean_fe, colour = functional_group),
               size = 4) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = group_cols, name = "Functional group") +
    theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = NULL, y = "Fold enrichment", title = title)
}

p_young  <- make_panel("young_only", "Enriched in young only")
p_shared <- make_panel("shared",     "Enriched in shared")
p_old    <- make_panel("old_only",   "Enriched in old only")

plot_trajectory <- (p_young | p_shared | p_old) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "ORA fold enrichment across old ages",
                  subtitle = "TMS Spleen macrophages") &
  scale_y_continuous(limits = c(0, max(plot_df$fold_enrichment)))


ggsave("droplet/plots/regnoise_vs2/ora_time.png", plot = plot_trajectory, width = 20, height = 10, dpi = 300)




plot_df |>
  filter(functional_group == "Genome maintenance & proliferation",
         hvg_class %in% c("old_only", "shared")) |>
  distinct(hvg_class, ID) |>
  arrange(hvg_class, ID) |>
  print(n = Inf)

pat <- "cell cycle|mitotic|mitosis|m phase|meiosis|checkpoint|cyclin|chromosome|chromatid|kinetochore|centromere|spindle|anaphase|pre replicative|replication|repair|dna damage|excision|mismatch|double.strand|homologous recombination|p53|translesion|lesion|abasic|ap site|strand|telomere|synthesis of dna|atr|atm"

stringr::str_extract("extracellular matrix organization", pat)




# do classic dot plot comparison
# ------------------------------------


# set old age groups
old_age_groups <- c('18m', '21m', '24m', '30m')

spleen_age <- spleen_full |>
  mutate(age = stringr::str_extract(unit_id, "(?<=_)[0-9]+m(?=_)"))

universe_sizes <- sapply(old_age_groups, function(a) {
  spleen_age |>
    filter(age %in% c("3m", a)) |>
    filter(mean_expr > 0) |>
    pull(gene) |> unique() |> length()
})
target_universe_n <- min(universe_sizes)


for (old_age in old_age_groups){


  # compare the hvgs in young and old
  spleen_age <- spleen_full |>
    mutate(age = stringr::str_extract(unit_id, "(?<=_)[0-9]+m(?=_)"))
  head(spleen_age)
  unique(spleen_age$age)

  female_young_hvgs <- spleen_age |>
      filter(age == "3m") |>
      filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
      filter(hvg) |>
      select(gene, mean_expr) |>
      mutate(age_group = 'young')
  female_young_hvgs


  female_old_hvgs <- spleen_age |>
      filter(age == old_age) |>
      filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
      filter(hvg) |>
      select(gene, mean_expr) |>
      mutate(age_group = 'old')
  female_old_hvgs


  # merge
  merged_df <- female_young_hvgs |>
    bind_rows(female_old_hvgs) |>
    filter(mean_expr > 0)
  merged_df


  # plot mean expression distribution

  plot  <- ggplot(merged_df, aes(x = mean_expr, fill = age_group)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 100) +
    scale_fill_manual(values = c("young" = "steelblue", "old" = "firebrick")) +
    labs(x = "Mean expression", y = "Count", fill = "Age group") +
    theme_classic()

  ggsave("droplet/plots/regnoise_vs2/histogram_hvg_comparison.png", plot = plot, width = 8, height = 5, dpi = 300)




  # match gmean of hvg and lvg

  set.seed(888)  

  cutoff <- 2  

  # bin the genes
  binned_df <- merged_df %>%
    mutate(
      raw_bin  = ggplot2::cut_interval(mean_expr, n = 100),
      expr_bin = if_else(mean_expr > cutoff, "high_tail", as.character(raw_bin)))

  # count values in bin
  bin_targets <- binned_df %>%
    count(expr_bin, age_group) %>%
    pivot_wider(names_from = age_group, values_from = n, values_fill = 0) %>%
    mutate(n_match = pmin(young, old))

  # 3. sample genes from each bin
  matched_df <- binned_df %>%
    inner_join(select(bin_targets, expr_bin, n_match), by = "expr_bin") %>%
    filter(n_match > 0) %>%
    group_by(expr_bin, age_group) %>%
    group_modify(~ slice_sample(.x, n = .x$n_match[1])) %>%
    ungroup()


  female_young_hvgs_binned <- matched_df |>
    filter(age_group == 'young') |>
    pull(gene) |>
    unique()
  female_young_hvgs_binned


  female_old_hvgs_binned <- matched_df |>
    filter(age_group == 'old') |>
    pull(gene) |>
    unique()
  female_old_hvgs_binned


  hvg_class_df <- tibble(gene = union(female_young_hvgs_binned, female_old_hvgs_binned)) |>
      mutate(hvg_class = case_when(
          gene %in% female_young_hvgs_binned & gene %in% female_old_hvgs_binned ~ "shared",
          gene %in% female_young_hvgs_binned                             ~ "young_only",
          TRUE                                                 ~ "old_only"))
  hvg_class_df


  write_tsv(hvg_class_df, 'droplet/data/spleen_youg_vs_old.tsv')








  #  ORA reactome
  # ---------------------------------------------------

  hvg_class_df <- read_tsv('droplet/data/spleen_youg_vs_old.tsv')

  universe_pool <- spleen_age |>
    filter(str_detect(cluster, regex("macrophage", ignore_case = TRUE))) |>
    filter(age %in% c("3m", old_age)) |>
    filter(mean_expr > 0) |>
    pull(gene) |> unique()

  extra_pool   <- setdiff(universe_pool, hvg_class_df$gene)
  n_extra      <- max(0, target_universe_n - length(unique(hvg_class_df$gene)))
  universe_all <- union(hvg_class_df$gene,
                        sample(extra_pool, min(n_extra, length(extra_pool))))

  cat("Universe size for", old_age, ":", length(universe_all), "\n")


  reactome <- msigdbr(species = "Mus musculus", collection = "C2",
                      subcollection = "CP:REACTOME") |>
    dplyr::select(gs_name, gene_symbol) 



  run_ora <- function(cls) {
    genes <- unique(hvg_class_df$gene[hvg_class_df$hvg_class == cls])   # foreground
    if (length(genes) < 5) return(NULL)

    res <- enricher(gene = genes, universe = universe_all,
                    TERM2GENE = reactome,
                    pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = 15, maxGSSize = 500)

    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    as_tibble(res@result) |> mutate(hvg_class = cls, .before = 1)}

  ora_all <- map_dfr(c("shared", "young_only", "old_only"), run_ora)

  ora_sig_reactome <- ora_all |>
    filter(p.adjust < 0.05) |>
    arrange(hvg_class, p.adjust) |>
    dplyr::select(hvg_class, ID, Description, GeneRatio, BgRatio,  p.adjust, Count, geneID)

  #ora_sig_reactome



  # plot
  top_terms <- ora_sig_reactome |>
    group_by(ID) |>
    summarise(best_p = min(p.adjust), .groups = "drop") |>
    slice_min(best_p, n = 30, with_ties = FALSE) |>
    pull(ID)

  plot_df <- ora_all |>
    filter(p.adjust < 0.05, ID %in% top_terms) |>
    separate(GeneRatio, c("k", "n"), sep = "/", convert = TRUE) |>
    mutate(
      gene_ratio = k / n,
      group = factor(hvg_class, levels = c("young_only", "shared", "old_only")),
      ID = ID |>
        gsub("^REACTOME_", "", x = _) |>
        gsub("_", " ", x = _) |>
        stringr::str_to_sentence())



  assign(paste0('plot', old_age),
    ggplot(plot_df, aes(x = group, y = reorder(ID, gene_ratio))) +
      geom_point(aes(size = Count, colour = p.adjust)) +
      scale_colour_gradient(low = "#8b0000", high = "#f4a582",
                            name = "adj. p", trans = "log10") +
      scale_size_continuous(name = "genes", range = c(1.5, 6)) +
      theme_bw(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
            axis.text.y = element_text(size = 13),
            panel.grid.minor = element_blank()) +
      labs(x = NULL, y = NULL,
          title = paste0("Comparison of 3 months vs ", old_age),
          subtitle = "TMS Spleen macrophages"))

}


# collect all plots 
combined_plot <- plot18m + plot21m + plot24m + plot30m


ggsave(paste0("droplet/plots/regnoise_vs2/ora_reactome_dotplot.png"), combined_plot, width = 20, height = 15)












