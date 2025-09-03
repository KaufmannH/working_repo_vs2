


# filter and prepare expression df (from P1_assemble_TMS_df.)

stratify_df <- function(df = NULL, data_source = NULL, cell_type_selection = NULL) {

  if (is.null(df)) { 
    print("No df given, loading saved df.")
      if (data_source == "facs") {
        print("Data source: FACS.")
        df <- readRDS("facs/data/gene_set_df.rds")
    } else if (data_source == "droplet") {
      print("Data source: Droplet.")
      df <- readRDS("droplet/data/gene_set_df.rds")
    } else {
      stop("Issue in data source or cell type renaming.")
    } 
  }
 
  if (!is.null(cell_type_selection)) {
    df <- df |> filter(cell_type %in% cell_type_selection)
  }

  df <- df |>
    group_by(cluster_name) |>
    mutate(
      gmean_rank = {
        tmp <- rep(NA_integer_, n())
        pos <- which(!is.na(gmean) & gmean > 0)
        tmp[pos] <- ntile(gmean[pos], 6)
        tmp },
      expression_bin = if_else(
        gmean == 0,"0", as.character(gmean_rank) )) |>
    ungroup() |> 
    mutate(category = paste(gene_variability, expression_bin, sep = " "))
    
  df$category <- ifelse(df$`Not_expressed`, "Not expressed 0", df$category)

  #df_renamed <- df # in if: rename before saving to prevent issues
  if (data_source == "facs") {
      print("Saved to: FACS.")
      #rename_cols <- endsWith(colnames(df), " gene")
      #colnames(df_renamed)[rename_cols] <- gsub(" ", "_", colnames(df_renamed)[rename_cols]) 
      saveRDS(df, "facs/data/strat_df.rds")
    } else if (data_source == "droplet") {
      print("Saved to: droplet.")
      #rename_cols <- endsWith(colnames(df_renamed), " gene")
      #colnames(df_renamed)[rename_cols] <- gsub(" ", "_", colnames(df)[rename_cols])
      saveRDS(df, "droplet/data/strat_df.rds")
    } else {
      stop("Issue when saving.")
    } 
  return(df)
}



# what is the range of expression in the expression bins?

plot_expression_bins <- function(df, data_source) {

  summary_table <- df
  # make all possible gene categories appear
  all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  full_grid <- expand_grid(
    category = all_cats,
    cluster_name = unique(summary_table$cluster_name),
    cell_type = unique(summary_table$cell_type))

  existing_combinations <- unique(summary_table[c("category", "cluster_name", "cell_type")])
  missing_combinations <- anti_join(full_grid, existing_combinations,
                                    by = c("category", "cluster_name", "cell_type"))
  artifical_genes <- missing_combinations |>
    mutate(
      gene = "placeholder_gene",
      gmean = 0,
      gene_count = 0)
  summary_table <- bind_rows(summary_table, artifical_genes)
  summary_table$category <- factor(summary_table$category, levels = all_cats)

  strat_box <- ggplot(summary_table, aes(x = category, y = gmean, fill = category)) +
    geom_boxplot(size = 0.2, outlier.alpha = 0.2) +
    facet_wrap(~paste(cell_type, cluster_name, sep = "\n"), scales = "free") +
      labs(
        title = "All genes expression groups",
        x = "Category",
        y = "Gmean",
        fill = "Category"
      ) +
    theme_classic() +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )


   if (data_source == "facs") {
        print("Saved to: FACS.")
       ggsave('facs/plots/3_m/Test_8/strat_expression_bins.png', plot = strat_box,  width = 10, height = 14 )
    } else if (data_source == "droplet") {
      print("Saved to: droplet.")
      ggsave('droplet/plots/3_m/Test_8/strat_expression_bins.png', plot = strat_box,  width = 10, height = 14 )
    } else {
      stop("Issue when saving.")
    } 
  
}
#plot_expression_bins <- plot_expression_bins(df)




# what is the range of LPS expression in the expression bins?

plot_expression_bins_lps <- function(df, data_source) {

  summary_table <- df |>
    filter(`LPS response gene` == TRUE) 
  
  all_cats <- c(
  "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")
  summary_table$category <- factor(summary_table$category, levels = all_cats)

  full_grid <- expand_grid(
    category = all_cats,
    cluster_name = unique(summary_table$cluster_name),
    cell_type = unique(summary_table$cell_type))

  existing_combinations <- unique(summary_table[c("category", "cluster_name", "cell_type")])
  missing_combinations <- anti_join(full_grid, existing_combinations,
                                    by = c("category", "cluster_name", "cell_type"))
  artifical_genes <- missing_combinations |>
    mutate(
      gene = "placeholder_gene",
      gmean = 0,
      gene_count = 0)
  summary_table <- bind_rows(summary_table, artifical_genes)

  strat_box_lps <- ggplot(summary_table, aes(x = category, y = gmean, fill = category)) +
    geom_boxplot(size = 0.2, outlier.alpha = 0.2) +
    facet_wrap(~paste(cell_type, cluster_name, sep = "\n"), scales = "free") +
    labs(
      title = "LPS genes expression groups",
      x = "Category",
      y = "Gmean",
      fill = "Category"
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )

   if (data_source == "facs") {
        print("Saved to: FACS.")
        ggsave('facs/plots/3_m/Test_8/strat_expression_bins_lps.png', plot = strat_box_lps,  width = 10, height = 14 )
    } else if (data_source == "droplet") {
      print("Saved to: droplet.")
      ggsave('droplet/plots/3_m/Test_8/strat_expression_bins_lps.png', plot = strat_box_lps,  width = 10, height = 14 )
    } else {
      stop("Issue when saving.")
    } 
  print(plot)
  dev.off()
}



# how many genes are in the expression bins? 

plot_expression_bin_numbers <- function(df, data_source) {

  summary_table <- df |>
    select(gene, cell_type, cluster_name, category, expression_bin) |>
    distinct() |>
    group_by(cluster_name, cell_type, expression_bin) |>
    summarise(gene_count = n(), category = first(category), .groups = "drop")

  plot <- ggplot(summary_table, aes(x = factor(expression_bin), y = gene_count, fill = factor(cluster_name))) +
    geom_col() +
    facet_wrap(~paste(cell_type, cluster_name, sep = "\n"), scales = "free") +
    labs(
      title = "All genes",
      x = "Category",
      y = "Number of genes per expression bin",
      fill = "Category"
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )

   
   if (data_source == "facs") {
        print("Saved to: FACS.")
        pdf("facs/plots/3_m/Test_8/strat_expression_numbers.pdf",  width = 12, height = 15 )
    } else if (data_source == "droplet") {
      print("Saved to: droplet.")
      pdf("droplet/plots/3_m/Test_8/strat_expression_numbers.pdf",  width = 12, height = 15 )
    } else {
      stop("Issue when saving.")
    } 
  print(plot)
  dev.off()
}
#plot_expression_bin_numbers(strat_df, "facs")





# how many genes are in each category?

plot_category_numbers <- function(df, data_source) {
 
  summary_table <- df |>
    select(gene, cell_type, cluster_name, category) |>
    distinct() |>
    group_by(cluster_name, cell_type, category) |>
    summarise(gene_count = n(), .groups = "drop")


  all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")
  summary_table$category <- factor(summary_table$category, levels = all_cats)
     

  full_grid <- expand_grid(
    category = all_cats,
    cluster_name = unique(summary_table$cluster_name),
    cell_type = unique(summary_table$cell_type))

  existing_combinations <- unique(summary_table[c("category", "cluster_name", "cell_type")])
  missing_combinations <- anti_join(full_grid, existing_combinations,
                                    by = c("category", "cluster_name", "cell_type"))
  artifical_genes <- missing_combinations |>
    mutate(
      gene = "placeholder_gene",
      gene_count = 0)
  summary_table <- bind_rows(summary_table, artifical_genes)
  summary_table$category <- factor(summary_table$category, levels = all_cats)

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

 
  plot <- ggplot(summary_table, aes(x = category, y = gene_count, fill = category)) +
    geom_col() +
    facet_wrap(~paste(cell_type, cluster_name, sep = "\n"), scales = "free") +
    labs(
      title = "All genes",
      x = "Category",
      y = "Number of genes per category",
      fill = "Category" ) +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )

   
   if (data_source == "facs") {
        print("Saved to: FACS.")
        pdf("facs/plots/3_m/Test_8/strat_bar_all.pdf",  width = 12, height = 15 )
    } else if (data_source == "droplet") {
      print("Saved to: droplet.")
      pdf("droplet/plots/3_m/Test_8/strat_bar_all.pdf",  width = 12, height = 15 )
    } else {
      stop("Issue when saving.")
    } 
  print(plot)
  dev.off()
}
#plot_category_numbers(strat_df, "facs")








# how many genes of each gene set are in each category? 

plot_gene_set_numbers <- function(df, data_source) {
 df <- strat_df
 data_source <- "facs"
  all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  df <- df |> 
    mutate(cluster_plot = paste(cell_type, cluster_name, sep = "\n")) #|>
  
  gene_set_order <- names(df)[grepl(" gene$", names(df))]

  counts_df <- df |>
    pivot_longer(
      cols      = ends_with(" gene"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    mutate(
      gene_set = factor(gene_set, levels = gene_set_order),
      category = factor(category, levels = all_cats)) |>
    count(gene_set, cluster_plot, category, name = "n_genes") |>
    complete( gene_set, cluster_plot, category, fill = list(n_genes = 1))


  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(counts_df, aes(category, n_genes, fill = category)) +
    geom_col() +
    ggh4x::facet_grid2(
      gene_set ~ cluster_plot,
      scales      = "free_y",
      independent = "y",
      axes        = "all",
      strip       = strip_split(position = c("right","top"))) +
    labs(x = "Expression category", y = "Number of genes") +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      strip.text.y.left  = element_blank(),
      strip.text.y.right = element_text(angle = 0),
      strip.background.y = element_blank(),
      legend.position    = "top")
      

  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_8/strat_bar_gene_sets.png", plot, width = 45, height = 45)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_8/strat_bar_gene_sets.png", plot, width = 45, height = 45)
  } else {
    stop("Issue when saving.")
  }    
}



# density to compare the dist of the gene sets
plot_gene_set_proportions <- function(df, data_source) {

  gene_set_order <- names(df)[grepl(" gene$", names(df))]

   all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  counts_df <- df |>
    pivot_longer(
      cols      = all_of(gene_set_order),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    mutate(
      gene_set    = factor(gene_set, levels = gene_set_order),
      category    = factor(category, levels = all_cats),
      cluster_plot = paste(cell_type, "\n", cluster_name)) |>
    count(gene_set, cluster_plot, category, name = "n_genes") |>
    group_by(gene_set, cluster_plot) |>
    mutate(density = n_genes / sum(n_genes)) |>
    ungroup() |>
    complete(gene_set, cluster_plot, category, fill = list(density = 0.0001))

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(counts_df, aes(category, density, fill = category)) +
    geom_col() +
    ggh4x::facet_grid2(
      gene_set    ~ cluster_plot,
      scales      = "free_y",   
      independent = "y",   
      axes        = "all",
      switch      = "y") +
    labs(
      x = "Expression category",
      y = "Fraction of genes") +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1),
      legend.position       = "none")

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_11/gene_set_proportions.png", plot, width = 45, height = 45)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_11/gene_set_proportions.png", plot, width = 45, height = 45)
    } else {
      stop("Issue when saving.")
    }   
}



gene_set_proportions_compact <- function(df, data_source) {

  gene_set_order <- names(df)[grepl(" gene$", names(df))]

   all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  counts_df <- df |>
    pivot_longer(
      cols      = all_of(gene_set_order),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    mutate(
      gene_set    = factor(gene_set, levels = gene_set_order),
      category    = factor(category, levels = all_cats),
      cluster_plot = paste(cluster_name,  "| ", cell_type )) |>
    count(gene_set, cluster_plot, category, name = "n_genes") 
   

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(counts_df, aes(y = n_genes, x = cluster_plot, fill = category)) +
    geom_col(position = "fill") +
    facet_wrap("gene_set") +
    labs(
      y = "Fraction of genes") +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 90, hjust = 1, size = 6),
      axis.title.x          = element_blank(),
      legend.position       = "top")

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_11/gene_set_proportions_compact.png", plot, width = 10, height = 10)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_11/gene_set_proportions_compact.png", plot, width = 10, height = 10)
    } else {
      stop("Issue when saving.")
    }   

}





rest <- function() {
  # correlation of LVG 1 genes across clusters

  df_lvg1 <- df_main |>
    filter(category == "LVG_1")

  cluster_genes <- df_lvg1 |>
    group_by(cluster_id) |>
    summarise(genes = list(unique(gene)), .groups = "drop")

  cluster_names <- cluster_genes$cluster_id
  n <- length(cluster_names)
  shared_genes_matrix <- matrix(0, nrow = n, ncol = n,
                                dimnames = list(cluster_names, cluster_names))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      genes_i <- cluster_genes$genes[[i]]
      genes_j <- cluster_genes$genes[[j]]
      shared_genes_matrix[i, j] <- length(intersect(genes_i, genes_j))
    }
  }

  col_fun <- colorRamp2(c(0, max(shared_genes_matrix)), c("white", "darkred"))

  png("facs/plots/3_m/Test_8/lvg1_corr_all.png", width = 1200, height = 1000)

  Heatmap(shared_genes_matrix,
          name = "Shared genes",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(shared_genes_matrix[i, j], x, y, gp = gpar(fontsize = 8))
          })

  dev.off()

  # there is a higher correlation between the tissues and sexes

  # is the correlation of lvg 1 genes that are LPSs stronger?


  df_lvg1_lps <- df_lps |>
    filter(category == "LVG_1") 

  cluster_genes <- df_lvg1_lps |>
    group_by(cluster_id) |>
    summarise(genes = list(unique(gene)), .groups = "drop")

  cluster_names <- cluster_genes$cluster_id
  n <- length(cluster_names)
  shared_genes_matrix <- matrix(0, nrow = n, ncol = n,
                                dimnames = list(cluster_names, cluster_names))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      genes_i <- cluster_genes$genes[[i]]
      genes_j <- cluster_genes$genes[[j]]
      shared_genes_matrix[i, j] <- length(intersect(genes_i, genes_j))
    }
  }

  col_fun <- colorRamp2(c(0, max(shared_genes_matrix)), c("white", "darkred"))

  png("facs/plots/3_m/Test_8/lvg1_corr_lps.png", width = 1200, height = 1000)

  Heatmap(shared_genes_matrix,
          name = "Shared genes",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(shared_genes_matrix[i, j], x, y, gp = gpar(fontsize = 8))
          })

  dev.off()








  rest_exploration <- function(df, cluster_name) {


    # there are clusters with high LVG_1 expression:
    high_var_cluster <- c("Lung_3m_male_14", "Lung_3m_male_12", "Lung_3m_male_17") #just lung
    high_var_cluster <- c("Lung_3m_male_14", "Lung_3m_male_12", "Lung_3m_male_17", "Kidney_3m_female_7", "Kidney_3m_male_11", "Marrow_3m_female_15","Marrow_3m_male_8")

    # all genes select these special clusters
    selected_LVG_1 <- df_main |>
      filter(cluster_id %in% high_var_cluster) |>
      filter(category == "LVG_1") 

    #write.csv(selected_LVG_1, 'data/strat_LVG1.csv') # just high var
    #write.csv(df_main, 'data/strat_LVG1_high_and_low.csv') # all in LVG 1


    # lps genes select special clusters
    selected_LVG_1_lps <- df_lps |>
      filter(cluster_id %in% high_var_cluster) |>
      filter(category == "LVG_1") 

    print(selected_LVG_1_lps, n = 100)


    # is there a big overlap of highly expressed genes in LVG 1 cluster and LPS response genes? 
    intersect <- intersect(unique(selected_LVG_1_lps$gene), unique(selected_LVG_1$gene))

    overlap <- length(intersect)
    all <- length(unique(selected_LVG_1$gene))
    perc <- overlap/all
    all
    perc
    # 5% Lung
    # 5.2% all




    # is it more than exprected statistically
    # hypergeometric test

    hyper_df <- df_main |>
      filter(cluster_id %in% high_var_cluster)

    all_clusters <- unique(hyper_df$cluster_id)
    query_genes <- intersect(unique(lps_df$gene), unique(df_main$gene))


    hypergeometric_test <- function(cluster_id, df, lps_response_genes) {
      # Filter data for the specific cluster
      df_cluster <- df[df$cluster_id == cluster_id, ]

      # get list of all genes
      all_genes_list <- df |>
          distinct(gene) |>
          pull(gene) 

      # get mouse innate genes in dataset
      mouse_innate_genes <- intersect(lps_response_genes, all_genes_list)

      # list of lvg_1s
      lvg_list <- df_cluster  |> 
          filter(category == "Intermediate_6") |> 
          distinct(gene) |> 
          pull(gene)



      # list of overlap between lvg and innate genes
      overlap_list <- intersect(mouse_innate_genes, lvg_list)


      # signature (overall pool) - category of interest = failures
      n <- length(all_genes_list) - length(lvg_list)
      # category of interest
      m <- length(lvg_list) 
      # number of draws
      k <- length(mouse_innate_genes) 
      # overlap 
      q <- 0:(length(overlap_list) - 1) 
      q2 <- length(overlap_list) 

      p_value <- phyper(q2 - 1, m, n, k, lower.tail = FALSE)
      
      return(p_value)
    }


    p_values_result <- sapply(all_clusters, function(cluster_id) {
      hypergeometric_test(cluster_id, hyper_df, query_genes)})
    p_values_result

    p_values_df <- data.frame(
      cluster_id = all_clusters,
      p_value = p_values_result,
      significant = p_values_result < 0.05
    )
    p_values_df


    # significant (lung)
    # siginificant (all)



    # negative control
    # is there a big overlap of highly expressed genes in LVG 1 cluster and house keeping genes? 

    hk_df <- read_lines("data/scHK_human.symbols.txt")


    df_hk <- df_main |> 
    mutate(hk = str_detect(gene, paste(hk_list, collapse = "|"))) |>
    filter(hk == TRUE)
    head(df_hk)

    # all genes select these special clusters
    selected_LVG_1 <- df_main |>
      filter(cluster_id %in% high_var_cluster) |>
      filter(category == "LVG_1") 
    unique(selected_LVG_1$gene)
    #write.csv(selected_LVG_1, 'data/strat_LVG1.csv')

    # lps genes select special clusters
    selected_LVG_1_hk <- df_hk |>
      filter(cluster_id %in% high_var_cluster) |>
      filter(category == "LVG_1") 

    print(selected_LVG_1_hk, n = 100)

    intersect <- intersect(unique(selected_LVG_1_hk$gene), unique(selected_LVG_1$gene))

    overlap <- length(intersect)
    all <- length(unique(selected_LVG_1_hk$gene))
    perc <- overlap/all
    all
    perc
    # 100%
  }
}






# here the not expressed are always extra
plot_gene_set_numbers_old <- function(df, data_source) {

  extracted_cols <- names(df)[endsWith(names(df), " gene")]
  #extracted_cols <- extracted_cols[extracted_cols != "Not_expressed_immune_response_gene"]


  if (data_source == "facs") {
      print("Saved to: FACS.")
      pdf("facs/plots/3_m/Test_8/strat_bar_gene_sets.pdf", width = 12, height = 15 )
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      pdf("droplet/plots/3_m/Test_8/strat_bar_gene_sets.pdf", width = 12, height = 15 )
    } else {
      stop("Issue when saving.")
    } 


  for (gene_set_name in extracted_cols) {
    print(gene_set_name)
    summary_table <- df |>
      filter(.data[[gene_set_name]] == TRUE) |>
      select(gene, cell_type, cluster_name, category) |>
      distinct() |>
      group_by(cluster_name, cell_type, category) |>
      summarise(gene_count = n(), .groups = "drop")

    # make all possible gene categories appear
    all_cats <- c(
    "LVG 0", "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", 
    "HVG 0", "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5",
    "Intermediate 0", "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", 
    "Not expressed 0")
    summary_table$category <- factor(summary_table$category, levels = all_cats)

    full_grid <- expand_grid(
      category = all_cats,
      cluster_name = unique(summary_table$cluster_name),
      cell_type = unique(summary_table$cell_type)
    )
    existing_combinations <- unique(summary_table[c("category", "cluster_name", "cell_type")])
    missing_combinations <- anti_join(full_grid, existing_combinations,
                                      by = c("category", "cluster_name", "cell_type"))
    artifical_genes <- missing_combinations |>
      mutate(
        gene = "placeholder_gene",
        gene_count = 0)
    summary_table <- bind_rows(summary_table, artifical_genes)
    print(nrow(summary_table))

    plot <- ggplot(summary_table, aes(x = (category), y = gene_count, fill = category)) +
      geom_col() +
      facet_wrap(~paste(cell_type, cluster_name, sep = "\n"), scales = "free") +
      labs(
        title = paste0("Gene set: ", gene_set_name),
        x = "Category",
        y = "Number of genes per category",
        fill = "Category"
      ) +
      theme_classic() +
      theme(
        legend.position = "top",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
      )

    print(plot)
  }
  dev.off()
}



# does not work
# needs long df
plot_gene_set_proportions_old <- function(long_df, data_source, cluster_name) {


  summary_table <- long_df |>
    filter(cluster_name == cluster_name) |>
    select(gene, cluster_name, cell_type, category, gene_set) |>
    distinct() |>
    group_by(category, gene_set) |>
    summarise(cluster_name = first(cluster_name), cell_type = first(cell_type), gene_count = n_distinct(gene), .groups = "drop") |>
    mutate(gene_count_density = gene_count / sum(gene_count)) |>
    ungroup()
dim(summary_table) # 171 rows


  # make all possible gene categories appear
  all_cats <- c(
  "LVG 0", "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", 
  "HVG 0", "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5",
  "Intermediate 0", "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", 
  "Not expressed 0")
  summary_table$category <- factor(summary_table$category, levels = all_cats)

  full_grid <- expand_grid(
    category = all_cats,
    cluster_name = unique(summary_table$cluster_name),
    cell_type = unique(summary_table$cell_type)
  )
  existing_combinations <- unique(summary_table[c("category", "cluster_name", "cell_type")])
  missing_combinations <- anti_join(full_grid, existing_combinations,
                                    by = c("category", "cluster_name", "cell_type"))
  artifical_genes <- missing_combinations |>
    mutate(
      gene = "placeholder_gene",
      gene_count = 0)
  summary_table <- bind_rows(summary_table, artifical_genes)
long_df

  
  density <- ggplot(summary_table, aes(x = category, y = gene_count_density, fill = gene_set)) +
    geom_col() +
    facet_grid(rows = vars(gene_set), scales = "free", labeller = label_value) +
    labs(
      title = cluster_name,
      x = "Category",
      y = "Scaled number of genes",
      fill = "Category"
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "bottom" )

extracted_cluster_list <- unique(summary_table$cluster_name)
pdf("facs/plots/3_m/Test_11/gene_set_proportions0.pdf",  width = 7, height = 14 )

for (name in extracted_cluster_list) {
  plot <- plot_gene_set_proportions(summary_table, name)
  print(plot)
}
dev.off()
}


