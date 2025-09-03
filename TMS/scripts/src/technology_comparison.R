
combine_technologies <- function(){
  # combine facs and droplet data with entropy etc. 

  df_facs_raw <- readRDS("facs//data/variability_dir_metadata.rds")
  df_droplet_raw <- readRDS("droplet/data/variability_dir_metadata.rds")

  df_facs <- df_facs_raw |> mutate(technology = 'facs')
  df_droplet <- df_droplet_raw |> mutate(technology = 'droplet') 

  bound_df <- bind_rows(df_facs, df_droplet)

  all_cats <- c(
  "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
  "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  gene_set_order <- names(bound_df)[grepl(" gene$", names(bound_df))]
  combined_df <- bound_df |>
    pivot_longer(
        cols      = all_of(gene_set_order),
        names_to  = "gene_set",
        values_to = "in_set") |>
      filter(in_set == 1) |>
      select(gene, gene_set, category, cluster_id, rao_q, technology, variability_direction) |>
      arrange(gene, category)

  return(combined_df)
}


plot_scatter_rao <- function(){
  # plot entropy and counts
  plot_df <- combined_df |>
    select(gene, technology, rao_q) |>
    distinct() |>
    pivot_wider(
      names_from = technology,
      values_from =  rao_q )

  hex_plot <- ggplot(plot_df, aes(x = droplet, y = facs)) +
    geom_hex() +
    scale_fill_viridis(option = "magma", direction = -1, trans = "log") + #limits = c(0, 300)
    geom_vline(xintercept = 4,  colour = "grey") + 
    geom_hline(yintercept = 4,  colour = "grey") + 
    labs(
        x = "Droplet",
        y = "FACS", 
        fill = "log(count)") +
    theme_classic()

  ggsave("comparison/plots/3_m/Test_14/scatter_combo_entropy.png", hex_plot, width = 7, height = 7)
}
# 

plot_scatter_var_dir <- function(){
  # plot variability direction and counts
  plot_df <- combined_df |>
    select(gene, variability_direction, technology) |>
    distinct() |>
    pivot_wider(
      names_from = technology,
      values_from =  variability_direction)

  hex_plot <- ggplot(plot_df, aes(x = droplet, y = facs)) +
    geom_hex() +
    scale_fill_viridis(option = "magma", direction = -1) +
    geom_vline(xintercept = 0,  colour = "grey") + 
    geom_hline(yintercept = 0,  colour = "grey") + 
    labs(
        x = "Droplet",
        y = "FACS", 
        fill = "Count") +
    theme_classic()

  ggsave("comparison/plots/3_m/Test_14/scatter_combo_var_direction.png", hex_plot, width = 7, height = 7)
}


  plot <- ggplot(plot_df, aes(x = facs, y = droplet)) +
    geom_point(size = 1, alpha = 0.5, color = 'orange') +
    labs(
      x = "Droplet",
      y = "FACS") +
    theme_classic()





hirarchical_clustering <- function(combined_df){
    # df to matrix
  cluster_matrix <- combined_df |>
    group_by(gene, category, technology) |>
    summarise(num_genes = n(), .groups = "drop") |>
    pivot_wider(names_from = category, values_from = num_genes, values_fill = 0) |>
    mutate(gene_tech = paste0(gene, "_", technology)) |>
    select(-technology, -gene) |>
    column_to_rownames("gene_tech") |>
    scale()
  head(cluster_matrix)

  # calc distance and clusters
  dist_mat <- dist(cluster_matrix, method = "euclidean")
  hclust_res <- hclust(dist_mat, method = "ward.D2")
  dynamic_clusters <- cutreeDynamic(dendro = hclust_res, 
                                    distM = as.matrix(dist_mat),
                                    deepSplit = FALSE, 
                                   cutHeight = 600)
  # fuse clusters to gene                                  
  clusters_df <- data.frame(
                  gene_tech = rownames(cluster_matrix),
                  heatmap_cluster = dynamic_clusters)

  # fuse clusters back to df
  cluster_df_annot <- as.data.frame(cluster_matrix) |>
    rownames_to_column(var = "gene_tech") |>
    left_join(clusters_df, by = 'gene_tech') |>
    mutate(gene_tech_copy = gene_tech) |>
    separate(gene_tech_copy, into = c("gene", "technology"), sep = "_") |>
    pivot_longer(
      cols = matches("^(LVG|HVG|Intermediate|Not expressed)"),
      names_to = "category",
      values_to = "vals" ) |>
    select(-vals)
unique(cluster_df_annot$heatmap_cluster)
# intermediate not there
  saveRDS(cluster_df_annot, paste0("comparison/data/cluster_annot_df.rds"))
  
}

#df <- hirarchical_clustering(combined_df)

 



heatmap_clustered <- function(combined_df) {

all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
   "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats)) |>
    arrange(col_category)
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  # load the cluster info
  combined_clustered_df <- readRDS(("comparison/data/cluster_annot_df.rds"))

  bin_df_big <- combined_df |>
    mutate(gene_tech = paste0(gene, "_", technology)) |>
    left_join(combined_clustered_df, by = c("gene_tech", "category", "technology")) |>
    mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    select(gene_tech, category, entropy_bin, rao_q, heatmap_cluster, technology) 
head(bin_df_big)

  # summary for matrix
  bin_df_small <- bin_df_big |>
    mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    group_by(gene_tech, category, heatmap_cluster) |>
    summarise(gene_count_clusters = n(), .groups = "drop") |> # entropy_bin = first(entropy_bin), 
    filter(!is.na(category)) |>
    arrange(heatmap_cluster, category) |>
    select(-heatmap_cluster)
head(bin_df_small)


  # pivot for heatmap
  heatmap_matrix <- bin_df_small |>
    pivot_wider(
      names_from = category,
      values_from = gene_count_clusters, 
      values_fill = 0) |>
    column_to_rownames("gene_tech") |>
    mutate(across(everything(), as.numeric)) |>
    as.matrix() 
  head(heatmap_matrix, 100)
  heatmap_matrix <- heatmap_matrix[, all_cats] 


cluster_annotation <- bin_df_big |>
  select(gene_tech, technology, heatmap_cluster) |>
  distinct(gene_tech, .keep_all = TRUE) |>
  filter(gene_tech %in% rownames(heatmap_matrix)) |>
  arrange(match(gene_tech, rownames(heatmap_matrix)))

row_annot <- rowAnnotation(
  Technology = cluster_annotation$technology,
  annotation_name_side = 'top',
  simple_anno_size_adjust = TRUE)

col_annot <- HeatmapAnnotation(
    Category = anno_simple(all_cats,
    col = col_vec),
    annotation_name_side = "left")


heaty <- Heatmap(
  heatmap_matrix,
  name = "Gene counts in category", 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = colorRamp2(range(heatmap_matrix, na.rm = TRUE), c("white", "darkblue")),
  show_row_names = FALSE,
  left_annotation = row_annot,
  bottom_annotation = col_annot)

    png("comparison/plots/3_m/Test_14/heatmap_3.png",  width = 1200, height = 1000, res = 150) 
    draw(heaty)
    dev.off()
}
#heatmap_clustered(combined_df)





heatmap_clustered_backup <- function(combined_df) {

all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
   "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  # load the cluster info
  conbined_clustered_df <- readRDS(("comparison/data/cluster_annot_df.rds"))

  bin_df_big <- combined_df |>
    left_join(conbined_clustered_df, by = c("gene", "technology", "category")) |>
    mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    mutate(gene_tech = paste0(gene, "_", technology)) |>
    select(gene_tech, category, entropy_bin, rao_q, cluster) 
head(bin_df_big)

  # summary for matrix
  bin_df_small <- bin_df_big |>
   mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    group_by(gene_tech, category) |>
    summarise(gene_count_clusters = n(), .groups = "drop") |> # entropy_bin = first(entropy_bin), 
    filter(!is.na(category)) |>
    arrange(category)
head(bin_df_small)

  # pivot for heatmap
  heatmap_matrix <- bin_df_small |>
    pivot_wider(
      names_from = category,
      values_from = gene_count_clusters, 
      values_fill = 0) |>
    column_to_rownames("gene_tech") |>
    as.matrix() 
  head(heatmap_matrix, 100)


cluster_annotation <- bin_df_big |>
  select(gene_tech, cluster) |>
  distinct(gene_tech, .keep_all = TRUE) |>
  filter(gene_tech %in% rownames(heatmap_matrix)) |>
  arrange(match(gene_tech, rownames(heatmap_matrix)))


  # reorder matrix
  ordered_genes <- cluster_annotation |>
    arrange(cluster) |>
    pull(gene_tech)
  heatmap_matrix <- heatmap_matrix[ordered_genes, ]

  cluster_annotation <- cluster_annotation %>%
    filter(gene_tech %in% ordered_genes) %>%
    arrange(match(gene_tech, ordered_genes))

row_annot <- rowAnnotation(
    Cluster = anno_simple(
    as.character(cluster_annotation$cluster) ),
    annotation_name_side = 'top',
    simple_anno_size_adjust = TRUE)


col_annot <- HeatmapAnnotation(
    Category = anno_simple(colnames(heatmap_matrix),
    col = col_vec[colnames(heatmap_matrix)]),
    annotation_name_side = "left")


  heaty <- Heatmap(
    heatmap_matrix,
    name = "Gene counts in category", 
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = colorRamp2(range(heatmap_matrix, na.rm = TRUE), c("white", "darkblue")),
    row_names_gp = gpar(fontsize = 0), 
    left_annotation = row_annot,
    bottom_annotation = col_annot)

    png("comparison/plots/3_m/Test_14/heatmap_3.png",  width = 1200, height = 1000, res = 150) 
    draw(heaty)
    dev.off()
}
#heatmap_clustered(combined_df)





heatmap_entropy_bins <- function(master_entropy_df, data_source) {

all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
   "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)


  bin_df_big <- master_entropy_df |>
    mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    select(gene, gmean, category, entropy_bin, rao_q) 

  # summary for matrix
  bin_df_small <- bin_df_big |>
   mutate(category = factor(category, levels = all_cats)) |>
    mutate(entropy_bin = cut(rao_q, breaks = c(0, 0.1, 1.3, 2.6, 3.6, 4.3, 5.5, Inf), right = FALSE)) |>
    select(gene, gmean, category, entropy_bin, rao_q) |>
    group_by(gene, category) |>
    summarise(gene_count_clusters = n(), .groups = "drop") # entropy_bin = first(entropy_bin),

  # pivot for heatmap
  heatmap_matrix <- bin_df_small |>
    pivot_wider(
      names_from = category,
      values_from = gene_count_clusters, 
      values_fill = 0) |>
    select(gene, all_of(all_cats)) |> 
    column_to_rownames("gene") |>
    as.matrix()
  head(heatmap_matrix, 100)


  # annotate the bins 
  annot_df <- bin_df_big |>
    select(gene, entropy_bin) |>
    filter(gene %in% rownames(heatmap_matrix)) |>
    distinct(gene, .keep_all = TRUE)  

  annot_df$entropy_bin <- factor(annot_df$entropy_bin, levels = sort(unique(annot_df$entropy_bin))  )

  gene_order <- annot_df |>
    arrange(entropy_bin) |>
    pull(gene)

  annot_df <- annot_df[match(gene_order, annot_df$gene), ] 
  heatmap_matrix <- heatmap_matrix[gene_order, ]

  entropy_levels <- levels(annot_df$entropy_bin)
  entropy_colors <- setNames(
    c("#111D4E", "#4F71BE", "#74B2E0", "#9ecae1", "#65BCAE", "#7EAB55", "#445E30"),
    entropy_levels)
    entropy_colors

  row_annot <- rowAnnotation(
    Entropy = anno_simple(
      as.character(annot_df$entropy_bin), 
      col = entropy_colors,),
    annotation_name_side = 'top',
    simple_anno_size_adjust = TRUE)

  col_annot <- HeatmapAnnotation(
    Category = anno_simple(all_cats, 
    col = col_vec),
    annotation_name_side = "left")


  heaty <- Heatmap(
    heatmap_matrix,
    name = "Gene counts in category", 
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = colorRamp2(range(heatmap_matrix, na.rm = TRUE), c("white", "darkblue")),
    row_names_gp = gpar(fontsize = 0), 
    left_annotation = row_annot,
    bottom_annotation = col_annot)

    if (data_source == "facs") {
      png("facs/plots/3_m/Test_14/heatmap_2.png", width = 1200, height = 1000, res = 150)
      draw(heaty)
      dev.off()
      print("Saved to: FACS.")
    } else if (data_source == "droplet") {
      png("droplet/plots/3_m/Test_14/heatmap_2.png", width = 1200, height = 1000, res = 150)
      draw(heaty)
      dev.off()
      print("Saved to: droplet")
    } else {
      stop("Issue when saving.")
    } 
}






# old plots
combine_num_genes <- function(){
    df_list <- list()
    data_source_list <- c("droplet", "facs")

    for (data in data_source_list) {
        temp_df <- read.csv(paste0(data, "/data/strat_df.csv"))
        temp_df$data_source <- data
        df_list[[data]] <- temp_df
    }
    final_df <- do.call(rbind, df_list)
    head(final_df)
    return(final_df)
}


plot_all_num_genes <- function(df){
    df <- final_df
    rename_cols <- endsWith(colnames(df), "_gene")
    colnames(df)[rename_cols] <- gsub("_", " ", colnames(df)[rename_cols])
    extracted_cols <- names(df)[endsWith(names(df), " gene")]
    extracted_cols <- extracted_cols[extracted_cols != "Not expressed immune response gene"]


    df <- df[df$cell_type == "Macrophage", ] # this should not be necessary
    df$cell_type <- paste(df$cell_type, paste0("(", as.character(df$data_source), ")"), sep = " ")
   

    #pdf("/home/hkaufm49/working_repo/TMS/comparison/plots/Test_8/comparison_num_genes.pdf", width = 12, height = 7 )
    data_source_list <- c("droplet", "facs")
    tissue_list <- unique(df$tissue)
    #tissue_list <- c("Lung")
    
    for (t in tissue_list){
        tissue_df <- df[df$tissue == t, ] 
        print(paste("Processing tissue:", t))

        tissue_folder <- file.path("/home/hkaufm49/working_repo/TMS/comparison/plots/Test_8", t)
        if (!dir.exists(tissue_folder)) {
        dir.create(tissue_folder, recursive = TRUE)
        }

        for (gene_set_name in extracted_cols) {
            print(paste("Gene set: ", gene_set_name))
            gene_df <-  tissue_df |> filter(.data[[gene_set_name]] == TRUE)

            for (data_s in data_source_list) {
                print(paste("Data source: ", data_s))
                source_df <- gene_df[gene_df$data_source == data_s, ]

                summary_table <- source_df |>
                    select(gene, cell_type, cluster_name, category) |>
                    distinct() |>
                    group_by(cluster_name, cell_type, category) |>
                    summarise(gene_count = n(), .groups = "drop")
                
                all_cats <- c(
                "LVG 0", "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", 
                "HVG 0", "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5",
                "Other 0", "Other 1",  "Other 2", "Other 3", "Other 4", "Other 5", 
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
                    gene_count = 0
                )
                summary_table <- bind_rows(summary_table, artifical_genes)

                if (nrow(summary_table) > 0) {
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
                            legend.position = "none",
                            strip.background = element_blank(),
                            strip.text.x = element_text(size = 10, face = "bold"),  
                            axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1), 
                            axis.text.y = element_text(size = 10),                  
                            axis.title.x = element_text(size = 12),  
                            axis.title.y = element_text(size = 12),  
                            plot.title = element_text(hjust = 0.5, size = 14, face = "bold") 
                        )
                    filename <- paste0(gsub(" ", "_", gene_set_name), "_", data_s, ".png")
                    filepath <- file.path(tissue_folder, filename)
                    
                    ggsave(filepath, plot, width = 12, height = 5, dpi = 300)
                    print(paste("Saved plot to:", filepath))
                 
                } else {
                        message(paste("Skipped plot with no data"))    
            }
            
          }
        }
    }

}




merge_plot_all_num_genes <- function(){
    base_path = "/home/hkaufm49/working_repo/TMS/comparison/plots/Test_8"


  tissue_folders <- list.dirs(base_path, recursive = FALSE)
  
  for (tissue_folder in tissue_folders) {
    tissue_name <- basename(tissue_folder)
    print(paste("Processing tissue:", tissue_name))
    
    png_files <- list.files(tissue_folder, pattern = "\\.png$", full.names = TRUE)

    gene_sets <- unique(gsub("_(droplet|facs)\\.png$", "", basename(png_files)))

    vertical_stack <- list()
    
    for (gene_set in gene_sets) {
      print(paste("Merging gene set:", gene_set))
      
      droplet_image_path <- file.path(tissue_folder, paste0(gene_set, "_droplet.png"))
      facs_image_path <- file.path(tissue_folder, paste0(gene_set, "_facs.png"))
      

      image_list <- list()
      
      if (file.exists(droplet_image_path)) {
        image_list <- append(image_list, image_read(droplet_image_path))
      }
      if (file.exists(facs_image_path)) {
        image_list <- append(image_list, image_read(facs_image_path))
      }
      

      if (length(image_list) > 0) {
        combined_image <- image_append(image_list, stack = FALSE)
        vertical_stack <- append(vertical_stack, combined_image)
      }
    }
    

    if (length(vertical_stack) > 0) {
      tissue_combined <- image_append(vertical_stack, stack = TRUE)
      

      output_filename <- file.path(base_path, paste0(tissue_name, "_combined.png"))
      image_write(tissue_combined, output_filename)
      
      print(paste("Saved combined tissue image to:", output_filename))
    }
  }

}