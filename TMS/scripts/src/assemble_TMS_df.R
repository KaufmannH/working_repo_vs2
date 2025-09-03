
# merge TMS data to one csv
# add cluster id and number of cells in each cluster

assemble_TMS_df_facs <- function(write = FALSE){
  # load data
  print("Loading data ...")
  gene_expression_data <- read.delim("facs/facs_raw_data/results_tibble_regnoise_facs_06_24.tsv", header = TRUE, sep = "\t")
  annotation_info_raw <- read.delim("facs/facs_raw_data/cluster_cell_types_match_final_facs.tsv", header = TRUE, sep = "\t")
  metadata_cell_numbers_raw <- read.delim("facs/facs_raw_data/cluster_cell_numbers_facs_final.tsv", header = TRUE, sep = "\t")
  print("Data loaded.")

  # prep dfs
  annotation_info <- annotation_info_raw |>
    mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters))

  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster))

  metadata_cell_numbers <- metadata_cell_numbers_raw |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) |>
    select(- age, -cluster, -tissue, -tissue_cluster)

  # join dfs
  numbers_joined <- gene_expression_data  |>
    left_join(metadata_cell_numbers, by = "cluster_id") |>
    select(-cell_ontology_class, -tissue) 
  colnames(numbers_joined)

  all_joined <- numbers_joined  |>
    left_join(annotation_info, by = c("cluster_id")) 
  colnames(all_joined)

  # tag LVGs
  df_incl_lvgs <- all_joined |>
    mutate(lvg = if_else(res_var < 1 & gmean != 0, TRUE, FALSE)) |>
    mutate(hvg = if_else(res_var >= 5 & gmean != 0, TRUE, FALSE))
  head(df_incl_lvgs)


  # select necessary columns
  df_finished <- df_incl_lvgs |>
    mutate(age = str_extract(age, "\\d+m")) |>
    mutate(age = as.numeric(str_remove(age, "m"))) |> # newly added, check if works? 
    rename(num_cells_in_cluster = num_cells,
           perc_hvg = perc.hvg,
           cell_type = cell_ontology_class) |>
    select(gene, gmean, cluster_id, cell_type, tissue, age, res_var, perc_hvg, hvg, lvg)
  colnames(df_finished)
  print("Df is done.")

  if (write) {
    # write combined df to a csv
    write.csv(df_finished, 'facs/data/combined_data.csv')
    print("Saved df to data.")

    # get cell type names, write to excel for manual selection of innate immune cells
    cell_type_list <- data.frame(unique(df_finished$cell_type))
    write_xlsx(cell_type_list, 'facs/data/cell_type_list.xlsx')
    print("Saved cell type list to data.")
  }

return(df_finished)
}



assemble_TMS_df_droplet <- function(write = FALSE){

  # load data
  print("Loading data ...")
  gene_expression_data <- read.delim("droplet/droplet_raw_data/full_results_10_22.tsv", header = TRUE, sep = "\t")
  annotation_info_raw <- read_excel("droplet/droplet_raw_data/manual_annotation.xlsx")
  metadata_cell_numbers_raw <- read.table('droplet/droplet_raw_data/combined_metadata_for_cell_numbers.txt', header = TRUE, sep = "\t") # there are cluster_ids that are not in the annotaiton info
  print("Data loaded.")


  # prep dfs
  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) 
  head(gene_expression_data)
  # every biological replicate and tissue have the same number of rows (genes)

  # prep dfs
  annotation_info <- annotation_info_raw  # cluster id exists already

  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster))

  metadata_cell_numbers <- metadata_cell_numbers_raw |> # this df contains one row per cell
    mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters)) |>
    group_by(cluster_id) |>
    summarise(num_cells_in_cluster = n(), .groups = 'drop') 


# join dfs
  numbers_joined <- gene_expression_data  |>
    left_join(metadata_cell_numbers, by = "cluster_id") |>
    select(-tissue) 
  colnames(numbers_joined)

  all_joined <- numbers_joined  |>
    left_join(annotation_info, by = c("cluster_id")) 
  colnames(all_joined)

  # tag LVGs
  df_incl_lvgs <- all_joined |>
    mutate(lvg = if_else(res_var < 1, TRUE, FALSE)) |>
    mutate(hvg = if_else(res_var >= 5 & gmean != 0, TRUE, FALSE))
  colnames(df_incl_lvgs)

  # select necessary columns
  df_finished <- df_incl_lvgs |>
    mutate(age = as.integer(str_extract(age, "\\d+(?=m)"))) |>
    rename(perc_hvg = perc.hvg,
           cell_type = manual_final) |>
           mutate(cell_type = str_remove(cell_type, "^(lymphocyte_|leukocyte_)")) |>
    select(gene, gmean, cluster_id, cell_type, tissue, age, res_var, perc_hvg, hvg, lvg)
  print("Df is done.")
 
  if (write) {
      # write combined df to a csv
    write.csv(df_finished, 'droplet/data/combined_data.csv')
    print("Saved df to data.")

    # get cell type names, write to excel for manual selection of innate immune cells
    cell_type_list <- data.frame(unique(df_finished$cell_type))
    write_xlsx(cell_type_list, 'droplet/data/cell_type_list.xlsx')
    print("Saved cell type list to data.")
  }

return(df_finished)
}


