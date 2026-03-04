# this is analogous to assemble_TMS_df.R which brings FACS and droplet into the same shape for the stratification pipeline



#-------------------------

# in the TMS the results were saved in a tibble but in one case not for pansci, so i added the code here
bs_stats <- function(tissue, bs_path){
  quantile_perc <- c(.025, .975) # this was specified in the pipeline_server script
  tissue <- 'Liver'
  #bs_path <- "RegulatedNoiseHel/data/processed/pansci_liver/bootstrapping_tables/"
  bs_files <- list.files(bs_path, pattern = "\\.rds$", full.names = TRUE)


  for (f in bs_files) {

   #f <- bs_files[6]
    message("Reading: ", f)
    age_group <- str_extract(f, "(?<=_)[0-9]+m_(female|male)(?=\\.rds)")
    print(paste('age_sex_group = ', age_group))
    bootstrap_res <- readRDS(f)
    res_var_cl <- readRDS(paste0('RegulatedNoiseHel/data/processed/pansci_liver/res_var_tables/res_var_cl_Liver_', age_group, '.rds' ))

    stats <- bootstrap_res |> 
      group_by(gene, cluster) |> 
      summarize(mean_resvar_bs = mean(ResVar),
                perc.hvg = sum(boolean_hvg) / n(),
                quant_low = quantile(ResVar, probs = quantile_perc[1]), 
                quant_high = quantile(ResVar, probs = quantile_perc[2]))
    
    results_tmp_tibble <- left_join(res_var_cl, stats, by = c("gene" = "gene", "cluster" = "cluster")) |>
      mutate(tissue = tissue, age = age_group)

    if (f == bs_files[1]) {
      results_tibble <- results_tmp_tibble
    } else {
      results_tibble <- rbind(results_tibble, results_tmp_tibble)
    }
  saveRDS(results_tibble, file = paste0("ImmuneNoise/pansci/data/prepped_for_strat/tmp/results_tibble_", age_group, ".rds"))
  results_tibble[1:10, 4:11]

  }
}



# if you want to conserve the cluster and cell type annotation this has to be added to the results table
transfer_cell_type_annotation <- function(so_path) {

# list the seurat objects from the pipeline
  so_path <- "RegulatedNoiseHel/data/processed/pansci_liver/seurat_objects/"
  so_files <- list.files(so_path, pattern = "\\.rds$", full.names = TRUE)


  all_results <- vector("list", length(so_files))
  for (i in seq_along(so_files)) {
    f <- so_files[i]
    #f <- so_files[1]

    # get metadata from seurat object after regnoise pipeline
    age_group <- str_extract(f, "[0-9]+m_(female|male)") # more general regex
    #age_group <- '12m_female'

    so <- readRDS(f)
    meta_df <- so@meta.data
    rm(so)
    gc()
    meta_df <- as_tibble(meta_df)
    meta_df <- tibble::remove_rownames(meta_df)


    # get df with og cluster names and new seurat clusters
    meta_df_selected <- meta_df |>
      dplyr::select(age_sex, Genotype, cell_ontology_class, Sub_cell_type, seurat_clusters ) |>
      distinct() |>
      mutate(og_clusters = str_extract(Sub_cell_type, "(?<=-).*")) |>
      mutate(cell_ontology_class = as.character(cell_ontology_class),
             seurat_clusters = as.character(seurat_clusters)) 


     # get results tibble from regnoise pipeline with res_var and bootstrapping
    results_tibble <- readRDS(paste0("ImmuneNoise/pansci/data/prepped_for_strat/tmp/results_tibble_", age_group, ".rds"))

    # match via seurat clusters
    results_merged <- results_tibble |>
      left_join(meta_df_selected, by = c(cluster = 'seurat_clusters', age = 'age_sex')) |>
      filter(cell_ontology_class == 'Hepatocytes') |> # REMOVE (select only hepatocytes)
      mutate(cell_ontology_class = 'hepatocyte') 
   
    # save
    write_tsv(results_merged, paste0("ImmuneNoise/pansci/data/prepped_for_strat/results_merged_", age_group, ".tsv"))
    all_results[[i]] <- results_merged
  }
  results_merged_all <- bind_rows(all_results)

   write_tsv(results_merged_all, ("ImmuneNoise/pansci/data/prepped_for_strat/results_with_og_clusters.tsv"))

   # remove og cluster info 
   main_result <- results_merged_all |>
    select(-Sub_cell_type, -og_clusters)  |>
    distinct() |>
    select(gene, age, cluster, gmean, res_var, cell_ontology_class, tissue, everything())

  write_tsv(main_result, ("ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages.tsv"))

}


annotate_gene_names <- function(){

  file <- read_tsv("ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages.tsv")
  df_clean <- file |>  mutate(ensembl_clean = sub("\\.\\d+$", "", gene))

  # specify gencode version (M27 is the one used in Pansci)
  gtf <- import("ImmuneNoise/pansci/data/gencode.vM27.annotation.gtf")

  gtf_map_raw <- mcols(gtf)[, c("gene_id", "gene_name")] |>
    as.data.frame() |>
    transmute(
      gene_id,
      ensembl_clean = sub("\\.\\d+$", "", gene_id),
      SYMBOL = gene_name) |>
    distinct(ensembl_clean, .keep_all = TRUE)

  # flag if different symbols have same gene names
  symbol_count <- gtf_map_raw |>
    distinct(ensembl_clean, SYMBOL) |>
    count(SYMBOL, name = "n_ens") 

  df_annotated <- df_clean |>
    left_join(gtf_map_raw, by = "ensembl_clean") |>
    left_join(symbol_count, by = "SYMBOL") |>
    mutate(
      gene = case_when(
        !is.na(SYMBOL) & SYMBOL != "" & n_ens == 1 ~ SYMBOL, 
        TRUE ~ gene) ) |> 
    select(-SYMBOL, -n_ens, -gene_id, -ensembl_clean) |>
    arrange(gene, age, cluster, gmean, res_var, cell_ontology_class, tissue, hvg, mean_resvar_bs) 
head(df_annotated)


  # optional: check how many % could not be mapped
  df_annotated_check <- df_annotated |>
    mutate(is_ensembl = str_detect(gene, "^ENSMUSG")) |>
    summarise(
        total_genes = n(),
        ensembl_genes = sum(is_ensembl),
        pct_ensembl = round(100 * ensembl_genes / total_genes, 2))
  df_annotated_check

  write_tsv(df_annotated, "ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages_annotated.tsv")
}




# this is the quick fix version but it was not used in the end because you cannot chose the same ensembl version as them
annotate_gene_names_old <- function(){


  file <- read_tsv("ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages.tsv")

  df_clean <- file |>
  mutate(ensembl_clean = sub("\\.\\d+$", "", gene))

  map_df <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = unique(df_clean$ensembl_clean),
    keytype = "ENSEMBL",
    columns = "SYMBOL") |>
    distinct(ENSEMBL, .keep_all = TRUE)

  df_annotated <- df_clean |>
    left_join(map_df, by = c("ensembl_clean" = "ENSEMBL")) |>
    mutate(gene = if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, gene) ) |>
    dplyr::select(-ensembl_clean, -SYMBOL, -og_clusters, -Sub_cell_type) |>
    distinct() |>
    arrange(cluster, gene)

  write_tsv(df_annotated, "ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages_annotated.tsv")
}





assemble_easysci_df <- function(write = FALSE){

  # load data
  print("Loading data ...")
  gene_expression_data <- read_tsv(("ImmuneNoise/pansci/data/prepped_for_strat/results_merged_all_ages_annotated.tsv"))
  print("Data loaded.")

  # prep df
  gene_expression_data_renamed <- gene_expression_data |>
    mutate(
    # no spaces any more
    cluster_id = paste(tissue, age, cluster, sep = "_"),
    organ   = stringr::str_extract(cluster_id, "^[A-Za-z]+"),
    months  = stringr::str_extract(cluster_id, "(?<=_)[0-9]+(?=m_)"),
    sex     = stringr::str_extract(cluster_id, "(?<=[0-9]m_)[A-Za-z]+"),
    sex     = stringr::str_to_lower(sex),
    cluster = stringr::str_extract(cluster_id, "(?<=_)[0-9]+$"),
    months  = paste0(as.integer(months), "m"),
    cluster_id = paste(organ, months, sex, cluster, sep = "_"))
head(gene_expression_data_renamed)


  # tag LVGs
  df_incl_lvgs <- gene_expression_data_renamed |>
    mutate(lvg = if_else(res_var < 1 & gmean != 0, TRUE, FALSE)) |>
    mutate(hvg = if_else(res_var >= 5 & gmean != 0, TRUE, FALSE))
  head(df_incl_lvgs)

  # select necessary columns
  df_finished <- df_incl_lvgs |>
    mutate(age = as.integer(str_extract(age, "\\d+(?=m)"))) |>
    dplyr::rename(perc_hvg = perc.hvg,
           cell_type = cell_ontology_class) |>
           mutate(cell_type = str_to_lower(cell_type)) |>
    dplyr::select(gene, gmean, cluster_id, cell_type, tissue, age, res_var, perc_hvg, hvg, lvg)
  print("Df is done.")
  #write <- TRUE

  if (write) {
      # write combined df to a csv
    write_csv(df_finished, paste0('ImmuneNoise/pansci/data/prepped_for_strat/combined_data.csv'))
    print("Saved df to data.")

   # get cell type names, write to excel for manual selection of innate immune cells
   # cell_type_list <- data.frame(unique(df_finished$cell_type))
   # write_xlsx(cell_type_list, pate0(out_path, '/cell_type_list.xlsx'))
    #print("Saved cell type list to data.")
  }

return(df_finished)
}


