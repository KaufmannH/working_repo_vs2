
# merge TMS data to one csv
# add cluster id and number of cells in each cluster


# in the TMS the results were saved in a tibble but in one case, so i added the code here
bs_stats <- function(tissue, bs_path){
  quantile_perc <- c(.025, .975) # this was specified in the pipeline_server script
  tissue <- 'Spleen'
  age_group <- '2m'
  bootstrap_res <- readRDS("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_additional/bootstrapping_tables/bootstrap_resSpleen_2m_male.rds")
  res_var_cl <- readRDS("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_additional/res_var_tables/res_var_cl_Spleen_2m_male.rds")

  stats <- bootstrap_res |> 
    group_by(gene, cluster) |> 
    summarize(mean_resvar_bs = mean(ResVar),
              perc.hvg = sum(boolean_hvg) / n(),
              quant_low = quantile(ResVar, probs = quantile_perc[1]), 
              quant_high = quantile(ResVar, probs = quantile_perc[2]))
  
  results_tibble <- left_join(res_var_cl, stats, by = c("gene" = "gene", "cluster" = "cluster")) |>
    mutate(tissue = tissue, age = age_group)
head(results_tibble)

write_tsv(results_tibble, file = "/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_additional/results_tibble.tsv")
results_tibble[1:10, 4:11]

}



assemble_additional_spleen <- function(write = FALSE){

  # load data
  print("Loading data ...")
  gene_expression_data <- read_tsv("/home/hkaufm49/analyses/NoiseControlCenter/RegulatedNoiseJoint/data/processed/spleen_additional/results_tibble.tsv")
  print("Data loaded.")

head(gene_expression_data)
  # prep dfs
  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) 
  head(gene_expression_data)
  # every biological replicate and tissue have the same number of rows (genes)

  # tag LVGs
  df_incl_lvgs <- gene_expression_data |>
    mutate(lvg = if_else(res_var < 1, TRUE, FALSE)) |>
    mutate(hvg = if_else(res_var >= 5 & gmean != 0, TRUE, FALSE))
  colnames(df_incl_lvgs)

  # select necessary columns
  df_finished <- df_incl_lvgs |>
    mutate(age = as.integer(str_extract(age, "\\d+(?=m)"))) |>
    rename(perc_hvg = perc.hvg) |>
    mutate( cell_type = 'Spleenocytes')  |>
    select(gene, gmean, cluster_id, cell_type, tissue, age, res_var, perc_hvg, hvg, lvg)
  print("Df is done.")
 
  if (write) {
      # write combined df to a csv
    write.csv(df_finished, '/home/hkaufm49/analyses/Spleen_10X/data/combined_data.csv')
    print("Saved df to data.")

  }

return(df_finished)
}


