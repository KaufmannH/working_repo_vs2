# prepare dataframe for analysis: load and tag genes for all conditions to be tested downstream. 



load_filter_og_df <- function(data_source) {

  df_main_all <- read.csv(paste0(data_source, "/data/combined_data.csv"))

  # get lists of existing immune cell types
  if (data_source == "facs") {
    print("Data source: FACS.")
    immune_cell_type_list <- read_xlsx("facs/data/cell_type_list_selected.xlsx") |> 
      filter(!is.na(innate_cell)) |> 
      pull(innate_cell)

    # filter df
    df_main_filtered <- df_main_all |>
      filter(age == 3) |>
      mutate(cell_type = str_to_title(cell_type)) |>
      filter(cell_type == "Macrophage")  # or %in% immune_cell_type_list

    # save
    write.csv(df_main_filtered, file = "facs/data/df_main_filtered.csv", row.names = FALSE)
    
  } else if (data_source == "droplet") {
    print("Data source: Droplet.")
    immune_cell_type_list <- read_xlsx("droplet/data/cell_type_list_selected.xlsx") |> 
      filter(!is.na(innate_cell)) |> 
      pull(innate_cell)
    
    # filter df
    df_main_filtered <- df_main_all |>
      filter(age == 3) |>
      mutate(cell_type = str_to_title(cell_type)) |>
      filter(cell_type == "Macrophage")  # or %in% immune_cell_type_list
    # save
    write.csv(df_main_filtered, file = "droplet/data/df_main_filtered.csv", row.names = FALSE)

  } else {
    stop("Issue in data source or cell type renaming.")
  }

  return(df_main_filtered)
}



tag_imm_dict_genes <- function(data_source){
#data_source <- "facs"
  df <- read.csv(paste0(data_source, "/data/df_main_filtered.csv"))
  imm_dict_gene_list <- readRDS("reference_gene_sets/ImmuneDict/marker_genes_condition.rds")
  imm_dict_gene_names <- names(imm_dict_gene_list)

  for (condition in imm_dict_gene_names){
      print(condition)
  
    df <- df |>
      mutate(!!paste0(condition, "_gene") := gene %in% imm_dict_gene_list[[condition]])

    # tag genes not in df
    all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
    all_combinations <- crossing(all_clusters, gene = imm_dict_gene_list[[condition]])
    missing_genes <- anti_join(all_combinations, df, 
                            by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
    missing_genes <- missing_genes |>
    mutate( !!paste0(condition, "_gene") := TRUE,
            Not_expressed  = TRUE,
            gmean = 0)
    df <- bind_rows(df, missing_genes)
  }
  return(df)
}



# add and tag housekeeping genes
tag_hk_genes <- function(data_source) {

 if (data_source == "facs") {
  df <- read.csv(paste0(data_source, "/data/df_main_filtered.csv"))
 } else if (data_source == "droplet") {
    df <- read.csv(paste0(data_source, "/data/df_main_filtered.csv"))
 } else {
   print("Issue loading data.")
 }
  # initiate the col
  df$Not_expressed <- FALSE

  # general house keeping genes
  housekeeping_df <- read.csv('reference_gene_sets/Housekeeping_TranscriptsMouse.csv', sep = ";") 
  housekeeping_list <- housekeeping_df |>
      pull(Genes)

  df <- df |>
    mutate(Housekeeping_gene = gene %in% housekeeping_list)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = housekeeping_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Housekeeping_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}



tag_hk_lin_genes <- function(df) {
  # housekeeping from Lin et al. 2019
  ## https://academic.oup.com/gigascience/article/8/9/giz106/5570567
  housekeeping_lin_list <- read_lines("reference_gene_sets/scHK_human.symbols.txt")

  # load the mouse/human conversions
  m2h <- read.delim("reference_gene_sets/20200307_ensembl/mouse.txt", header = FALSE, sep = "\t")
  colnames(m2h) <- c("mouse", "human")
  h2m <- read.delim("reference_gene_sets/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
  colnames(h2m) <- c("human", "mouse")
  head(h2m)

  mouse_lin <- h2m |>
    filter(human %in% housekeeping_lin_list) |>
    pull(mouse) |>
    unique()
  length(mouse_lin) # 1026

  # add and tag housekeeping genes
  df <- df |>
      mutate(Housekeeping_Lin_gene = gene %in% mouse_lin)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = housekeeping_lin_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Housekeeping_Lin_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}




# add and tag innate response genes 
tag_innate_response_genes <- function(df) {

  #load the innate immunity genes of mice database # 646 genes
  mouse_innate_genes_df_raw <- read.csv('reference_gene_sets/mouse_innate_genes.csv')
  mouse_innate_genes <- mouse_innate_genes_df_raw |> 
      select(Gene.Symbol) |>
      rename(gene = Gene.Symbol) |>
      unique() |>
      pull(gene)

  df <- df |>
     mutate(Immune_response_gene = gene %in% mouse_innate_genes)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = mouse_innate_genes)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Immune_response_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  #create one response gene for every cell type
  #selection <- df_main_filtered |>
  #  select(cluster_id, cell_type) |>
  #  distinct()
  #mouse_innate_genes <- crossing(
  #selection,
  #gene = mouse_innate_genes$gene)

  # add to main df
  #df <- df |>
    #full_join(mouse_innate_genes, by = c("gene", "cluster_id", "cell_type")) |>
   
   # mutate(Not_expressed_immune_response_gene = is.na(gmean)) |>
    #mutate(Immune_response_gene = !Not_expressed_immune_response_gene) |>
   # mutate(across(c( hvg, lvg), ~replace_na(., FALSE))) |>
   # mutate(gmean = ifelse(Not_expressed_immune_response_gene, 0, gmean)) 

  return(df)
}



 # probably not functional
  # add markers: get Marker Genes
tag_marker_genes <- function(df) {
  markers_raw  <- read.delim("reference_gene_sets/lineage_specific_markers.csv", header = TRUE, sep = "\t")
  available_cell_types_big = c('Alveolar macrophages', 'Kupffer cells', 'Macrophages', 'Langerhans cells', 'Microglia', 'Monocytes', 'Basophils', 'Dendritic cells', 'Eosinophils', 'NK cells', 'Neutrophils')
  available_cell_types_small = c('Dendritic cells','Macrophages','Monocytes', 'Alveolar macrophages') 

  markers_df <- markers_raw |>
    filter(grepl("Mm", species)) |>
    rename(gene = official.gene.symbol) |>
    rename(cell_type = cell.type) |>
    mutate(gene = str_to_lower(gene),           
          gene = str_replace(gene, "^(.)", str_to_upper)) |>
    filter(cell_type %in% available_cell_types_small) |>
    select(gene, cell_type) |>
    mutate(cell_type = case_when(
      cell_type == "Macrophages" ~ "Macrophage",
      cell_type == "Alveolar macrophages" ~ "Macrophage",
      cell_type == "Monocytes" ~ "Monocyte",
      cell_type == "Dendritic cells" ~ "DC")) |>
    unique() |>
    #group_by(cell_type) |>
    #summarise(count = n()) |>
    mutate(is_marker = TRUE)
  head(markers_df)
  #write.csv(markers_df, "data/lineage_specific_markers.csv", row.names = FALSE)

  # get the lists extra
  macrophage_unique <- markers_df |>
    filter(cell_type == 'Macrophage') |>
    distinct(gene) |>
    pull(gene)

  monocyte_unique <- markers_df |>
    filter(cell_type == 'Monocyte') |>
    distinct(gene) |>
    pull(gene)

  dc_unique <- markers_df |>
    filter(cell_type == 'DC') |>
    distinct(gene) |>
    pull(gene)


  macrophage_unique <- setdiff(macrophage_unique, union(monocyte_unique, dc_unique))
  monocyte_unique <- setdiff(monocyte_unique, union(macrophage_unique, dc_unique))
  dc_unique <- setdiff(dc_unique, union(macrophage_unique, monocyte_unique))

  length(intersect(all_genes, macrophage_unique))
  length(intersect(all_genes, monocyte_unique))
  length(intersect(all_genes, dc_unique))

  # make new df markers
  macrophage_df <- data.frame(gene = macrophage_unique, cell_type = "Macrophage", is_marker = TRUE)
  monocyte_df <- data.frame(gene = monocyte_unique, cell_type = "Monocyte", is_marker = TRUE)
  dc_df <- data.frame(gene = dc_unique, cell_type = "DC", is_marker = TRUE)
  markers_df <- bind_rows(macrophage_df, monocyte_df, dc_df)


  # add to main df
  df_main_markers <- df |>
    left_join(markers_df, by = c("gene", "cell_type")) |>
    mutate(Marker_gene = ifelse(is.na(is_marker), FALSE, is_marker))
  head(df_main_markers)

  # check outcome
  marker_summary <- df_main_markers |> 
    group_by(cell_type, is_marker) |>
    summarise(count = n(), .groups = 'drop')
  marker_summary

  return(df_main_markers)
}
#df_marker_genes <- tag_marker_genes(df_hk_lin_genes)

tag_mac_marker_genes <- function(df) {
 mac_marker_genes <- c('Fcgr1', 'Mertk',  'Spi1', 'Maf', 'Mafb', 'Cebpa')

  df <- df |>
     mutate(Macrophage_marker_gene = gene %in% mac_marker_genes)
  
  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = mac_marker_genes)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Macrophage_marker_gene = TRUE,
          Not_expressed  = TRUE, 
          gmean = 0)
  df <- bind_rows(df, missing_genes)
  
  return(df)
}


# load LPS stimulated genes in monocytes
tag_lps_genes <- function(df) {

  lps_df <- read_excel("reference_gene_sets/Bhatt_2012_data.xlsx")
  colnames(lps_df)[colnames(lps_df) %in% c("Probe...1", "Probe...2", "Promoter Class", "...4")] <- c("gene", "gene_id", "promoter_class", "group")
  lps_df_early <- lps_df |> filter(group <= 5)
  lps_df_late <- lps_df |> filter(group > 5)
  lps_gene_list_early <- lps_df_early |> pull(gene)
  lps_gene_list_late <- lps_df_late |> pull(gene)

  # number of genes in there
  length(unique(lps_df$gene))
  length(unique(lps_df_early$gene))
  length(unique(lps_df_late$gene))

  df <- df |> 
  # dont ever do that again!
  # mutate(LPS_response_early_gene = str_detect(gene, paste(lps_gene_list_early, collapse = "|"))) |>
  # mutate(LPS_response_late_gene = str_detect(gene, paste(lps_gene_list_late, collapse = "|"))) |>  
  # mutate(lps_stim = LPS_response_early_gene | lps_stim_late) 
    mutate(
      LPS_response_early_gene = gene %in% lps_gene_list_early,
      LPS_response_late_gene = gene %in% lps_gene_list_late,
      LPS_response_gene = LPS_response_early_gene | LPS_response_late_gene)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  # early
  all_combinations <- crossing(all_clusters, gene = lps_gene_list_early)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( LPS_response_early_gene = TRUE,
          LPS_response_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)
# late
 all_combinations <- crossing(all_clusters, gene = lps_gene_list_late)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( LPS_response_late_gene = TRUE,
          LPS_response_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)


  return(df)  
}

# load chemokine sinalling pathway
tag_chemokine_genes <- function(df) {
  cytokine_df <- read.csv("reference_gene_sets/InnateDB_cytokine_signalling.csv")
  cytokine_gene_list <- cytokine_df |>
    pull(name) |>
    unique()

  df <- df |>
     mutate(Cytokine_response_gene = gene %in% cytokine_gene_list)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = cytokine_gene_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Cytokine_response_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df) 
}


# load TLR pathway
tag_tlr_genes <- function(df) {
  tlr_df <- read.csv("reference_gene_sets/InnateDB_genes_TLR.csv")
  tlr_gene_list <- tlr_df |>
    pull(name) |>
    unique()

  df <- df |>
    mutate(TLR_response_gene = gene %in% tlr_gene_list)

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = tlr_gene_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( TLR_response_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}



# load xue et al. response genes
tag_xue_genes <- function(df) {
  xue_response_raw <- readRDS("reference_gene_sets/xue_response_genes.rds")
  xue_response <- xue_response_raw |> select(Condition, Upregulated)

  # load the mouse/human conversions
  m2h <- read.delim("reference_gene_sets/20200307_ensembl/mouse.txt",  header = FALSE, sep = "\t")
  colnames(m2h) <- c("mouse", "human")
  h2m <- read.delim("reference_gene_sets/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
  colnames(h2m) <- c("human", "mouse")

  # see if the gene is in def and then pull the ortholog
  xue_response_orthologs <- xue_response |>
    mutate(Orthologs = map(Upregulated, ~ h2m$mouse[match(.x, h2m$human)]))

  xue_gene_list <- xue_response_orthologs |> 
    select(Condition, Orthologs) |> 
    unnest(Orthologs) |> 
    filter(!is.na(Orthologs)) |>    
    group_by(Condition) |> 
    summarise(GeneList = list(unique(Orthologs))) |> 
    deframe()


  names(xue_gene_list)[names(xue_gene_list) == "IFNg"] <- "IFNy_response_gene"
  names(xue_gene_list)[names(xue_gene_list) == "IL10"] <- "IL10_response_gene"
  names(xue_gene_list)[names(xue_gene_list) == "IL4"] <- "IL4_response_gene"
  names(xue_gene_list)[names(xue_gene_list) == "TNF"] <- "TNF_response_gene"
  length(xue_gene_list[["IFNy_response_gene"]])
  length(xue_gene_list[["IL10_response_gene"]])
  length(xue_gene_list[["IL4_response_gene"]])
  length(xue_gene_list[["TNF_response_gene"]])


  df <- df |>
    mutate(IFNy_response_gene = gene %in% xue_gene_list[["IFNy_response_gene"]]) |>
    mutate(IL10_response_gene = gene %in% xue_gene_list[["IL10_response_gene"]]) |>
    mutate(IL4_response_gene = gene %in% xue_gene_list[["IL4_response_gene"]]) |>
    mutate(TNF_response_gene = gene %in% xue_gene_list[["TNF_response_gene"]]) 

  # tag genes not in df
  for (gene_set_name in names(xue_gene_list)) {
    gene_list <- xue_gene_list[[gene_set_name]]
    
    all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
    all_combinations <- crossing(all_clusters, gene = gene_list)

    missing_genes <- anti_join(all_combinations, df, 
                              by = c("cluster_id", "cell_type", "gene", "tissue", "age"))

    missing_genes <- missing_genes |>
      mutate(
        !!gene_set_name := TRUE,
        Not_expressed = TRUE,
        gmean = 0)
    df <- bind_rows(df, missing_genes)
  }
  return(df)
}




tag_autoimmune_genes <- function(df){
  # add autimmune pleiotropic loci
  autimmune_df <- read_excel('reference_gene_sets/Marquez_pleiotropic_genes.xlsx')
  autoimmune_list <- autimmune_df |> pull(Gene) |> unique()

  # load the mouse/human conversions
  m2h <- read.delim("reference_gene_sets/20200307_ensembl/mouse.txt",  header = FALSE, sep = "\t")
  colnames(m2h) <- c("mouse", "human")
  h2m <- read.delim("reference_gene_sets/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
  colnames(h2m) <- c("human", "mouse")

  # conversion from human to mouse
  mouse_auto <- h2m |>
    filter(human %in% autoimmune_list) |>
    pull(mouse) |> 
    unique()

  df <- df |>
    mutate(Autoimmunity_gene = gene %in% mouse_auto) 
  # here i start with 34 human genes, which gives me 40 mouse genes of which 34 were detected in the atlas

   # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = mouse_auto)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Autoimmunity_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}



# add negative control gene sets

# sperm DNA condensation
tag_sperm_genes <- function(df) {
  sperm_table <- read_excel('reference_gene_sets/GO_gene_lists/GO_sperm_dna_condensation.xlsx')
  sperm_list <- sperm_table |> pull(Symbol) |> unique()
  length(sperm_list)

  df <- df |>
    mutate(Sperm_DNA_condensation_gene = gene %in% sperm_list) 

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = sperm_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Sperm_DNA_condensation_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)
  
  return(df)
}



# meiosis
tag_meiosis_genes <- function(df) {
  meiosis_table <- read_excel('reference_gene_sets/GO_gene_lists/GO_G1_M_transition_meiosis.xlsx')
  meiosis_list <- meiosis_table |> pull(Symbol) |> unique()
  length(meiosis_list)

  df <- df |>
    mutate(Meiosis_gene = gene %in% meiosis_list) |> unique()

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = meiosis_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
  mutate( Meiosis_gene = TRUE,
          Not_expressed  = TRUE,
          gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}


tag_oocyte_genes <- function(df) {
  # add the oocyte maturation signature for basal level
  oocyte_table <- read_excel('reference_gene_sets/GO_gene_lists/GO_oocyte_maturation.xlsx')
  oocyte_list <- oocyte_table |> pull(Symbol) |> unique()
  length(oocyte_list)

  df <- df |>
    mutate(Oocyte_maturation_gene = gene %in% oocyte_list) 

  # tag genes not in df
  all_clusters <- df |> distinct(cluster_id, cell_type, tissue, age)
  all_combinations <- crossing(all_clusters, gene = oocyte_list)
  missing_genes <- anti_join(all_combinations, df, 
                           by = c("cluster_id", "cell_type", "gene", "tissue", "age"))
  missing_genes <- missing_genes |>
    mutate(Oocyte_maturation_gene = TRUE,
            Not_expressed  = TRUE,
            gmean = 0)
  df <- bind_rows(df, missing_genes)

  return(df)
}


# tag not expressed genes (form either FACS or droplet all detected genes)

tag_not_expresssed_genes <- function(df, data_source) {

  data_source = "facs"
   if (data_source == "facs") {
        print("Data source: FACS.")
        all_genes_list <- unique(readRDS("reference_gene_sets/all_genes_facs.rds"))
    } else if (data_source == "droplet") {
      print("Data source: Droplet.")
      all_genes_list <- unique(readRDS("reference_gene_sets/all_genes_droplet.rds"))
    } else {
      stop("Issue in data source or cell type renaming.")
    } 
  
 
  #create df with every gene in gene list for every cell type and cluster
  selection <- df |>
    select(cluster_id, cell_type) |>
    distinct()
 
  all_genes_df <- crossing(selection, gene = all_genes_list)
  df$expressed <- TRUE

  all_genes_df_joined <- all_genes_df |>
   left_join(df, by = c("gene", "cluster_id", "cell_type")) |>
    mutate(expressed = replace_na(expressed, FALSE),   # false if not in df
          Not_expressed = !expressed,
          gmean = if_else(Not_expressed, 0, gmean))|>
    select(-expressed) 
  
  return(all_genes_df_joined)
}


#-------
# renaming etc. 

# get rid of NA in the gene sets
set_NA_false <- function(df) {
  df <- df |>
    mutate(across(ends_with("_gene"), ~replace_na(.x, FALSE)))
  return(df)
}

set_gene_variability <- function(df) {
  # make not expressed genes 0
  df <- df |> 
    mutate(Not_expressed = gmean == 0) |>
    mutate(gmean = if_else(Not_expressed, 0, gmean))

  # put genes into gene set categories
  df_main_with_markers <- df |>
    mutate(
      gene_variability = case_when(
        hvg == TRUE ~ "HVG",  
        lvg == TRUE ~ "LVG", 
        Not_expressed == TRUE ~ "Not expressed", 
       # Not_expressed_immune_response_gene == TRUE ~ "Not expressed", 
        TRUE ~ "Intermediate"))
  return(df_main_with_markers)
}


rename_cluster_id <- function(df) {
  df <- df |>
    mutate(cluster_name = str_replace(cluster_id,
      "^(.+?)_3m_(female|male)_(\\d+)$",
      "\\1 (\\2, cluster \\3)"
    ) |> str_replace_all("_", " "))
  return(df)
}


rename_gene_sets <- function(df) {
  df <- df |> rename_with(~ gsub("_", " ", .x), ends_with("_gene"))
  existing_gene_sets <- names(df)[endsWith(names(df), " gene")]

  return(list(df = df, existing_gene_sets = existing_gene_sets))
}



tag_no_gene_set_genes <- function(df, existing_gene_sets, data_source) {
  df <- df |>
     mutate(`Other gene` = if_all(all_of(existing_gene_sets), ~ .x == 0))

  saveRDS(df, paste0(data_source, '/data/gene_set_df.rds'))
  print("Saved df to data.")
  return(df) 
}

