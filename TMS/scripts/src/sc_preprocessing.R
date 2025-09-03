
library(Seurat)

# facs preprocessing is still missing


assemble_TMS_seurat_droplet <- function(){
    path <- "droplet_raw_data/cell_level/processed_10_22/3m"
    files <- list.files(path, pattern = "\\processed.rds$", full.names = TRUE)
    object_list <- lapply(files, readRDS)
    sample_ids <- gsub("\\_processed.rds$", "", basename(files))

    # rename cell ids
    for (names in seq_along(object_list)) {
    seurat_list[[names]] <- RenameCells(object_list[[names]], add.cell.id = sample_ids[names])
    }
    # merge objects
    main_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)



    # add needed metadata cols
    main_seurat$age <- as.integer(sub("m_.*", "", main_seurat$age_sex))
    main_seurat$cluster_id <- (paste0(main_seurat$tissue , "_" , main_seurat$age_sex , "_" , main_seurat$seurat_clusters))
    head(main_seurat@meta.data)
    unique(main_seurat$cluster_id)


    # import the df with metadata info
    df_droplet <- read.csv("droplet_analysis/data/combined_data.csv") 

    df_droplet_filtered <- df_droplet |>
    mutate(age = str_extract(age, "\\d+m")) |>  
    mutate(age = as.numeric(str_remove(age, "m"))) |> 
    filter(age == 3) |>
    rename(cell_type = manual_final) |> 
    filter(!(is.na(cell_type)))
    colnames(df_droplet_filtered)
    unique(df_droplet_filtered$cell_type)


    # map cell type
    cluster_map <- df_droplet_filtered[, c("cluster_id", "cell_type")]
    cluster_map <- unique(cluster_map)

    # assign to seurat object
    main_seurat$cell_type <- cluster_map$cell_type[match(main_seurat$cluster_id, cluster_map$cluster_id)]
    main_seurat$cluster_id <- paste0(main_seurat$tissue, " (", main_seurat$sex, ", cluster ", main_seurat$seurat_clusters, ")")
    main_seurat$ident <- as.factor(paste0(main_seurat$cell_type , " " , main_seurat$cluster_id))

    # set Idents
    Idents(main_seurat) <- main_seurat$ident
    Idents(main_seurat)
    unique(main_seurat@meta.data$cell_type)
    table(main_seurat@meta.data$cell_type)

    # save
    #saveRDS(main_seurat, "droplet/data/fused_seurat_droplet_3m.rds")



    # subset for just immune cell types

    # load cell types
    cell_type_names <-  read_xlsx('droplet_analysis/data/cell_type_list.xlsx')
    cell_type_name_list <- cell_type_names |> pull(unique.combined_data.manual_final.) |> unique() 
    immune_cell_type_name_list <- c( "leukocyte_macrophage", "leukocyte_progenitor", "lymphocyte_B", "lymphocyte_T", 
    "leukocyte_monocyte" , "myeloid" , "lymphocyte" , "lymphocyte_plasma" , 
    "lymphocyte_NK" , "leukocyte_dc" ,  "lymphocyte_NKT" , "leukocyte_neutrophil", 
    "leukocyte_basophil"  ,  "leukocyte_granulocyte" , "granulocyte_progenitor", 
    "leukocyte" ,  "lymphocyte_thymocyte"   )
    immune_cell_type_name_list

    # filter seurat for immune cells
    immune_seurat <- subset(x = main_seurat, subset = cell_type %in% immune_cell_type_name_list)
    unique(immune_seurat@meta.data$cell_type)
    #saveRDS(immune_seurat, "droplet/data/immune_seurat_droplet_3m.rds")

}