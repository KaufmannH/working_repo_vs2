
# This is for getting a feeling for the data quality. 


mac_marker_expression <- function(data_source, df){
  marker_genes <- c('Fcgr1', 'Mertk',  'Spi1', 'Maf', 'Mafb', 'Cebpa', 'Cd68', 'Cd14', 'Adgre1', 'Itgam', 'Siglec1')
  filtered_df <- df |>
    filter(gene %in% marker_genes) 

  plot <- ggplot(filtered_df, aes(x = cluster_name, y = gene)) +
    geom_point(aes(colour = gmean), size = 6) +
    scale_colour_gradient(low = "lightgreen", high = "purple") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12)) +
    labs(
      x = "Cluster",
      y = "Macrophage marker gene",
      colour = "gmean")

 if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_13/makro_markers.png", plot, width = 20, height = 8)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_13/makro_markers.png", plot, width = 20, height = 8)
    } else {
      stop("Issue when saving.")
    }   
}



genes_per_cluster <- function(data_source, df){
  # check how many clusters the genes are expressed in
  num_cells_express_gene_df <- df |>
    group_by(gene) |>
    mutate(gene_count = n()) |>
    summarise(gene, gene_count, cluster_id)
  print(num_cells_express_gene_df, n= 100)

  qc <- ggplot(num_cells_express_gene_df, aes(x = gene_count)) +
    geom_histogram(binwidth = 1, fill = "grey", colour = "white") +
    labs(x = "Number of clusters a gene is expressed in",
        y = "Number of genes") +
    theme_classic()

    print("Saving plot: num_clusters_gene_is_expressed.png")

   if (data_source == "facs") {
    ggsave('facs/plots/3_m/num_clusters_gene_is_expressed.png', plot = qc, width = 8, height = 11)
    
  } else if (data_source == "droplet") {
    ggsave('droplet/plots/3_m/num_clusters_gene_is_expressed.png', plot = qc, width = 8, height = 11)
    
  } else {
    stop("Issue in data source.")
  }
  
  return(df)
}



t <- function(df){
    

immune_droplet <- readRDS("droplet/data/immune_seurat_droplet_3m.rds")
# number of mice in data

metadata_raw <- read.table('droplet/droplet_raw_data/combined_metadata_for_cell_numbers.txt', header = TRUE, sep = "\t")
metadata <-  metadata_raw |> mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters)) 
head(metadata)
# mit mouse id

met_raw <- read_xlsx("/home/hkaufm49/working_repo/TMS/droplet/droplet_raw_data/manual_annotation.xlsx")
# mit manual final
head(met_raw)
joined <- metadata   |>
    left_join(met_raw, by = c("cluster_id")) 
unique(joined$manual_final)


immune_cell_type_list <- c(
"leukocyte_progenitor", "lymphocyte_B",  "myeloid", "lymphocyte", "leukocyte", "leukocyte_mast", "lymphocyte_NK", "leukocyte_dc", 
"lymphocyte_T", "leukocyte_monocyte", "lymphocyte_plasma", "lymphocyte_NKT", "granulocyte_progenitor", 
"leukocyte_neutrophil", "leukocyte_granulocyte", "leukocyte_macrophage", "leukocyte_basophil" )

metadat_immune <- joined |> 
  filter(age == "3m" ) |>
  filter(manual_final %in% immune_cell_type_list)
head(metadat_immune)

unique((metadat_immune$mouse.id))

#  get avg number of genes per cell
immune_droplet$gene_count <- Matrix::colSums(immune_droplet@assays$originalexp@counts > 0)
mean_gene_count <- mean(immune_droplet$gene_count)
mean(immune_droplet$gene_count)
}

