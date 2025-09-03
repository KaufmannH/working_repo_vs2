
# merge TMS data to one csv

library(data.table)
library(dplyr)
library(readxl)
library(writexl)
library(Seurat)



# load data
obj <- readRDS("/home/hkaufm49/working_repo/TMS/droplet/cell_level/processed_10_22/Lung_3m_female_processed.rds")


# prep dfs
annotation_info <- annotation_info_raw %>%
  mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters))

gene_expression_data <- gene_expression_data %>%
  mutate(cluster_id = paste0(tissue , "_", age, "_", cluster))


metadata_cell_numbers <- metadata_cell_numbers_raw %>%
  mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) %>%
  select(- age, -cluster, -tissue, -tissue_cluster)
head(metadata_cell_numbers)

# join dfs
numbers_joined <- gene_expression_data  %>%
  left_join(metadata_cell_numbers, by = "cluster_id") %>%
  select(-cell_ontology_class, -tissue) 
dim(numbers_joined)

all_joined <- numbers_joined  %>%
  left_join(annotation_info, by = c("cluster_id")) 
dim(all_joined)

# tag LVGs
df_finished <- all_joined %>%
  mutate(lvg = if_else(res_var < 1, TRUE, FALSE))
colnames(df_finished)

test <- df_finished |>
  filter(cluster_id  == "Lung_3m_male_13" ) |>
  group_by(cluster_id, gene, num_cells) |>
  select(cluster_id, gene, gmean, num_cells)
print(test, n = 500)

# write combined df to a csv
#write.csv(df_finished, 'facs_analysis/data/combined_data.csv')

# get cell type names, write to excel for manual selection of innate immune cells
cell_type_list <- data.frame(unique(df_finished$manual_final))
write_xlsx(cell_type_list, 'facs_analysis/data/cell_type_list.xlsx')

