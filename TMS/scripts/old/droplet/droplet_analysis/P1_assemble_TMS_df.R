
# merge TMS data to one csv

library(data.table)
library(dplyr)
library(readxl)
library(writexl)


# load data
gene_expression_data <- read.delim("./droplet/full_results_10_22.tsv", header = TRUE, sep = "\t")
annotation_info <- read_excel("./droplet/manual_annotation.xlsx")
metadata_cell_numbers_raw <- read.table('./droplet/combined_metadata_for_cell_numbers.txt', header = TRUE, sep = "\t") # there are cluster_ids that are not in the annotaiton info


# prep dfs
gene_expression_data <- gene_expression_data %>%
  mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) 
head(gene_expression_data)
# every biological replicate and tissue have the same number of rows (genes)

metadata_cell_numbers <- metadata_cell_numbers_raw %>%
  mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters)) %>%
  group_by(cluster_id) %>%
  summarise(num_cells_per_cluster = n(), .groups = 'drop') 
(metadata_cell_numbers)
sum(metadata_cell_numbers$num_cells_per_cluster)

# join dfs
numbers_joined <- gene_expression_data  %>%
  left_join(metadata_cell_numbers, by = "cluster_id") 
dim(numbers_joined)

all_joined <- numbers_joined  %>%
  left_join(annotation_info, by = c("cluster_id", "tissue")) 
dim(all_joined)

# in previous version no cell numbers and merge used:
#combined_data <- merge(gene_expression_data, meta_data_joined, by = c("cluster_id", "tissue"))

# tag LVGs
df_finished <- all_joined %>%
  mutate(lvg = if_else(res_var < 1, TRUE, FALSE))

# write combined df to a csv
#write.csv(df_finished, './droplet_analysis/data/combined_data.csv')

# get cell type names, write to excel for manual selection of innate immune cells
cell_type_list <- data.frame(unique(df_finished$manual_final))
write_xlsx(cell_type_list, 'droplet_analysis/data/cell_type_list.xlsx')

