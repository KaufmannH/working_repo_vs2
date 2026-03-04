# What is the expression level of immune response genes in comparison to lineage specific genes?
# only in 3m and lung

# here PanSci
library(readxl)
library(patchwork)
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)
library(viridis)


# colors
palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")



## 1. Data preparation

# already lung and 3 months data only - no filtering required
df_main_all <- read_csv("./data/pansci_lung_3m_small.csv")
(df_main_all$gmean)
# get only certain immune cells

# get the innate immune response gene list
# get the house keepign gene list
# get the lineage specific gene list


#load the innate immunity genes of mice database
mouse_innate_genes_df <- read.csv('data/mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df %>% 
    pull(Gene.Symbol)

# load lineage specific genes
file_paths <- list.files(path = "data/lineage_spec_markers", pattern = "\\.tsv$", full.names = TRUE)

data_list <- lapply(file_paths, function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  cell_type <- sub("\\.tsv$", "", basename(file))  
  df$cell_type <- cell_type 
  return(df)})
markers_df <- bind_rows(data_list)
head(markers_df$species) 

markers_df <- markers_df %>%
    select(official.gene.symbol, cell_type, species, organ, specificity_mouse, sensitivity_mouse, gene.type) %>%
    filter(grepl("Mm", species)) %>%
    rename(gene = official.gene.symbol) %>%
    rename(cell_type_for_marker = cell_type) 
head(markers_df)

# house keeping
housekeeping_df <- read.csv('data/Housekeeping_TranscriptsMouse.csv', sep = ";") 
housekeeping_list <- housekeeping_df %>%
    pull(Genes)
housekeeping_list


# load LPS stimulated genes in monocytes
lps_df <- read_excel("./data/Bhatt_2012_data.xlsx")
colnames(lps_df)[colnames(lps_df) %in% c("Probe...1", "Probe...2", "Promoter Class", "...4")] <- c("gene", "gene_id", "promoter_class", "group")
lps_gene_list <- lps_df %>%
  pull(gene)
lps_gene_list

colnames(df_main_all)

# add info on genes
df_main <- df_main_all  %>% 
    # tag innate immune genes
     mutate(inner_gene = (str_detect(gene, paste(mouse_innate_genes, collapse = "|")))) %>% 
    # tag lps response genes
    mutate(lps_stim = str_detect(gene, paste(lps_gene_list, collapse = "|"))) %>%
    # tag house keeping genes
    mutate(housekeeping_gene = (str_detect(gene, paste(housekeeping_list, collapse = "|")))) %>% 
    relocate(cell_type, gmean, inner_gene, housekeeping_gene, .after = 4)# %>%
df_main$lps_stim




# vector for renaming the cell types (to fit the cell type specific marker info)
cell_type_name_change = c( 
    'Myeloid_cells_Alveolar macrophages' =  'leukocyte_macrophage', 
     'Myeloid_cells_Monocytes'= 'leukocyte_monocyte',
      'Myeloid_cells_Neutrophils'= 'leukocyte_neutrophil',
       'Myeloid_cells_Basophils'= 'leukocyte_basophil',
        'Myeloid_cells_Interstitial macrophages'= 'leukocyte_macrophage',
         'Myeloidcells_Dendritic cells' = 'leukocyte_dc'
    )

# Rename values in the 'cell_type' column
df_main <- df_main  %>%
  mutate(cell_type_for_marker = recode(cell_type, !!!cell_type_name_change)) %>%
  #filter(cell_type_for_marker %in% c(leukocyte_basophil, leukocyte_dc, leukocyte_macrophage, leukocyte_monocyte, leukocyte_neutrophil))#%>%
  filter(!is.na(cell_type_for_marker))
head(df_main)

cell_type_summary <- df_main  %>%
  group_by(cell_type_for_marker) %>%
  summarise(n())
cell_type_summary


# join dfs to tag markers
df_main_with_markers <- df_main %>%
  # add lineage specific markers
  left_join(markers_df %>% 
              mutate(is_marker = TRUE), 
            by = c("cell_type_for_marker"))   %>%
  mutate(is_marker = ifelse(is.na(is_marker), FALSE, is_marker))  %>%
  mutate(gene_set = case_when(
  inner_gene & is_marker ~ "Innate & Lineage Specific",  # housekeeping not yet in there
  inner_gene ~ "Innate Response Gene",
  is_marker ~ "Lineage Specific",
  housekeeping_gene ~ "Housekeeping Gene", 
  TRUE ~ "Other")) %>%
  group_by(cell_type, gene_set) %>%
  mutate(mean_expression = mean(gmean)) %>%
  ungroup() %>%
  select(mean_expression, gmean, cell_type, gene_set)
print(df_main_with_markers, n = 20)


# check how many are both innate and lineage specific
innate_lineage_count <- df_main_with_markers %>%
  filter(gene_set == "Innate & Lineage Specific") %>%
  summarise(count = n())
innate_lineage_count

# reporting
heatmap_plot <- ggplot(df_main_with_markers, aes(x = cell_type, y = gene_set, fill = mean_expression)) +
  geom_tile(color = "black") +
  scale_fill_viridis(option = "C", name = "Mean gene expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 
ggsave('./plots/heatmap_lineage_expression_3m_lung_pansci.png', plot = heatmap_plot)

# to make it wider, add to ggsave
#width = 8, height = 4, dpi = 150      

