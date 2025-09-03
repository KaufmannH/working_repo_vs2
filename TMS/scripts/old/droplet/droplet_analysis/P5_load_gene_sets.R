# prepare dataframe for analysis: load and tag genes for all conditions to be tested downstream. 

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(ggvenn)
library(purrr)
library(tibble)

# innate immune response genes ()
# housekeeping genes ()
# ...



# add and filter main df from TMS
df_main_all <- read.csv("droplet_analysis/data/combined_data.csv") 
colnames(df_main_all)
unique(df_main_filtered$cell_type)

num_cells <- df_main_all %>%
  group_by(cluster_id, num_cells_per_cluster) %>%
  select (cluster_id, num_cells_per_cluster) %>%
  unique()
tail(num_cells)
sum(num_cells$num_cells_per_cluster, na.rm = TRUE) # one cluster without annotation (non expressed innate genes)

# filter for age and cell types
df_main_filtered <- df_main_all  %>% 
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(age = as.numeric(str_remove(age, "m"))) %>% 
  filter(age == 3) %>%
  #filter(tissue == "Lung") %>% # until here 3460 cells
  #filter(manual_final %in% c("leukocyte_dc", "leukocyte_macrophage", "leukocyte_monocyte")) %>%
  filter(manual_final %in% c("leukocyte_macrophage")) %>%
  mutate(cell_type = manual_final) %>% 
  mutate(cell_type = case_when(
      str_detect(cell_type, "leukocyte_dc") ~ "DC",
      str_detect(cell_type, "leukocyte_macrophage") ~ "Macrophage",
      str_detect(cell_type, "leukocyte_monocyte") ~ "Monocyte",
      TRUE ~ as.character(cell_type))) %>%
  select(cluster_id, cell_type, gene, res_var, gmean, hvg, lvg, num_cells_per_cluster)
head(df_main_filtered)
sum(df_main_filtered$inner_gene)
dim(df_main_filtered)
unique(df_main_filtered$cell_type)

# all genes list
all_genes <- unique(df_main_filtered$gene)


num_cells <- df_main_filtered %>%
  group_by(cluster_id, num_cells_per_cluster) %>%
  select (cluster_id, num_cells_per_cluster) %>%
  unique()
sum(num_cells$num_cells_per_cluster, na.rm = TRUE)


# load the mouse/human conversions
m2h <- read.delim("droplet_analysis/data/20200307_ensembl/mouse.txt",  header = FALSE, sep = "\t")
colnames(m2h) <- c("mouse", "human")
h2m <- read.delim("droplet_analysis/data/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
colnames(h2m) <- c("human", "mouse")
head(m2h)


# add and tag innate response genes 

#load the innate immunity genes of mice database # 646 genes
mouse_innate_genes_df_raw <- read.csv('droplet_analysis/data/mouse_innate_genes.csv')
mouse_innate_genes <- mouse_innate_genes_df_raw %>% 
    select(Gene.Symbol) %>%
    rename(gene = Gene.Symbol) %>%
    unique()
dim(mouse_innate_genes)
# add to main df
df_response_genes <- df_main_filtered %>%
  full_join(mouse_innate_genes, by = "gene") %>%
  mutate(response_not_expressed = is.na(gmean)) %>%
  mutate(response_expressed = !response_not_expressed) %>%
  mutate(across(c( hvg, lvg), ~replace_na(., FALSE)))
tail(df_response_genes)
sum(df_response_genes$response_not_expressed)
sum(!df_response_genes$response_not_expressed) # 104 are not expressed
dim(df_response_genes)


# add and tag housekeeping genes
# general house keeping genes
housekeeping_df <- read.csv('droplet_analysis/data/Housekeeping_TranscriptsMouse.csv', sep = ";") 
housekeeping_list <- housekeeping_df %>%
    pull(Genes)
housekeeping_list

# housekeeping from Lin et al.
## https://academic.oup.com/gigascience/article/8/9/giz106/5570567

housekeeping_lin_list <- read_lines("droplet_analysis/data/scHK_human.symbols.txt")

mouse_lin <- h2m |>
  filter(human %in% housekeeping_lin_list) |>
  pull(mouse)
unique(mouse_lin)

# add and tag housekeeping genes
df_main <- df_response_genes %>%
    # tag house keeping genes
    mutate(housekeeping_gene_lin = (str_detect(gene, paste(mouse_lin, collapse = "|")))) %>% 
    mutate(housekeeping_gene = (str_detect(gene, paste(housekeeping_list, collapse = "|")))) %>% 
    mutate(mouse_id = str_extract(cluster_id, ".*(?=_[^_]+$)")) 
head(df_main)




# add markers: get Marker Genes

markers_raw  <- read.delim("droplet_analysis/data/lineage_specific_markers.tsv", header = TRUE, sep = "\t")
available_cell_types_big = c('Alveolar macrophages', 'Kupffer cells', 'Macrophages', 'Langerhans cells', 'Microglia', 'Monocytes', 'Basophils', 'Dendritic cells', 'Eosinophils', 'NK cells', 'Neutrophils')
available_cell_types_small = c('Dendritic cells','Macrophages','Monocytes', 'Alveolar macrophages') 


markers_df <- markers_raw %>%
  filter(grepl("Mm", species)) %>%
  rename(gene = official.gene.symbol) %>%
  rename(cell_type = cell.type) %>%
  mutate(gene = str_to_lower(gene),           
         gene = str_replace(gene, "^(.)", str_to_upper)) %>%
  filter(cell_type %in% available_cell_types_small) %>%
  select(gene, cell_type) %>%
  mutate(cell_type = case_when(
    cell_type == "Macrophages" ~ "Macrophage",
    cell_type == "Alveolar macrophages" ~ "Macrophage",
    cell_type == "Monocytes" ~ "Monocyte",
    cell_type == "Dendritic cells" ~ "DC")) %>%
  unique() %>%
  #group_by(cell_type) %>%
  #summarise(count = n()) %>%
  mutate(is_marker = TRUE)
head(markers_df)
#write.csv(markers_df, "data/lineage_specific_markers.csv", row.names = FALSE)

# get the lists extra
macrophage_unique <- markers_df %>%
  filter(cell_type == 'Macrophage') %>%
  distinct(gene) %>%
  pull(gene)

monocyte_unique <- markers_df %>%
  filter(cell_type == 'Monocyte') %>%
  distinct(gene) %>%
  pull(gene)

dc_unique <- markers_df %>%
  filter(cell_type == 'DC') %>%
  distinct(gene) %>%
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
df_main_markers <- df_main %>%
  left_join(markers_df, by = c("gene", "cell_type")) %>%
  mutate(is_marker = ifelse(is.na(is_marker), FALSE, is_marker))
head(df_main_markers)

# check outcome
marker_summary <- df_main_markers %>% 
  group_by(cell_type, is_marker) %>%
  summarise(count = n(), .groups = 'drop')
marker_summary
 
#_____

# load LPS stimulated genes in monocytes
lps_df <- read_excel("droplet_analysis/data/Bhatt_2012_data.xlsx")
colnames(lps_df)[colnames(lps_df) %in% c("Probe...1", "Probe...2", "Promoter Class", "...4")] <- c("gene", "gene_id", "promoter_class", "group")
lps_df_early <- lps_df %>% filter(group <= 5)
lps_df_late <- lps_df %>% filter(group > 5)
lps_gene_list_early <- lps_df_early %>% pull(gene)
lps_gene_list_late <- lps_df_late %>% pull(gene)

#stopped here now too
df_lps <- df_main %>% 
 mutate(lps_stim_early = str_detect(gene, paste(lps_gene_list_early, collapse = "|"))) %>%
 mutate(lps_stim_late = str_detect(gene, paste(lps_gene_list_late, collapse = "|"))) %>%
   mutate(lps_stim = lps_stim_early | lps_stim_late) 
head(df_lps) 


# check outcome
lps_summary <- df_lps %>%  
  group_by(cell_type, lps_stim) %>%
  summarise(count = n(), .groups = 'drop')
lps_summary
#saveRDS(df_lps, "data/df_all_cells_lps_tag.rds") 


# load chemokine sinalling pathway
cytokine_df <- read.csv("droplet_analysis/data/InnateDB_cytokine_signalling.csv")
head (cytokine_df)
cytokine_gene_list <- cytokine_df %>%
  pull(name)
cytokine_gene_list

df_main_cytokine <- df_lps %>%
 mutate(cytokine_gene = str_detect(gene, paste(cytokine_gene_list, collapse = "|"))) #%>%
head(df_main_cytokine)



# load TLR pathway
tlr_df <- read.csv("droplet_analysis/data/InnateDB_genes_TLR.csv")
head (tlr_df)
tlr_gene_list <- tlr_df %>%
  pull(name)
tlr_gene_list

df_main_tlr <- df_main_cytokine %>%
 mutate(tlr_gene = str_detect(gene, paste(tlr_gene_list, collapse = "|"))) #%>%
sum(df_main_tlr$tlr_gene)


# load xue et al. response genes
xue_response_raw <- readRDS("droplet_analysis/data/xue_response_genes.rds")
xue_response <- xue_response_raw %>% select(Condition, Upregulated)

# see if the gene is in def and then pull the ortholog
xue_response_orthologs <- xue_response %>%
  mutate(Orthologs = map(Upregulated, ~ h2m$mouse[match(.x, h2m$human)]))
xue_response_orthologs

xue_gene_list <- xue_response_orthologs %>% 
  select(Condition, Orthologs) %>% 
  unnest(Orthologs) %>% 
  group_by(Condition) %>% 
  summarise(GeneList = list(unique(Orthologs))) %>% 
  deframe()  

df_main_xue <- df_main_tlr %>%
 mutate(ifng_gene = gene %in% xue_gene_list[["IFNg"]]) %>%
 mutate(il10_gene = gene %in% xue_gene_list[["IL10"]]) %>%
 mutate(il4_gene = gene %in% xue_gene_list[["IL4"]]) %>%
mutate(tnf_gene = gene %in% xue_gene_list[["TNF"]]) 



# add negative control gene sets

# add the spern DNA condensation signature for basal level
sperm_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_sperm_dna_condensation.xlsx')
sperm_list <- sperm_table |> pull(Symbol)
unique(sperm_list)

df_neg_contrl <- df_main_xue %>%
  mutate(sperm_gene = gene %in% sperm_list) 


# add the meiosis signature for basal level
meiosis_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_G1_M_transition_meiosis.xlsx')
meiosis_list <- meiosis_table |> pull(Symbol)
unique(meiosis_list)

df_neg_contrl <- df_neg_contrl %>%
  mutate(meiosis_gene = gene %in% meiosis_list) 


# add the oocyte maturation signature for basal level
oocyte_table <- read_excel('droplet_analysis/data/GO_gene_lists/GO_oocyte_maturation.xlsx')
oocyte_list <- oocyte_table |> pull(Symbol)

df_neg_contrl <- df_neg_contrl %>%
  mutate(oocyte_gene = gene %in% oocyte_list) 

# check how many clusters the genes are expressed in
num_cells_express_gene_df <- df_neg_contrl %>%
  group_by(gene) %>%
  mutate(gene_count = n()) %>%
  #arrange(gene) %>%
  select(gene, gene_count, cluster_id)
  #pull(gene_count)
print(num_cells_express_gene_df)



qc <- ggplot(num_cells_express_gene_df, aes(x = (gene_count))) +
  geom_histogram(binwidth = 1, fill = "grey", colour = "white") +
  labs(x = "Number of clusters a gene is expressed in",
       y = "Number of genes") +
  theme_classic()
ggsave('droplet_analysis/plots/3_m/num_clusters_gene_is_expressed.png', plot = qc, width = 8, height = 11)


# add autimmune pleiotropic loci
autimmune_df <- read_excel('droplet_analysis/data/Marquez_pleiotropic_genes.xlsx')

autoimmune_list <- autimmune_df |> pull(Gene)
autoimmune_list

# conversion from human to mouse
mouse_auto <- h2m |>
  filter(human %in% autoimmune_list) |>
  pull(mouse)

df_neg_contrl <- df_neg_contrl %>%
  mutate(autoimmune_gene = gene %in% mouse_auto) 


# put genes into gene set categories
df_main_with_markers <- df_neg_contrl %>%
  mutate(
    gene_set = case_when(
      hvg ~ "HVG",  
      lvg ~ "LVG", 
      response_not_expressed ~ "Response gene not expressed", 
      TRUE ~ "Other")) %>%
  mutate(Other = !hvg & !lvg & !response_not_expressed &  !housekeeping_gene_lin & !housekeeping_gene) %>%
  select(gmean, gene, cell_type, gene_set, res_var, cluster_id, num_cells_per_cluster, response_expressed, housekeeping_gene_lin, housekeeping_gene, response_not_expressed, lps_stim, lps_stim_early,lps_stim_late, cytokine_gene, tlr_gene, ifng_gene, il10_gene, il4_gene, tnf_gene, sperm_gene, meiosis_gene, autoimmune_gene, Other)
head(df_main_with_markers, n = 100)

table(df_main_with_markers$gene, df_main_with_markers$cluster_id) %>% as.data.frame() %>% as_tibble()





high_threshold <- quantile(df_main_with_markers$gmean, 0.75, na.rm = TRUE)
low_threshold <- quantile(df_main_with_markers$gmean, 0.25, na.rm = TRUE)

df_main_with_markers_expr_level <- df_main_with_markers %>%
  mutate(
    expression_level = case_when(
      is.na(gmean) ~ "NA",  
      gmean >= high_threshold ~ "Highly Expressed",
      gmean > low_threshold ~ "Medium Expressed",
      TRUE ~ "Lowly Expressed"
    )
  )


check_cell_types <- (df_main_with_markers_expr_level) %>%
  group_by(cell_type) %>%
  summarise(count = n()) 
check_cell_types # there are the same number of DC and macrophage genes bc the same bio repl and tissue have the same # genes 

#saveRDS(df_main_with_markers_expr_level, "droplet_analysis/data/df_with_markers_2.rds") 
