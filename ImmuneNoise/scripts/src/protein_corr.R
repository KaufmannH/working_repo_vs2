
library(readxl)
library(dplyr)
library(stringr)
library(matrixStats)
library(purrr)
library(tidyr) 
library(ggplot2)
library(GGally)
library(patchwork)
library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggridges)
# 


# 1. QC of bone marrow mice proteomics from natoli (confidential)

data <- read_xlsx('ImmuneNoise/proteomics_marrow_natoli.xlsx', sheet = 'LPS_Prote', col_names = TRUE )

# select only the steady state data
data_filtered <- data |>
    mutate(gene = str_split_i(ID...1, pattern = "_", i = 1)) |>
    mutate(id = str_split_i(ID...1, pattern = "_", i = 2)) |>
    dplyr::select(c(gene, id, UT_R1, UT_R2, UT_R3, UT_R4, UT_R5))
data_filtered



# Quantile normalization (Bolstad et al. 2003) according to https://ab604.github.io/docs/bspr_workshop_2018/transform.html
quantile_normalisation <- function(df){
  
  # Find rank of values in each column
  df_rank <- map_df(df,rank,ties.method="average")
  # Sort observations in each column from lowest to highest 
  df_sorted <- map_df(df,sort)
  # Find row mean on sorted columns
  df_mean <- rowMeans(df_sorted)
   
  # Function for substiting mean values according to rank 
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  # Replace value in each column with mean according to rank 
  df_final <- map_df(df_rank,index_to_mean, my_mean= df_mean)
  
  return(df_final)
}

data_normalized <- data_filtered |> 
    dplyr::select(-c(gene:id)) |> 
    quantile_normalisation() |> 
    bind_cols(data_filtered[,1:2]) |>
    dplyr::select(c(gene, id), everything())
data_normalized


data_log <- data_normalized |>
   mutate(across(c(UT_R1, UT_R2, UT_R3, UT_R4, UT_R5), ~ log1p(.x),  .names = "{.col}_log")) 
data_log

data_qc <- data_log |>              
    mutate(mean = rowMeans(across(ends_with("_log")))) |>
    mutate(var = rowVars(as.matrix(across(ends_with("_log"))), na.rm = TRUE)) |>
    mutate(sd  = sqrt(var)) |>
    dplyr::select(c(gene, id, mean, var, sd), everything())
data_qc



# save df
write_tsv(data_qc, 'ImmuneNoise/proteomics_marrow_natoli_prepped.tsv')



# distribution of expression between batches

df_long <- data_qc |>
  pivot_longer(cols = ends_with("_log"),
               names_to = "batch",
               values_to = "log_intensity")
df_long[, c('batch', 'log_intensity')]

# after normalization
plot_after <- ggplot(df_long, aes(x = log_intensity, fill = batch, color = batch)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  labs(
    x = "log intensity",
    y = "Count",
    title = "Distribution of expression across batches after normalization") +
  theme_classic(base_size = 15)


# before normalization
data_log_before <- data_filtered |>
 mutate(across(c(UT_R1, UT_R2, UT_R3, UT_R4, UT_R5), ~ log1p(.x),  .names = "{.col}_log")) 
data_log_before


df_long_before <- data_log_before |>
  pivot_longer(cols = ends_with("_log"),
               names_to = "batch",
               values_to = "log_intensity")
df_long_before[, c('batch', 'log_intensity')]

# after normalization
plot_before <- ggplot(df_long_before, aes(x = log_intensity, fill = batch, color = batch)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  labs(
    x = "log intensity",
    y = "Count",
    title = "Distribution of expression across batches before normalization") +
  theme_classic(base_size = 15)


plot <- plot_before | plot_after

ggsave('comparison/plots/proteomics/protein_marrow_dist.png', plot, width = 15, height = 7)






# 2. correlation of turnover and variability of proteomics

# EB cells: early differentiation
turnover <- read_xlsx('ImmuneNoise/protein_turnover_sabatier.xlsx', col_names = TRUE, skip = 1)
head(turnover)


# proteomics: 

# annotate the protein that belongs to the gene name and entrezgene id
mouse_ids <- as.character(data_qc$id)
data_annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = mouse_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENTREZID")

head(data_annot)



# convert mouse to human

mm <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

mouse_syms <- unique(na.omit(data_annot$SYMBOL))

# 1) mouse core IDs
m1 <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
  filters    = "mgi_symbol",
  values     = mouse_syms,
  mart       = mm
)

# 2) human homolog info
m2 <- getBM(
  attributes = c("ensembl_gene_id",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_orthology_type"),
  filters    = "mgi_symbol",
  values     = mouse_syms,
  mart       = mm
)

orth <- merge(m1, m2, by = "ensembl_gene_id", all.x = TRUE)
names(orth) <- c("mouse_ensembl","mouse_symbol","mouse_entrez",
                 "human_symbol","human_ensembl","orthology_type")

# prefer one-to-one
orth <- orth[order(orth$mouse_symbol, orth$orthology_type != "ortholog_one2one"), ]
orth_best <- orth[!duplicated(orth$mouse_symbol), ]

orth_selected <- orth_best |>
  dplyr::select(mouse_symbol, human_symbol)

# fuse
merged_df <- turnover |>
  dplyr::rename(gene_names = query_gene_names) |> 
  left_join(orth_selected, by = c("gene_names" = 'human_symbol')) |>
  left_join(data_qc, by = c('mouse_symbol' = "gene") )
merged_df

ridge_df <- merged_df |>
   mutate(mean_protein = factor(ntile(mean, 5), labels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"))) |>
   mutate(var_protein = factor(ntile(sd, 5), labels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")))




  
  col_vec_mean <- c("Group 1" = "#1E6057", "Group 2" = "#2D8F82", "Group 3" = "#3CBFAE", "Group 4" = "#89DACF", "Group 5" = "#B0E6DF")

plot <- ggplot(ridge_df, aes(x = turnover_value, y = as.factor(mean_protein), fill = as.factor(mean_protein))) +
  geom_density_ridges(colour = "white", size = 0.7) +
  labs(
    x = "Mean protein turnover (log2 (L/H))",
    y = "Intensity group", 
    fill = "Expression groups") +
     geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = col_vec_mean) +
    theme_classic() +
   theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16))  


ggsave('comparison/plots/proteomics/protein_marrow_ridge_mean.png', plot, width = 12, height = 8)



col_vec_var <- c("Group 1" = "#1E6057", "Group 2" = "#2D8F82", "Group 3" = "#3CBFAE", "Group 4" = "#89DACF", "Group 5" = "#B0E6DF")

plot <- ggplot(ridge_df, aes(x = turnover_value, y = as.factor(var_protein), fill = as.factor(var_protein))) +
  geom_density_ridges(colour = "white", size = 0.7) +
  labs(
    x = "Mean protein turnover (log2 (L/H))",
    y = "Intensity group", 
    fill = "Expression groups") +
     geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = col_vec_var) +
    theme_classic() +
   theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16))  
    
    ggsave('comparison/plots/proteomics/protein_marrow_ridge_var.png', plot, width = 12, height = 8)
