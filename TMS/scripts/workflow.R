# here funcitons are loaded and executed
# environments: tms_env (local), single_cell_r_env on VM, .libPaths(c("/home/hkaufm49/anaconda3/envs/single_cell_r_env", .libPaths()))

# load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(ggridges)
library(RColorBrewer)
library(writexl)
library(ggh4x)
library(clusterCrit)
library(dynamicTreeCut)
library(hexbin)
library(viridis)
library(future)
#library(progressr)
#library(beepr)


# load functions
source("scripts/src/assemble_TMS_df.R")
source("scripts/src/load_gene_sets.R")
source("scripts/src/stratification.R")
source("scripts/src/gene_distribution.R")
source("scripts/src/qc.R")
source("scripts/src/technology_comparison.R")
 
#----------


# 1. FACS

#basic_df_facs <- assemble_TMS_df_facs(write = TRUE)

# takes saved dfs and adds gene sets etc. 
#df_main_filtered <- load_filter_og_df(data_source = "facs")


# tag genes in which gene sets they are a part of
dict_df <- tag_imm_dict_genes(data_source = "facs")

df_hk <- tag_hk_genes(data_source = "facs")
df_hk_lin <- tag_hk_lin_genes(df_hk)
#df_response_genes <- tag_innate_response_genes(df_hk_lin) # too many
df_lps <- tag_lps_genes(df_hk_lin)
#df_marker <- tag_mac_marker_genes(df_lps) # too few
df_chemokine <- tag_chemokine_genes(df_lps)
df_tlr <- tag_tlr_genes(df_chemokine)
df_xue <- tag_xue_genes(df_tlr)
df_autoimmune <- tag_autoimmune_genes(df_xue)
df_sperm <- tag_sperm_genes(df_autoimmune)
df_meiosis <- tag_meiosis_genes(df_sperm)
df_oocyte <- tag_oocyte_genes(df_meiosis)
#df_all <- tag_not_expresssed_genes(df_oocyte, data_source = "facs")

# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(dict_df)
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "facs")

# stratification
strat_df <- stratify_df(data_source = "facs", cell_type_selection = "Macrophage") # also take df as arg
strat_df <- readRDS("facs/data/strat_df.rds")


colnames(strat_df)
# expression bins
e <- plot_expression_bins(strat_df, "facs")
#l <- plot_expression_bins_lps(strat_df, "facs")
n <- plot_expression_bin_numbers(strat_df, "facs")

# categories
h <- plot_category_numbers(strat_df, "facs")
g <- plot_gene_set_numbers(strat_df, "facs")

d <- plot_gene_set_proportions(strat_df, "facs")
j <- gene_set_proportions_compact(strat_df, "facs")

# qc
mac_marker_expression(df = strat_df, data_source = "facs")

# gene distribution
selection_gene_distr_df <- gene_distribution(strat_df,'facs') 
df_dummy <- dummy_entropy("facs")
all_genes_rao_prep <- all_genes_entropy_prep(strat_df, 'facs')
subsampled_df <- subsample_genes(all_genes_rao_prep) # faster entropy calc
df_entropy <- quadratic_entropy(all_genes_rao_prep, 'facs', write = TRUE) # or dummy data or selected genes data or subsampled data 
plot_raos_gene_selection(df_entropy, "facs")
plot_raos_all_genes(df_entropy, "facs")


master_entropy_df <- merge_entropy_results("facs")
variability_df <- calc_variability_direction(master_entropy_df, "facs")
#scatter_entropy_hvg(variability_df, "facs")
heatmap_entropy_bins(master_entropy_df, "facs")

# turnover
#add the mean data

# combine technologies
combined_df <- combine_technologies()
plot_scatter_rao()
plot_scatter_var_dir()
hirarchical_clustering(combined_df)


tail(master_entropy_df)
# test if hvg is there
t <- strat_df |>
  filter(hvg == TRUE) |>
  arrange(perc_hvg)
head(t)
dim(t)


#_____








# 2. DROPLET

#basic_df_droplet <- assemble_TMS_df_droplet(write = TRUE)

# takes saved dfs and adds gene sets etc. 
#df_main_filtered <- load_filter_og_df(data_source = "droplet")

# tag genes in which gene sets they are a part of
df_hk <- tag_hk_genes(data_source = "droplet")
df_hk_lin <- tag_hk_lin_genes(df_hk)
#df_response_genes <- tag_innate_response_genes(df_hk_lin) # too many
df_lps <- tag_lps_genes(df_hk_lin)
#df_marker <- tag_mac_marker_genes(df_lps) # too few
df_chemokine <- tag_chemokine_genes(df_lps)
df_tlr <- tag_tlr_genes(df_chemokine)
df_xue <- tag_xue_genes(df_tlr)
df_autoimmune <- tag_autoimmune_genes(df_xue)
df_sperm <- tag_sperm_genes(df_autoimmune)
df_meiosis <- tag_meiosis_genes(df_sperm)
df_oocyte <- tag_oocyte_genes(df_meiosis)
#df_all <- tag_not_expresssed_genes(df_oocyte, data_source = "droplet")


# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(df_oocyte)
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "droplet")


# stratification
strat_df <- stratify_df(data_source = "droplet", cell_type_selection = "Macrophage") # also take df as arg
strat_df <- readRDS("droplet/data/strat_df.rds")

# expression bins
e <- plot_expression_bins(strat_df, "droplet")
#l <- plot_expression_bins_lps(strat_df, "droplet")
n <- plot_expression_bin_numbers(strat_df, "droplet")

# categories
h <- plot_category_numbers(strat_df, "droplet")
g <- plot_gene_set_numbers(strat_df, "droplet")

d <- plot_gene_set_proportions(strat_df, "droplet")
j <- gene_set_proportions_compact(strat_df, "droplet")

# qc
mac_marker_expression(df = strat_df, data_source = "droplet")


# gene distribution
selection_gene_distr_df <- gene_distribution(strat_df,'droplet') 
df_dummy <- dummy_entropy("droplet")
all_genes_rao_prep <- all_genes_entropy_prep(strat_df, 'droplet')
subsampled_df <- subsample_genes(all_genes_rao_prep) # faster entropy calc
df_entropy <- quadratic_entropy(all_genes_rao_prep, 'droplet', write = TRUE) # or dummy data or selected genes data or subsampled data 
plot_raos_gene_selection(df_entropy, "droplet")
plot_raos_all_genes(df_entropy, "droplet")

master_entropy_df <- merge_entropy_results("droplet")
variability_df <- calc_variability_direction(master_entropy_df, "droplet")
#scatter_entropy_hvg(variability_df, "droplet")
heatmap_entropy_bins(master_entropy_df, "droplet")

colnames(strat_df)
unique(master_entropy_df$cluster_id)

t <- master_entropy_df |>
  filter(gene_variability == "HVG") |>
  filter(res_var < 5)
t

# test how many hvg are there
t <- master_entropy_df |>
  group_by(gene_variability) |>
  count()
t

