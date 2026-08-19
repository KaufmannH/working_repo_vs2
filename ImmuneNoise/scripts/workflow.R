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
library(writexl)
#library(org.Mm.eg.db)
#library(rtracklayer)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(ggridges)
library(RColorBrewer)
library(ggh4x)
library(clusterCrit)
library(dynamicTreeCut)
library(hexbin)
library(viridis)
library(future)
#library(progressr)
#library(beepr)


# load functions
source("ImmuneNoise/scripts/src/assemble_TMS_df.R")
source("ImmuneNoise/scripts/src/assemble_pansci_df.R")
source("ImmuneNoise/scripts/src/load_gene_sets.R")
source("ImmuneNoise/scripts/src/stratification.R")

source("ImmuneNoise/scripts/src/gene_distribution.R")
source("ImmuneNoise/scripts/src/qc.R")
source("ImmuneNoise/scripts/src/compare_technologies.R")
 
#----------




# 1. DROPLET 

#test
d <- read_csv("ImmuneNoise/droplet/data/df_main_filtered.csv")
m <- read_tsv("ImmuneNoise/droplet/droplet_raw_data/full_results_10_22.tsv")
m

# what genes still have super low var
low_resvar_genes <- m |>
    filter(res_var == 0) |>
    select(gene, cluster_id, res_var, gmean) 
low_resvar_genes



# load functions
source('ImmuneNoise/scripts/src/assemble_additional_spleen.R')
source("ImmuneNoise/scripts/src/load_gene_sets.R")
source("ImmuneNoise/scripts/src/stratification.R")
#basic_df_droplet <- assemble_TMS_df_droplet(write = TRUE)
df_main_all <- read.csv(paste0('ImmuneNoise/droplet/data/combined_data.csv'))
unique(df_main_all$cell_type)

# takes saved dfs and adds gene sets etc. 
selected_cell_type_list <- c( 'NKT','NK')
selected_cell_type_list <- c('epithelial_basal', 'endothelial', 'fibroblast', 'epithelial', 'epithelial_luminal', 'epithelial_bile_duct')
selected_cell_type_list <- c( 'macrophage','myeloid', 'monocyte', 'dc' )
df_main_filtered <- load_filter_og_df(data_source = "droplet", selected_cell_type_list = selected_cell_type_list) #tissue_selection = 'Spleen'




# tag genes in which gene sets they are a part of
df_hk <- tag_hk_genes(data_source = "droplet") #selected_cell_type = "Hepatocyte"
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


# short version (for doc): 
df_hk <- tag_hk_genes(data_source = "droplet")
df_mac <- tag_mac_immune_response(df_hk)
df_control <- tag_control_genes(df_mac)
#df_fasted <- tag_fasted_genes(df_control)
df_prot <- tag_protein_coding_genes(df_control)
df_mirna <- tag_miRNA_genes(df_prot)
df_lncrna <- tag_lncRNA_genes(df_mirna)
df_ribo <- tag_ribozyme_genes(df_lncrna)



# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(df_ribo) # or oocyte
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "droplet")


# stratification
strat_df <- stratify_df(data_source = "droplet", cell_type_selection = NULL) # REVERSE # also take df as arg 
strat_df <- readRDS("ImmuneNoise/droplet/data/strat_df.rds")

subsampled_strat_df <- readRDS("droplet/data/strat_subsampled_df.rds")
unique(strat_df$cell_type)

# expression bins
e <- plot_expression_bins(strat_df, "droplet")
#l <- plot_expression_bins_lps(strat_df, "droplet")
n <- plot_expression_bin_numbers(strat_df, "droplet")

# categories
h <- plot_category_numbers(strat_df, "droplet")
g <- plot_gene_set_numbers(strat_df, "droplet")
f <- plot_all_clusters(strat_df, "droplet")

d <- plot_gene_set_proportions(strat_df, "droplet")
j <- gene_set_proportions_compact(strat_df, "droplet")
i <- gene_set_proportions_big_groups(strat_df, "droplet")

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









# 2. FACS (right now: spleen)

# load functions
source("ImmuneNoise/scripts/src/assemble_TMS_df.R")
source("ImmuneNoise/scripts/src/assemble_pansci_df.R")
source("ImmuneNoise/scripts/src/load_gene_sets.R")
source("ImmuneNoise/scripts/src/stratification.R")

#basic_df_facs <- assemble_TMS_df_facs(write = TRUE)

# takes saved dfs and adds gene sets etc. 
selected_cell_type_list <- c('monocyte', 'macrophage', 'myeloid cell', 'Kupffer cell', 'myeloid leukocyte', 'myeloid dendritic cell', 'classical monocyte', 'dendritic cell')
df_main_filtered <- load_filter_og_df(data_source = "facs", selected_cell_type_list = selected_cell_type_list)


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


# short version (for doc): 
df_hk <- tag_hk_genes(data_source = "facs")
df_mac <- tag_mac_immune_response(df_hk)
df_control <- tag_control_genes(df_mac)


# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(df_control)
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "facs")

# stratification
strat_df <- stratify_df(data_source = "facs") # also take df as arg
strat_df <- readRDS("ImmuneNoise/facs/data/strat_df.rds")


# expression bins
e <- plot_expression_bins(strat_df, "facs")
#l <- plot_expression_bins_lps(strat_df, "facs")
n <- plot_expression_bin_numbers(strat_df, "facs")

# categories
h <- plot_category_numbers(strat_df, "facs")
g <- plot_gene_set_numbers(strat_df, "facs")
f <- plot_all_clusters(strat_df, "facs")

d <- plot_gene_set_proportions(strat_df, "facs")
j <- gene_set_proportions_compact(strat_df, "facs")
i <- gene_set_proportions_big_groups(strat_df, "facs")

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





# 3. ADDITIONAL SPLEEN 

basic_df_spleen_10X <-  assemble_spleen_10X(write = TRUE)

# takes saved dfs and adds gene sets etc. 
df_main_filtered <- load_filter_og_df(data_source = "spleen_10X", selected_cell_type_list = c("Myeloid")) #selected_cell_type = "Hepatocyte"


# tag genes in which gene sets they are a part of
df_hk <- tag_hk_genes(data_source = "spleen_10X", selected_tissue = "Spleen" ) #selected_cell_type = "Hepatocyte"
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


# short version (for doc): 
df_hk <- tag_hk_genes(data_source = "spleen_10X")
df_mac <- tag_mac_immune_response(df_hk)
df_control <- tag_control_genes(df_mac)


# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(df_control) # or oocyte
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "spleen_10X")

  unique(tagged_df$cell_type)


# stratification
strat_df <- stratify_df(data_source = "spleen_10X", cell_type_selection = NULL) # REVERSE # also take df as arg 
strat_df <- readRDS("ImmuneNoise/spleen_10X/data/strat_df.rds")

subsampled_strat_df <- readRDS("spleen_10X/data/strat_subsampled_df.rds")
unique(strat_df$cell_type)

# expression bins
e <- plot_expression_bins(strat_df, "spleen_10X")
#l <- plot_expression_bins_lps(strat_df, "droplet")
n <- plot_expression_bin_numbers(strat_df, "spleen_10X")

# categories
h <- plot_category_numbers(strat_df, "spleen_10X")
g <- plot_gene_set_numbers(strat_df, "spleen_10X")
f <- plot_all_clusters(strat_df, "spleen_10X")

d <- plot_gene_set_proportions(strat_df, "spleen_10X")
j <- gene_set_proportions_compact(strat_df, "spleen_10X")
gene_set_proportions_big_groups(strat_df, "spleen_10X")

# qc
mac_marker_expression(df = strat_df, data_source = "spleen_10X")


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


# testing CCL2 in old vs young spleen

strat_df <- read_csv("/home/hkaufm49/analyses/working_repo/ImmuneNoise/droplet/data/combined_data.csv")

t <- strat_df |>
  filter(if_any(starts_with('Ccl'), ~ .x > 0)) |>
  select(age, lvg, gmean, starts_with('Ccl')) |>
  group_by(age) |>
  summarize(mean(gmean))
t

unique(strat_df$cell_type)



plot <- ggplot(strat_df, aes(x = gmean, y = res_var, color = age)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_smooth(method = "loess", color = "firebrick", linewidth = 0.8, se = TRUE) +
  scale_x_log10() +
  labs(
    x = "Mean expression per cluster (log10 scale)",
    y = "Residual variance",
    title = "Mean expression vs. residual variance") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 20)

  ggsave('ImmuneNoise/droplet/plots/all_ages/gmean_resvar_corr.png', plot)





# 4. EASYSCI
# right now only hepatocytes

# assemble df
bs_stats(tissue = 'Liver', bs_path = "RegulatedNoiseHel/data/processed/pansci_liver/bootstrapping_tables")
transfer_cell_type_annotation( so_path = "RegulatedNoiseHel/data/processed/pansci_liver/seurat_objects")
annotate_gene_names()
assemble_easysci_df(write = TRUE)

# takes saved dfs and adds gene sets etc. 
df_main_filtered <- load_filter_og_df(data_source = "pansci")
head(df_main_filtered)

# tag genes in which gene sets they are a part of

# short version (for doc): 
df_hk <- tag_hk_genes(data_source = "pansci")
df_mac <- tag_mac_immune_response(df_hk)
df_control <- tag_control_genes(df_mac)

# cleaning, renmaing etc. 
#report_genes_percluster <- genes_per_cluster(df_autoimmune)
no_na <- set_NA_false(df_control) # or oocyte
prep <- set_gene_variability(no_na)
renamed <- rename_cluster_id(prep)
renamed_df <- rename_gene_sets(renamed)$df
gene_sets <- rename_gene_sets(renamed)$existing_gene_sets
tagged_df <- tag_no_gene_set_genes(renamed_df, gene_sets, "pansci")


source("ImmuneNoise/scripts/src/stratification.R")
# stratification
strat_df <- stratify_df(data_source = "pansci") # also takes df as arg
strat_df <- readRDS("ImmuneNoise/pansci/data/strat_df.rds")
#subsampled_strat_df <- readRDS("ImmuneNoise/pansci/data/strat_subsampled_df.rds")

head(strat_df)
# see how many hvgs there are
t <- strat_df |>
  group_by( hvg) |>
  count()
t


# expression bins
e <- plot_expression_bins(strat_df, "pansci")
#l <- plot_expression_bins_lps(strat_df, "pansci")
n <- plot_expression_bin_numbers(strat_df, "pansci")

# categories
h <- plot_category_numbers(strat_df, "pansci")
g <- plot_gene_set_numbers(strat_df, "pansci") #todo
f <- plot_all_clusters(strat_df, "pansci") #todo

d <- plot_gene_set_proportions(strat_df, "pansci")
j <- gene_set_proportions_compact(strat_df, "pansci")


# qc
#mac_marker_expression(df = strat_df, data_source = "pansci")

# TODO: implement this part
# gene distribution
selection_gene_distr_df <- gene_distribution(strat_df,'pansci') 
df_dummy <- dummy_entropy("pansci")
all_genes_rao_prep <- all_genes_entropy_prep(strat_df, 'pansci')
subsampled_df <- subsample_genes(all_genes_rao_prep) # faster entropy calc
df_entropy <- quadratic_entropy(all_genes_rao_prep, 'pansci', write = TRUE) # or dummy data or selected genes data or subsampled data 
plot_raos_gene_selection(df_entropy, "pansci")
plot_raos_all_genes(df_entropy, "pansci")

master_entropy_df <- merge_entropy_results("pansci")
variability_df <- calc_variability_direction(master_entropy_df, "pansci")
#scatter_entropy_hvg(variability_df, "pansci")
heatmap_entropy_bins(master_entropy_df, "pansci")




