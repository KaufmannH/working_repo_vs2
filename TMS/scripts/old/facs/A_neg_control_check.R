# take df from exploring distributions (test 6) and investigate the genes from the negative control.
# ok stefan said its not necessary

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(purrr)
library(gprofiler2)

df_main <- readRDS("facs_analysis/data/df_with_markers.rds")
colnames(df_main)

sum(df_main$sperm_gene)

# filter for all genes of the negative control that are expressed in the data set
genes_expressed <- df_main |>
  select(-expression_level) |>
  filter(sperm_gene == TRUE | meiosis_gene == TRUE) |>
  filter(gmean > -1 ) |>
  pull(gene) |>
  unique()
genes_expressed

all_sperm_sign  <- read_excel('facs_analysis/data/GO_gene_lists/GO_sperm_dna_condensation.xlsx')
sperm_list <- all_sperm_sign |> pull(Symbol)

all_meiosis_signature <- read_excel('facs_analysis/data/GO_gene_lists/GO_G1_M_transition_meiosis.xlsx')
meiosis_list <- all_meiosis_signature |> pull(Symbol)

fuse_sign <- unique(union(sperm_list, meiosis_list))
length(fuse_sign)
# there are 28 genes in the signatures

overlap <- intersect(fuse_sign, genes_expressed)
length(overlap)
# of thee 28 available genes there are 21 detected in the data set
print(overlap)


# now check them for how general they are with gprofiler
overlap_genes_pathways <- gost(  query = list(overlap = overlap), organism = "mmusculus", sources = c("KEGG", "REAC"),  evcodes = TRUE)
res <- overlap_genes_pathways$result
res

# they are all associated with cell cycle, and replication so not sperm specific. 

