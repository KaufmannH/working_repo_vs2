  
# comparing literature derived gene sets with Immune dict (Cui et al. 2024) gene sets




# immune dict conditions
imm_dict_gene_list <- readRDS("reference_gene_sets/ImmuneDict/marker_genes_condition.rds")
names(imm_dict_gene_list)

# summarize the cytokines (bc in gene set there are )
cytokine_list_raw <- imm_dict_gene_list[c("41BBL", "APRIL", "BAFF", "Cardiotrophin-1", "CD27L", "CD30L", "CD40L",
  "FasL", "Flt3l", "G-CSF", "GITRL", "GM-CSF", "IL-Y", "IL10", "IL11", "IL12", "IL13", "IL15",
  "IL17A", "IL17B", "IL17C", "IL17D", "IL17E", "IL17F", "IL18", "IL19", "IL1a", "IL1b", "IL1ra",
  "IL2", "IL20", "IL21", "IL22", "IL23", "IL24", "IL27", "IL3", "IL30", "IL31", "IL33", "IL34",
  "IL36a", "IL36RA", "IL4", "IL5", "IL6", "IL7", "IL9", "IFNa1", "IFNb", "IFNe", "IFNg", "IFNk",
  "IFNl2", "LIF", "LIGHT", "LTA1-B2", "LTA2-B1", "M-CSF", "OSM", "OX40L")]
cytokine_list <- unlist(cytokine_list_raw, use.names = FALSE)
cytokine_list
subset_imm_dict <- list()
subset_imm_dict[["cytokines"]] <- cytokine_list
subset_imm_dict['IFNg'] <- imm_dict_gene_list['IFNg']
subset_imm_dict['IL4'] <- imm_dict_gene_list['IL4']
subset_imm_dict['IL10'] <- imm_dict_gene_list['IL10']
subset_imm_dict


# chemokine
cytokine_df <- read.csv("reference_gene_sets/InnateDB_cytokine_signalling.csv")
cytokine_gene_list <- cytokine_df |> pull(name) |> unique()
cytokine_gene_list


# xue
xue_response_raw <- readRDS("reference_gene_sets/xue_response_genes.rds")
xue_response <- xue_response_raw |> select(Condition, Upregulated)

# load the mouse/human conversions
m2h <- read.delim("reference_gene_sets/20200307_ensembl/mouse.txt",  header = FALSE, sep = "\t")
colnames(m2h) <- c("mouse", "human")
h2m <- read.delim("reference_gene_sets/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
colnames(h2m) <- c("human", "mouse")

# see if the gene is in def and then pull the ortholog
xue_response_orthologs <- xue_response |>
mutate(Orthologs = map(Upregulated, ~ h2m$mouse[match(.x, h2m$human)]))

xue_gene_list <- xue_response_orthologs |> 
select(Condition, Orthologs) |> 
unnest(Orthologs) |> 
filter(!is.na(Orthologs)) |>    
group_by(Condition) |> 
summarise(GeneList = list(unique(Orthologs))) |> 
deframe()

names(xue_gene_list)[names(xue_gene_list) == "IFNg"] <- "IFNy_response_gene"
names(xue_gene_list)[names(xue_gene_list) == "IL10"] <- "IL10_response_gene"
names(xue_gene_list)[names(xue_gene_list) == "IL4"] <- "IL4_response_gene"
names(xue_gene_list)

# summarize the literature gene sets
lit_sets <- list()
lit_sets[['cytokines']] <- cytokine_gene_list
lit_sets['IFNg'] <- xue_gene_list['IFNy_response_gene']
lit_sets['IL4'] <- xue_gene_list['IL4_response_gene']
lit_sets['IL10'] <- xue_gene_list['IL10_response_gene']
lit_sets


# overlap
conds <- intersect(names(subset_imm_dict), names(lit_sets))
conds

overlap_df <- tibble(
  Condition = conds,
  set1 = subset_imm_dict[conds],
  set2 = lit_sets[conds]) |>
  rowwise() |>
  mutate(
    OverlapPercent = length(intersect(set1, set2)) / length(union(set1, set2)) * 100) |>
  ungroup() |>
  select(Condition, OverlapPercent)
print(overlap_df)



plot <- ggplot(overlap_df, aes(x = Condition, y = OverlapPercent, fill = Condition)) +
  geom_bar(stat = "identity" ) +
  labs(y = "Percent Overlap", x = "Gene sets") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_classic()

ggsave("comparison/plots/3_m/Test_14/gene_set_overlap.png", plot = plot, width = 6, height = 4, dpi = 300)



