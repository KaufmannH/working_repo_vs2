# old = light label
# young = heavy label
# positive values: more old (slow), neg: more new (fast)


tunrover_df_raw <- read_excel("protein_turnover_sabatier.xlsx", skip = 1)
head(tunrover_df_raw)
tunrover_df_raw$PG_group_quantity <- as.numeric(tunrover_df_raw$PG_group_quantity)

# conversion of human to mouse genes
h2m <- read.delim("reference_gene_sets/20200307_ensembl/human.txt", header = FALSE, sep = "\t")
colnames(h2m) <- c("human", "mouse")
# human has many mouse matches
unique_h2m <- h2m |>
  group_by(human) |>
  slice(1) |>
  ungroup()

# calc weighted mean turnover per gene/cell_line 
# PG_group_quantity: gives weight to each measurement by protein abundance.
turnofer_df <- tunrover_df_raw |>
  left_join(unique_h2m, by = c("query_gene_names" = "human")) |>
  rename(gene = mouse) |>
  filter(!is.na(gene)) |>
  mutate(cell_line = case_when(
      str_detect(sample, "HeLa")      ~ "HeLa",
      str_detect(sample, "iPSCs?")    ~ "iPSC",
      str_detect(sample, "EBs?")      ~ "EB",
      TRUE                            ~ "other")) |>
   filter(!is.na(turnover_value), !is.na(PG_group_quantity), PG_group_quantity > 0) |>
    group_by(gene, cell_line) |> 
    summarise(mean_turnover = sum(turnover_value * PG_group_quantity, na.rm= TRUE) /
                    sum(PG_group_quantity, na.rm=TRUE),   .groups = "drop")
head(turnofer_df)
tail(turnofer_df)
sum(is.na(turnofer_df$cell_line))

# fuse with previous df
fused_df <- variability_df |>
    left_join(turnofer_df, by = 'gene') |>
    filter(!is.na(gene),
            !is.na(cell_line),
            !is.na(mean_turnover))|> # not all genes were in turnover data
     mutate(binary_var_direction = if_else(
                variability_direction > 0,
                "more_hvg",
                "more_lvg"))
head(fused_df)




data_source <- 'facs'
# are the cell line turnovers the same? 
plot <- ggplot(fused_df, aes(x = mean_turnover, y = rao_q, color = factor(cell_line))) +
  geom_jitter(width = 0.001, height = 0, size = 1, alpha = 0.5) +
  labs(
    x = "Protein turnover",
    y = "Rao's entropy") +
  theme_classic()

  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_15/scatter_turnover.png", plot, width = 7, height = 7)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_15/scatter_turnover.png", plot, width = 7, height = 7)
  } else {
    stop("Issue when saving.")
  } 



# check out gene sets and their turnover


# get ribosomal proteins and see how many of htem house keeping and then link with protein stability vs not in housekeeping

# read ribo genes
ribo_df_raw <- read_excel('GO_ribosome.xlsx')
ribo_df <- ribo_df_raw |>
    rename(gene = Symbol,
    annotation =`Annotated Term` ) |>
    select(gene, annotation)
unique(ribo_df$annotation)

ribo_list <- ribo_df |>
    pull(gene) |>
    unique()


# tag ribo genes in df
ribo_df <- strat_df |>
      mutate(`Ribosomal gene` = gene %in% ribo_list) |>
  filter(`Ribosomal gene` == TRUE) |>
  mutate(category = factor(category, levels = all_cats)) |>
  unique() 
  colnames(strat_df)
  
# add ribo and gene set data to turnover
gene_set_turnover_df <- ribo_df |>
    left_join(fused_df) |> 
    filter(!is.na(mean_turnover)) |>
    #select(gene, rao_q, mean_turnover,variability_direction, `IFNy response gene`, `LPS response gene`,`TNF response gene`, `TLR response gene` ) |>
    unique() |>
    arrange(gene) 
#(gene twice bc of cell lines)


   all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")

gene_set_order <- names(gene_set_turnover_df)[grepl(" gene$", names(gene_set_turnover_df))]

long_gene_set_df <- gene_set_turnover_df |>
     pivot_longer(
      cols      = all_of(gene_set_order),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    select(gene, rao_q, mean_turnover, gene_set, category) |>
    distinct() |>
    arrange(gene, gene_set, category)

head(long_gene_set_df)
unique(long_gene_set_df$gene_set)
t <- long_gene_set_df |>
  group_by(gene, gene_set, category) |>
  count()
t


# gene sets and their turnover
 df_l <- long_gene_set_df |>
    select(gene, mean_turnover, gene_set) |>
    unique()

plot <- ggplot(df_l, aes( x = mean_turnover, y = gene_set)) +
  geom_jitter(height = 0.2, size = 1.5, alpha = 0.7, color = "orange") +
  labs(
    x = "Mean protein turnover (log₂ L/H)",
    y = "Gene set" ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank())



  plot <- ggplot(long_gene_set_df, aes( x = mean_turnover, y = gene_set , group = gene_set)) +
   geom_density_ridges(colour = "white", size = 0.7) +
  labs(
    x = "Mean protein turnover (log₂ L/H)",
    y = "Gene set" ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank())  


  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_15/ridge_turnover_gene_set.png", plot, width = 7, height = 7)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_15/ridge_turnover_gene_set.png", plot, width = 7, height = 7)
  } else {
    stop("Issue when saving.")
  } 


# categories and their turnover

# there is for every gene there are two turnover values (cell lines)
plot_df <-long_gene_set_df |>
  select(gene, mean_turnover, category) |>
  unique()

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(plot_df, aes( x = mean_turnover, y = category , fill = category)) +
   geom_density_ridges(colour = "white", size = 0.7) +
  labs(
    x = "Mean protein turnover (log₂ L/H)",
    y = "Category" ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = col_vec) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank())  


  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_15/ridge_turnover_category.png", plot, width = 7, height = 7)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_15/ridge_turnover_category.png", plot, width = 7, height = 7)
  } else {
    stop("Issue when saving.")
  } 






 



