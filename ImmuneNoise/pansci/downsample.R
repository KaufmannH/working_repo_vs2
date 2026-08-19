
library(tidyverse)
library(Seurat)
library(rtracklayer)

# for run new seurate vs2 pipeline (26.07.2026)
# changed the gene names form the start, not after


data <- readRDS('/share/NoiseControlCenter_data/data/pansci/regnoise_upstream/data/prepped_for_pipeline/Liver/Liver_3m_female_prepped.rds')

head(data@meta.data)
dim(data@meta.data)
head(rownames(data)) 
Assays(data)



annotate_gene_names_seurat <- function(seu,
                                       assay    = DefaultAssay(seu),
                                       gtf_path = "/share/NoiseControlCenter_data/data/pansci/regnoise_upstream/data/gencode.vM27.annotation.gtf") {

  # 1. get current features
  feat <- tibble(gene = rownames(seu[[assay]])) |>
    mutate(ensembl_clean = sub("\\.\\d+$", "", gene))

  # 2. ensembl  symbol map from GTF 
  gtf <- rtracklayer::import(gtf_path)
  gtf_map_raw <- mcols(gtf)[, c("gene_id", "gene_name")] |>
    as.data.frame() |>
    transmute(
      gene_id,
      ensembl_clean = sub("\\.\\d+$", "", gene_id),
      SYMBOL        = gene_name) |>
    distinct(ensembl_clean, .keep_all = TRUE)

  symbol_count <- gtf_map_raw |>
    distinct(ensembl_clean, SYMBOL) |>
    count(SYMBOL, name = "n_ens")

  # 3. get new names
  new_names <- feat |>
    left_join(gtf_map_raw, by = "ensembl_clean") |>
    left_join(symbol_count, by = "SYMBOL") |>
    mutate(new_name = case_when(
      !is.na(SYMBOL) & SYMBOL != "" & n_ens == 1 ~ SYMBOL,
      TRUE                                        ~ gene)) |>
    pull(new_name)

  # 4. guarantee uniqueness
  dup_n <- sum(duplicated(new_names))
  if (dup_n > 0) {
    message("make.unique() adjusted ", dup_n, " duplicate symbol(s)")
    new_names <- make.unique(new_names, sep = "-")
  }
  stopifnot(length(new_names) == nrow(seu[[assay]]))

  # % still unmapped
  pct_ens <- round(100 * mean(str_detect(new_names, "^ENSMUSG")), 2)
  message(sprintf("%.2f%% of features still carry an Ensembl ID (unmapped)", pct_ens))

  # 5. rebuild from renamed counts
  counts <- LayerData(seu, assay = assay, layer = "counts")
  rownames(counts) <- new_names

  seu_renamed <- CreateSeuratObject(counts    = counts,
                                    meta.data = seu[[]],
                                    assay     = assay,    
                                    project   = Project(seu))
    return(seu_renamed)                  
}


annotated_so <- annotate_gene_names_seurat(data)
head(rownames(annotated_so)) 



# adjust columns 

annotated_so@meta.data <- annotated_so@meta.data |>
  mutate(
    sex = tolower(sex),                                   
    age = paste0(as.integer(sub("_months$", "", age)), "m"),  
    age_sex = paste(age, sex, sep = "_") )

head(annotated_so@meta.data)

Assays(annotated_so) 



# downsample 

total_target <- 30000
meta <- annotated_so@meta.data |> tibble::rownames_to_column("barcode")

counts    <- meta |> count(mouse.id, name = "n_cells")
n_mice    <- nrow(counts)
per_mouse <- floor(total_target / n_mice)


keep <- meta |>
  group_by(mouse.id) |>
  slice_sample(n = per_mouse) |>      
  ungroup() |>
  pull(barcode)

seu_down <- subset(annotated_so, cells = keep)

dim(seu_down)
table(seu_down@meta.data$mouse.id)


# check assay type
Assays(seu_down) 


saveRDS(seu_down, '/home/hkaufm49/analyses/Regnoise_vs2/data/rawdata/Liver_pansci/Liver_3m_female_downsampled.rds')


