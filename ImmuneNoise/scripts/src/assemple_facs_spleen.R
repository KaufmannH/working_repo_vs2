# this is for checking if facs analysis was fine

library(Seurat)
library(dplyr)
library(data.table)


data_dir <- "/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/Spleen_facs/raw"   
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
length(files)  


head(read.table(files[1], sep = ",", header = FALSE))

read_cell <- function(f) {
  dt <- fread(f, header = FALSE, col.names = c("gene", "count"))
  setNames(dt$count, dt$gene)
}


cell_list <- lapply(files, read_cell)
names(cell_list) <- tools::file_path_sans_ext(basename(files))


gene_sets <- lapply(cell_list[1:min(10, length(cell_list))], names)
all_same <- all(sapply(gene_sets[-1], identical, gene_sets[[1]]))
all_same  

if (all_same) {
  counts <- do.call(cbind, cell_list)
  rownames(counts) <- names(cell_list[[1]])
} else {
  all_genes <- Reduce(union, lapply(cell_list, names))
  counts <- sapply(cell_list, function(x) {
    v <- x[all_genes]; v[is.na(v)] <- 0; v
  })
  rownames(counts) <- all_genes
}
dim(counts)  # genes x cells

meta <- data.frame(
  cell_id  = colnames(counts),
  well     = sub("^([A-P][0-9]+)-.*",                "\\1", colnames(counts)),
  plate    = sub("^[A-P][0-9]+-([^-]+)-.*",          "\\1", colnames(counts)),
  mouse.id = sub("^[A-P][0-9]+-[^-]+-([^-]+)-.*",    "\\1", colnames(counts)),
  sex      = sub(".*_([MF])-.*",                     "\\1", colnames(counts)),
  tissue   = "Spleen",
  cell_ontology_class = "Splenocytes",
  stringsAsFactors = FALSE
)

# Derived columns
meta$age     <- paste0(sub("^([0-9]+)_.*", "\\1", meta$mouse.id), "m")   # <-- mouse.id
meta$sex     <- ifelse(meta$sex == "F", "female",
                ifelse(meta$sex == "M", "male", NA))
meta$age_sex <- paste0(meta$age, "_", meta$sex)

rownames(meta) <- meta$cell_id
head(meta)
table(meta$plate)
table(meta$age_sex)


spleen <- CreateSeuratObject(
  counts       = counts,
  project      = "TabulaMuris_Spleen_SS2",
  meta.data    = meta,
  min.cells    = 3,
  min.features = 500)
spleen

head(spleen@meta.data)



# change assay for pipeline

spleen <- SeuratObject::RenameAssays(spleen, RNA = "originalexp")
DefaultAssay(object = spleen) <- "originalexp"

spleen[["originalexp"]] <- as(spleen[["originalexp"]], "Assay")
class(spleen[["originalexp"]])   # should be "Assay" (not "Assay5")
DefaultAssay(spleen)             # should still be "originalexp"
spleen

saveRDS(spleen, '/home/hkaufm49/analyses/working_repo/ImmuneNoise/facs/data/Spleen_facs/prepped_for_pipeline/spleen_facs.rds')



dim(spleen@meta.data)