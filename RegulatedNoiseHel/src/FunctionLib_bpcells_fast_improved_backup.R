# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
#library(scDblFinder)
library(scater)
library(philentropy)
library(foreach)
library(doFuture)
library(doRNG)
library(doParallel)


 

MakeEmbsClusters.df <- function(SeuratObj , EmbsType= 'pca')
{
  # It makes a data.frame which contains Embeddings & Clusters
  Clusters <- SeuratObj[['seurat_clusters']]
  Embs     <- Embeddings(object = SeuratObj[[EmbsType]])
  EmbsClusters <- merge(Clusters,Embs,by=0, all=TRUE)[,-1]
  rownames(EmbsClusters) <- rownames(Embs)
  return(EmbsClusters)
}


FindCentroids      <- function(EmbsClusters)
{
  # It compute the centroids in each cluster
  Centroids <- aggregate(as.matrix(EmbsClusters[,-1]) ~ seurat_clusters,EmbsClusters,mean)
  return(Centroids)
}


FindUnstableCells  <- function(EmbsClusters, Centroids)
{
  # It compute the distance between each cell and each cluster, 
  # if the minimal distance is NOT with the belonging cluster, than it is an unstable cell.
  
  Embs_euDis = data.frame( cell_id = c(),
                           cell_cluster = c(),
                           euclidian = c(), 
                           centroid_cluster=c())
  
  rowCounter = 0
  
  for (cell in 1: nrow(EmbsClusters)){
    
    if (cell %% 100 == 0) {
      
      print(paste(cell, "cells processed"))
      
    }
    
    
    
    cell_emb <- EmbsClusters[cell,-1]
    
    df <- data.frame( cell_id = c(),
                      cell_cluster = c(),
                      euclidian = c(), 
                      centroid_cluster=c())
    
    for (c in 1:nrow(Centroids[,-1]))
    {
      
      rowCounter <- rowCounter + 1
      philentropy:: distance(rbind(cell_emb, Centroids[c,-1]), mute.message = TRUE) -> euDis
      
      d <- data.frame( cell_id = rownames(cell_emb),
                       cell_cluster = EmbsClusters[cell,1],
                       euclidian = euDis ,
                       centroid_cluster = Centroids[c,1],
                       row.names = rowCounter)
      
      df <- rbind(df,d)
      
    }
    df <- df[df$euclidian == min(df$euclidian),]
    Embs_euDis <- rbind(Embs_euDis,df)
  }
  UnstableCells <- Embs_euDis[Embs_euDis$cell_cluster != Embs_euDis$centroid_cluster,]
  return(UnstableCells)
}


AssignUnstableCellFlag <- function(so){
  
  EmbsClusters  <- MakeEmbsClusters.df(so , EmbsType= 'pca')
  Centroids     <- FindCentroids(EmbsClusters)
  
  UnstableCells <- FindUnstableCells(EmbsClusters, Centroids)
  
  so@meta.data$UnstableCells <- row.names(so@meta.data) %in% UnstableCells$cell_id
  return(so)
  
}



# =============================================================================
# Seurat 5 compatible res_var functions
# Requires: Seurat, SeuratObject, sctransform, Matrix, tibble, dplyr, foreach
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: reconstruct a vst-like list from a Seurat5 SCTModel
# sctransform::get_residual_var() expects a list with specific slots;
# Seurat 5 stores the model info in SCTModel objects instead of misc$vst.out
# -----------------------------------------------------------------------------
.make_vst_like_from_seurat5 <- function(so, assay = "SCT", model = NULL) {

  if (is.null(model)) {
    model <- levels(so[[assay]])[1]
  }

  cell_attr <- Seurat::SCTResults(so, assay = assay, slot = "cell.attributes", model = model)
  feat_attr <- Seurat::SCTResults(so, assay = assay, slot = "feature.attributes", model = model)
  model_str <- Seurat::SCTResults(so, assay = assay, slot = "model", model = model)
  arguments <- Seurat::SCTResults(so, assay = assay, slot = "arguments", model = model)
  clips     <- Seurat::SCTResults(so, assay = assay, slot = "clips", model = model)

  # sctransform expects "gene" / "log_gene" in cell_attr;
  # Seurat 5 stores these as "feature" / "log_feature"
  renames <- c(
    "feature"         = "gene",
    "log_feature"     = "log_gene",
    "nfeature"        = "ngene",
    "umi_per_feature" = "umi_per_gene",
    "log_umi_per_feature" = "log_umi_per_gene"
  )
  cn <- colnames(cell_attr)
  for (old_name in names(renames)) {
    cn[cn == old_name] <- renames[old_name]
  }
  colnames(cell_attr) <- cn

  # Split feature.attributes into gene_attr (metadata) vs model_pars_fit

  known_feat_cols <- c(
    "detection_rate", "gmean", "variance",
    "residual_mean", "residual_variance",
    "genes_log_gmean_step1"
  )
  step1_cols <- grep("^step1_", colnames(feat_attr), value = TRUE)

  model_par_cols <- setdiff(colnames(feat_attr), c(known_feat_cols, step1_cols))

  gene_attr      <- feat_attr[, intersect(colnames(feat_attr), known_feat_cols), drop = FALSE]
  model_pars_fit <- as.matrix(feat_attr[, model_par_cols, drop = FALSE])

  model_pars <- NULL
  if (length(step1_cols) > 0) {
    model_pars <- as.matrix(feat_attr[, step1_cols, drop = FALSE])
    colnames(model_pars) <- sub("^step1_", "", colnames(model_pars))
  }

  vst_like <- list(
    cell_attr        = cell_attr,
    gene_attr        = gene_attr,
    model_pars_fit   = model_pars_fit,
    model_pars       = model_pars,
    model_pars_nonreg = NULL,
    model_str        = model_str,
    model_str_nonreg = NULL,
    arguments        = if (is.null(arguments)) list() else arguments
  )

  # preserve clip ranges if present
  if (is.null(vst_like$arguments$res_clip_range) && !is.null(clips$vst)) {
    vst_like$arguments$res_clip_range <- clips$vst
  }
  if (is.null(vst_like$arguments$sct.clip.range) && !is.null(clips$sct)) {
    vst_like$arguments$sct.clip.range <- clips$sct
  }

  return(vst_like)
}


# -----------------------------------------------------------------------------
# Helper: coerce BPCells / Assay5 layer to a matrix type sctransform accepts
# -----------------------------------------------------------------------------
.as_umi_matrix <- function(x) {
  if (inherits(x, "dgCMatrix") || is.matrix(x)) return(x)

  x_sparse <- tryCatch(as(x, "dgCMatrix"), error = function(e) NULL)
  if (!is.null(x_sparse)) return(x_sparse)

  x_dense <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (!is.null(x_dense)) return(x_dense)

  stop("Could not coerce counts slice to matrix/dgCMatrix.")
}


# =============================================================================
# res_var_subset
# Calculate residual variance and geometric mean per gene for a subset of cells
# using the SCT model already fit on the full dataset.
# =============================================================================
res_var_subset <- function(so, cells_filter, count_slot_name,
                           sct_assay = "SCT") {

  # resolve cell names
  sel_cells <- if (is.logical(cells_filter)) {
    colnames(so)[cells_filter]
  } else {
    as.character(cells_filter)
  }

  if (length(sel_cells) < 2L) {
    stop("Need at least 2 cells in subset.")
  }

  # --- raw counts for subset (BPCells-safe) ---
  counts_layer <- SeuratObject::LayerData(so[[count_slot_name]], layer = "counts")
  umi.subset   <- .as_umi_matrix(counts_layer[, sel_cells, drop = FALSE])

  # --- reconstruct vst-like object from Seurat5 SCT model ---
  vst_out.subset <- .make_vst_like_from_seurat5(so, assay = sct_assay)

  # subset cell_attr to selected cells
  missing_cells <- setdiff(sel_cells, rownames(vst_out.subset$cell_attr))
  if (length(missing_cells) > 0L) {
    stop("Some cells missing in SCT cell.attributes: ",
         paste(head(missing_cells, 5), collapse = ", "))
  }
  vst_out.subset$cell_attr <- vst_out.subset$cell_attr[sel_cells, , drop = FALSE]

  # --- residual variance via sctransform (same as OG) ---
  getResVar.subset <- sctransform::get_residual_var(
    vst_out       = vst_out.subset,
    umi           = umi.subset,
    res_clip_range = c(-sqrt(ncol(umi.subset)), sqrt(ncol(umi.subset)))
  )

  # --- geometric mean: mean of log1p-normalised data across subset ---
  data_layer <- SeuratObject::LayerData(so[[sct_assay]], layer = "data")
  gmean_sub  <- data_layer[, sel_cells, drop = FALSE]

  if (inherits(gmean_sub, "Matrix")) {
    getgMean.subset <- Matrix::rowMeans(gmean_sub)
  } else {
    getgMean.subset <- rowMeans(as.matrix(gmean_sub))
  }
  names(getgMean.subset) <- rownames(gmean_sub)

  # --- output ---
  res_var_mean <- cbind(
    res_var = getResVar.subset,
    gmean   = getgMean.subset[names(getResVar.subset)]
  )

  return(res_var_mean)
}


# =============================================================================
# res_var_by_group
# Re-estimate variability and gmean per gene across cells in each group
# (defined by a metadata column). Returns a tibble.
# =============================================================================
res_var_by_group <- function(so, meta_col, hvg_cutoff = 5,
                             count_slot_name, n_cores = 1,
                             sct_assay = "SCT") {

  cl_vec <- so@meta.data[[meta_col]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)

  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    `%op%` <- foreach::`%dopar%`
  } else {
    `%op%` <- foreach::`%do%`
  }

  res_var_tbl <- foreach::foreach(
    cl = levels(cl_vec),
    .combine  = "rbind",
    .packages = c("tibble", "dplyr", "SeuratObject", "Seurat",
                   "sctransform", "Matrix"),
    .export   = c("res_var_subset", ".make_vst_like_from_seurat5", ".as_umi_matrix")
  ) %op% {

    cells_filter <- cl_vec == cl
    n_cells      <- sum(cells_filter)

    if (n_cells < 2L) {
      message("WARNING: Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
      return(NULL)
    }

    res_var_mean_subs <- res_var_subset(
      so              = so,
      cells_filter    = cells_filter,
      count_slot_name = count_slot_name,
      sct_assay       = sct_assay
    )

    tibble::as_tibble(res_var_mean_subs, rownames = "gene") %>%
      dplyr::mutate(cluster = cl)
  }

  # empty fallback
  if (is.null(res_var_tbl) || nrow(res_var_tbl) == 0L) {
    return(tibble::tibble(
      gene    = character(),
      res_var = numeric(),
      gmean   = numeric(),
      cluster = character(),
      hvg     = logical()
    ))
  }

  res_var_tbl <- res_var_tbl %>%
    dplyr::mutate(hvg = res_var > hvg_cutoff)

  return(res_var_tbl)
}


# =============================================================================
# res_var_bootstrap
# Bootstrap residual variance by resampling cells WITH REPLACEMENT within each
# cluster, then calling sctransform::get_residual_var() on the resampled raw
# counts each time (faithful to the original Seurat4 approach).
# Parallelisation is across clusters.
# =============================================================================
res_var_bootstrap <- function(so, clusters_slot, n_cycles = 1000,
                              n_cores = 1, hvg_cutoff = 5,
                              res_var_cl, keep_all = FALSE,
                              count_slot_name = "originalexp",
                              sct_assay = "SCT") {

  # --- build the vst-like object once (shared across all clusters) ---
  vst_out_full <- .make_vst_like_from_seurat5(so, assay = sct_assay)

  # --- pull full raw counts once (materialise if BPCells) ---
  counts_full <- SeuratObject::LayerData(so[[count_slot_name]], layer = "counts")
  counts_full <- .as_umi_matrix(counts_full)

  # --- cluster vector ---
  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # --- parallel setup (across clusters) ---
  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    `%op%` <- foreach::`%dopar%`
  } else {
    `%op%` <- foreach::`%do%`
  }

  set.seed(42)

  out.tibble <- foreach::foreach(
    cl = cl_levels,
    .combine  = "rbind",
    .packages = c("tibble", "dplyr", "sctransform", "Matrix"),
    .export   = c("vst_out_full", "counts_full", "cl_vec",
                   "n_cycles", "hvg_cutoff", "res_var_cl", "keep_all")
  ) %op% {

    message("Processing cluster ", cl)

    # cells in this cluster
    cell_names <- colnames(so)[cl_vec == cl]
    n_cells    <- length(cell_names)

    if (n_cells < 2L) {
      message("WARNING: Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
      return(NULL)
    }

    # hvgs for this cluster (used for filtering output if keep_all = FALSE)
    hvgs_cl <- res_var_cl %>%
      dplyr::filter(cluster == cl, hvg) %>%
      dplyr::pull(gene)

    # --- run n_cycles bootstrap iterations ---
    iter_results <- vector("list", n_cycles)

    for (iter in seq_len(n_cycles)) {

      # resample cells with replacement (same as OG)
      selected.cells <- sample(cell_names, size = n_cells, replace = TRUE)

      # slice raw counts for resampled cells
      umi.subset <- counts_full[, selected.cells, drop = FALSE]

      # slice cell_attr for resampled cells
      vst_out.iter           <- vst_out_full
      vst_out.iter$cell_attr <- vst_out_full$cell_attr[selected.cells, , drop = FALSE]

      # compute residual variance from the model on resampled counts
      ResVar <- sctransform::get_residual_var(
        vst_out = vst_out.iter,
        umi     = umi.subset
      )

      # build per-iteration result
      features    <- names(ResVar)
      boolean_hvg <- ResVar >= hvg_cutoff

      iter_results[[iter]] <- tibble::tibble(
        gene        = features,
        boolean_hvg = boolean_hvg,
        ResVar      = as.numeric(ResVar),
        iter        = iter
      )
    }

    # combine all iterations for this cluster
    tmp.tibble <- dplyr::bind_rows(iter_results)
    tmp.tibble$cluster <- cl

    # filter to hvgs only (unless keep_all)
    if (!keep_all) {
      tmp.tibble <- tmp.tibble %>%
        dplyr::filter(gene %in% hvgs_cl)
    }

    tmp.tibble
  }

  # empty fallback
  if (is.null(out.tibble) || nrow(out.tibble) == 0L) {
    return(tibble::tibble(
      gene        = character(),
      boolean_hvg = logical(),
      ResVar      = numeric(),
      iter        = integer(),
      cluster     = character()
    ))
  }

  return(out.tibble)
}



plot_resvar <- function(res_var_cl, bootstrap_res, test_gene, housek_gene, plot_cluster) {
  
  resvar_test_gene <- res_var_cl %>% 
    filter(cluster == plot_cluster, gene == test_gene) %>% pull(res_var)
  
  resvar_plot <- ggplot(subset(bootstrap_res, gene %in% c(test_gene, housek_gene) & cluster == plot_cluster), 
                        aes(x = ResVar, y = gene)) +
    geom_boxplot() +
    geom_point(aes(x = resvar_test_gene, y = test_gene), colour = "red", size = 5, shape = 16)
  
  return(resvar_plot)
}




# An adaptation of res_var_bootstrap to calculate stability for stable genes
# This function takes as an input a seurat object, clusters and a number of cycles and will
# then perform that many cycles of bootstrapping by choosing random subsets of cells
# and re-calculating the variability for each gene in those subsets
res_var_bootstrap_low <- function(so, clusters_slot, n_cycles = 1000, n_cores, lvg_cutoff = 1, res_var_cl, keep_all = FALSE){
  
  n_cl <- 0
  
  for (cl in levels(so@meta.data[[clusters_slot]])){
    
    print(cl)
    
    # Get Number of Hyper Variable Genes
    features_tbl <- rownames(so@assays$originalexp@meta.features) %>% 
      as_tibble() %>% dplyr::rename(features = value)
    
    # 1. Cell filter for cluster
    cells_filter <- so@meta.data[[clusters_slot]] == cl 
    cell_names <- colnames(so)[cells_filter]
    
    # 2. Sample Cells in each Cluster 1000 with replacement
    
    set.seed(42) # Set the seed to be able to reproduce the same set of sampling
    
    # setup a cluster for parallelization 
    #cluster <- makeCluster(n_cores)
    #registerDoParallel(cluster)
    
    #registerDoFuture()
    #plan(multisession, workers = n_cores)
    
    out.matrix <- foreach(iter = 1:n_cycles, .combine = 'rbind', .export = ls(globalenv())) %do% {
      
      selected.cells <- sample(cell_names, size = length(cell_names), replace = TRUE)
      
      umi.subset     <- so@assays$originalexp@counts[,selected.cells]
      vst_out.subset <- so@assays$SCT@misc$vst.out
      vst_out.subset$cell_attr <- vst_out.subset$cell_attr[selected.cells,]
      
      ResVar         <- sctransform::get_residual_var(vst_out = vst_out.subset, umi = umi.subset)
      
      # build output matrix
      features <- names(ResVar)
      boolean_lvg <- features %in% names(ResVar[ResVar <= lvg_cutoff])
      
      cbind(boolean_lvg, ResVar, iter)
      
    }
    
    #stopCluster(cluster)
    
    # lvgs for cluster
    lvgs <- res_var_cl %>% filter(cluster == cl,
                                  res_var <= lvg_cutoff)
    
    # filter results for lvgs
    if (keep_all) {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene")
    } else {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene") %>%
        filter(gene %in% lvgs$gene)
    }
    
    rm(out.matrix)
    gc()
    
    tmp.tibble$cluster <- cl
    
    # combine result tibbles
    if (n_cl < 1) {
      out.tibble <- tmp.tibble
    } else {
      out.tibble <- bind_rows(out.tibble, tmp.tibble)
    }
    n_cl <- n_cl + 1
    
  }
  
  return(out.tibble)
}



# function for sampling LVGs to compare with HVGs (still needs weighting adaptation)
sample_LVG <- function(res_var_cl, HVG_cutoff = 5, LVG_cutoff = 1) {
  
  count = 0
  
  for (cl in unique(res_var_cl$cluster)) {
    
    hvg <- res_var_cl %>% filter(cluster == cl, res_var >= HVG_cutoff)
    lvg <- res_var_cl %>% filter(cluster == cl, res_var <= LVG_cutoff)
    
    lvg_sample <- sample(lvg$gene, size = nrow(hvg), replace = FALSE)
    lvg_subset <- lvg %>% filter(gene %in% lvg_sample)
    
    if (count == 0) {
      lvg_subset_tbl <- lvg_subset
    } else {
      lvg_subset_tbl <- rbind(lvg_subset_tbl, lvg_subset)
    }
    
    count = count + 1
    
  }
  
  return(lvg_subset_tbl)
  
}