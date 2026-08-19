# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(scDblFinder)
library(scater)
library(philentropy)
library(foreach)
library(doFuture)
library(doRNG)
library(doParallel)


detect_doublets <- function(count_matrix, cluster_vector) {
  
  # 1. density approach
  dbl_dens.out <- computeDoubletDensity(count_matrix)
  dbl_dens.score <- tibble(cell = colnames(count_matrix), 
                           dbl_dens.score = dbl_dens.out)
  
  
  # 2. intermediate cluster approach
  dbl_cl.out <- findDoubletClusters(count_matrix, clusters = cluster_vector)
  dbl_cl.outtib <- as_tibble(dbl_cl.out)
  dbl_cl.outtib$query_cluster <- rownames(dbl_cl.out)
  
  chosen.doublet.cluster <- rownames(dbl_cl.out)[isOutlier(dbl_cl.out$num.de, type = "lower", log = TRUE)]
  
  
  # 3. output
  doublet_results <- list(doublet_density = dbl_dens.score, 
                          cluster_table = dbl_cl.outtib, 
                          cluster_outlier = chosen.doublet.cluster)
  
  return(doublet_results)
  
  
}


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
# Seurat 5 compatible res_var functions (standard matrix, no BPCells)
# Matches the Seurat4 OG structure: uses sctransform::get_residual_var()
# on raw counts for both res_var_subset and bootstrap.
# Requires: Seurat, SeuratObject, sctransform, Matrix, tibble, dplyr,
#           foreach, doParallel, doRNG
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
    "feature"             = "gene",
    "log_feature"         = "log_gene",
    "nfeature"            = "ngene",
    "umi_per_feature"     = "umi_per_gene",
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
    cell_attr         = cell_attr,
    gene_attr         = gene_attr,
    model_pars_fit    = model_pars_fit,
    model_pars        = model_pars,
    model_pars_nonreg = NULL,
    model_str         = model_str,
    model_str_nonreg  = NULL,
    arguments         = if (is.null(arguments)) list() else arguments
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


# =============================================================================
# res_var_subset
# Calculate residual variance and geometric mean per gene for a subset of cells
# using the SCT model already fit on the full dataset.
# Uses sctransform::get_residual_var() on raw counts (same as Seurat4 OG).
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

  # --- raw counts for subset (Seurat5 layer access) ---
  umi.subset <- SeuratObject::LayerData(so[[count_slot_name]], layer = "counts")
  umi.subset <- umi.subset[, sel_cells, drop = FALSE]

  # --- reconstruct vst-like object from Seurat5 SCT model ---
  vst_out.subset <- .make_vst_like_from_seurat5(so, assay = sct_assay)

  # subset cell_attr to selected cells
  missing_cells <- setdiff(sel_cells, rownames(vst_out.subset$cell_attr))
  if (length(missing_cells) > 0L) {
    stop("Some cells missing in SCT cell.attributes: ",
         paste(head(missing_cells, 5), collapse = ", "))
  }
  vst_out.subset$cell_attr <- vst_out.subset$cell_attr[sel_cells, , drop = FALSE]

  # --- residual variance via sctransform (same as Seurat4 OG) ---
  getResVar.subset <- sctransform::get_residual_var(
    vst_out        = vst_out.subset,
    umi            = umi.subset,
    res_clip_range = c(-sqrt(ncol(umi.subset)), sqrt(ncol(umi.subset)))
  )

  # --- geometric mean: mean of log1p-normalised data across subset ---
  data_layer <- SeuratObject::LayerData(so[[sct_assay]], layer = "data")
  gmean_sub  <- data_layer[, sel_cells, drop = FALSE]

  if (inherits(gmean_sub, "Matrix")) {
    getgMean.subset <- Matrix::rowMeans(gmean_sub)
  } else {
    getgMean.subset <- rowMeans(gmean_sub)
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

  res_var_tbl <-
    foreach::foreach(
      cl = levels(cl_vec),
      .combine = "rbind",
      .export  = ls(globalenv())
    ) %do% {

      cells_filter <- so@meta.data[[meta_col]] == cl
      print(paste("cluster: ", cl))

      n_cells <- sum(cells_filter)
      print(paste("number of cells: ", n_cells))

      if (n_cells < 2) {
        message("WARNING: Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
        return(NULL)
      }

      res_var_mean_subs <- res_var_subset(
        so              = so,
        cells_filter    = cells_filter,
        count_slot_name = count_slot_name,
        sct_assay       = sct_assay
      )

      # convert to tibble and add cluster info
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
    dplyr::mutate(hvg = ifelse(res_var > hvg_cutoff, TRUE, FALSE))

  return(res_var_tbl)
}


# =============================================================================
# res_var_bootstrap
# Bootstrap residual variance by resampling cells WITH REPLACEMENT within each
# cluster, then calling sctransform::get_residual_var() on the resampled raw
# counts each time (faithful to the Seurat4 OG approach).
# Structure: sequential over clusters, parallel over iterations.
# =============================================================================
res_var_bootstrap <- function(so, clusters_slot, n_cycles = 1000,
                              n_cores = 10, hvg_cutoff = 5,
                              res_var_cl, keep_all = FALSE,
                              count_slot_name = "originalexp",
                              sct_assay = "SCT",
                              seed = 42) {

  # --- build the vst-like object once (replaces so@assays$SCT@misc$vst.out) ---
  vst_base <- .make_vst_like_from_seurat5(so, assay = sct_assay)

  # --- cluster vector ---
  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # --- parallel setup (across iterations, same as Seurat4 OG) ---
  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    `%op%` <- foreach::`%dopar%`
  } else {
    `%op%` <- foreach::`%do%`
  }
  set.seed(seed); doRNG::registerDoRNG(seed)

  # helper: restrict cell_attr to sampled cells
  .subset_vst <- function(vst_obj, sel_cells) {
    vst_local <- vst_obj
    vst_local$cell_attr <- vst_local$cell_attr[sel_cells, , drop = FALSE]
    vst_local
  }

  # --- sequential loop over clusters (same as Seurat4 OG) ---
  out_list <- vector("list", length(cl_levels))
  names(out_list) <- cl_levels

  for (cl in cl_levels) {

    message("Cluster ", cl)

    # cells in this cluster
    cells_filter <- cl_vec == cl
    cell_names   <- colnames(so)[cells_filter]
    n_cells      <- length(cell_names)

    if (n_cells < 2L) {
      message("  Skip: only ", n_cells, " cell(s)")
      next
    }

    # pre-slice counts for this cluster once (Seurat5 layer access)
    umi_cl <- SeuratObject::LayerData(so[[count_slot_name]], layer = "counts")
    umi_cl <- umi_cl[, cell_names, drop = FALSE]

    # hvgs for this cluster (for filtering output)
    hvgs <- res_var_cl %>%
      dplyr::filter(cluster == cl & hvg %in% TRUE)
    keep_genes <- if (keep_all) rownames(umi_cl) else hvgs$gene

    # parallel over iterations (reproducible via doRNG, same as Seurat4 OG)
    iter_tbl <-
      foreach::foreach(
        iter = seq_len(n_cycles),
        .combine  = "rbind",
        .packages = c("tibble", "dplyr", "sctransform"),
        .export   = c()
      ) %op% {

        # bootstrap sample: integer indices (faster than character indexing)
        sel <- sample.int(n_cells, size = n_cells, replace = TRUE)
        sel_cells <- cell_names[sel]

        # slice pre-sliced counts by integer position
        umi_subset     <- umi_cl[, sel, drop = FALSE]
        vst_out_subset <- .subset_vst(vst_base, sel_cells)

        # residual variance per gene for this bootstrap
        ResVar <- sctransform::get_residual_var(
          vst_out = vst_out_subset,
          umi     = umi_subset
        )

        boolean_hvg <- ResVar >= hvg_cutoff

        tb <- tibble::tibble(
          gene        = names(ResVar),
          boolean_hvg = unname(boolean_hvg),
          ResVar      = unname(ResVar),
          iter        = iter
        )

        # filter inside the loop to keep memory small (same as Seurat4 OG)
        if (!keep_all) {
          tb <- dplyr::filter(tb, gene %in% keep_genes)
        }

        tb
      }

    if (!is.null(iter_tbl) && nrow(iter_tbl)) {
      iter_tbl$cluster <- cl
      out_list[[cl]] <- iter_tbl
    }

    # free per-cluster memory (same as Seurat4 OG)
    rm(umi_cl); gc()
  }

  # --- combine results ---
  out <- dplyr::bind_rows(out_list)

  if (is.null(out) || nrow(out) == 0L) {
    return(tibble::tibble(
      gene        = character(),
      boolean_hvg = logical(),
      ResVar      = numeric(),
      iter        = integer(),
      cluster     = character()
    ))
  }

  return(out)
}


res_var_bootstrap_slow <- function(so, clusters_slot, n_cycles = 1000, n_cores, hvg_cutoff = 5, res_var_cl, keep_all = FALSE){
  #sct_assay <- "SCT"
 # clusters_slot <- "seurat_clusters"
 # n_cores <- 1
 # n_cycles <- 1


  # pull residuals and expression
  res_mat <- GetAssayData(so, assay = sct_assay, layer = "scale.data")
  if (is.null(res_mat) || nrow(res_mat) == 0L) stop("No residuals in SCT/scale.data.")
  gmn_mat <- GetAssayData(so, assay = sct_assay, layer = "data")
  if (is.null(gmn_mat) || nrow(gmn_mat) == 0L) stop("No data layer in SCT/data.")

  # align genes 
  genes <- intersect(rownames(res_mat), rownames(gmn_mat))
  res_mat <- res_mat[genes, , drop = FALSE]
  gmn_mat <- gmn_mat[genes, , drop = FALSE]

  # get clusters
  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # list for result
  out_list <- vector("list", length(cl_levels))
  names(out_list) <- cl_levels

  for (cl in cl_levels) {
    idx_all <- which(cl_vec == cl)
    n_cells <- length(idx_all)
    message("Processing cluster ", cl, " (", n_cells, " cell(s)).")
    if (n_cells < 2L) {
      message("Skip cluster ", cl, " (", n_cells, " cell(s)).")
      next
    }

    # slice matrices for cluster
    res_cl <- res_mat[, idx_all, drop = FALSE]
    gmn_cl <- gmn_mat[, idx_all, drop = FALSE]

    # parallelizing 
    if (n_cores > 1L) {
      doParallel::registerDoParallel(n_cores)
      `%op%` <- `%dopar%`
    } else {
      `%op%` <- `%do%`
    }
    set.seed(42); doRNG::registerDoRNG(42)

    # Bootstrap iterations for this cluster
    tbl_cl <-
      foreach(iter = seq_len(n_cycles), .combine = "rbind",
              .packages = c("matrixStats","tibble")) %op% {
        # sample same number of cells with replacement
        idx <- sample.int(n_cells, size = n_cells, replace = TRUE)

        rv <- rowVars(res_cl[, idx, drop = FALSE], na.rm = TRUE)
        gm <- rowMeans(gmn_cl[, idx, drop = FALSE])

        tibble(
          gene        = genes,
          res_var     = unname(rv),
          gmean       = unname(gm),
          cluster     = cl,
          boolean_hvg = rv > hvg_cutoff,
          iter        = iter
        )
      }

    out_list[[cl]] <- tbl_cl
    gc(FALSE)
  }

  # bind results
  out <- dplyr::bind_rows(Filter(Negate(is.null), out_list))
  if (nrow(out) == 0L) {
    tibble(
      gene = character(), res_var = numeric(), gmean = numeric(),
      cluster = character(), boolean_hvg = logical(), iter = integer()
    )
  } else {
    out
  }

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