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


# This function takes as input a seurat object and a logical vector for defining a subset of cells
# and calculates the residual variance and geometric mean for each gene across the subset of cells.
# It requires an SCT model and a populated vst.out slot in so[["SCT"]]@misc$vst.out
res_var_subset <- function(so, cells_filter, count_slot_name, vst_out = NULL, integrated = FALSE){
  
  # raw counts for subset
  umi.subset <- so@assays[[count_slot_name]]@counts[,cells_filter]

    # vst.out for subset
    if (is.null(vst_out)) {
      vst_out.subset           <- so@assays$SCT@misc$vst.out
      vst_out.subset$cell_attr <- vst_out.subset$cell_attr[cells_filter,]
    } else {
      vst_out.subset           <- vst_out
      vst_out.subset$cell_attr <- vst_out.subset$cell_attr[cells_filter,]
    }
    
    
    # recalculate residual variance for genes across the subset
    getResVar.subset  <- sctransform::get_residual_var(vst_out = vst_out.subset, umi = umi.subset,
                                                      res_clip_range = c(-sqrt(ncol(umi.subset)), 
                                                                          sqrt(ncol(umi.subset))))

  # calculate geometric mean (mean of log1p of corrected counts)
  if(!integrated) {
    getgMean.subset <- apply(so@assays$SCT@data[,cells_filter], 1, mean)
  } else {
    getgMean.subset <- apply(so@assays$integrated@data[,cells_filter], 1, mean)
  }
  
  
  # output matrix
  res_var_mean <- cbind(res_var = getResVar.subset, 
                        gmean = getgMean.subset[names(getResVar.subset)])
  
  return(res_var_mean)
}


# This function takes as input a seurat object and a metadata column. It will reestimate
# the variability and gmean of all genes across the cells in each group (defined by metadata column)
# and return results as a tibble. meta_col should be a factor
res_var_by_group <- function(so, meta_col, hvg_cutoff = 5, count_slot_name, n_cores, vst_out = NULL, integrated = FALSE){
  #lvg_cutoff = 1

  #cluster <- makeCluster(n_cores)
  #registerDoParallel(cluster)
  
  res_var_tbl <- 
    foreach(cl = levels(so@meta.data[[meta_col]]), .combine = 'rbind', .export = ls(globalenv())) %do% {
      cells_filter <- so@meta.data[[meta_col]] == cl
    
      #H
      n_cells <- sum(cells_filter)
      if (n_cells < 2) {
        message("WARNING: Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
        return(NULL)}

      res_var_mean_subs <- res_var_subset(so, cells_filter, count_slot_name, vst_out = vst_out, integrated = integrated)
      
      # convert to tibble and add cluster info
      res_var_tbl <- as_tibble(res_var_mean_subs, rownames = "gene") %>%
        mutate(cluster = cl)
        
    }
  
  #stopCluster(cluster)
  
  res_var_tbl <- res_var_tbl %>%
    mutate(
      hvg = res_var > hvg_cutoff,
     # lvg = res_var <= lvg_cutoff
    )

  return(res_var_tbl)
  
}

# This function takes as an input a seurat object, clusters and a number of cycles and will
# then perform that many cycles of bootstrapping by choosing random subsets of cells
# and re-calculating the variability for each gene in those subsets.
# Structure: sequential over clusters, parallel over iterations.
res_var_bootstrap <- function(
  so,
  clusters_slot,
  n_cycles         = 1000,
  n_cores          = 10,
  hvg_cutoff       = 5,
  res_var_cl,
  keep_all         = FALSE,
  count_slot_name  = "originalexp",
  vst_out          = NULL,
  seed             = 42
){


  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # parallel  
  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    `%op%` <- `%dopar%`
  } else {
    `%op%` <- `%do%`
  }
  set.seed(seed); doRNG::registerDoRNG(seed)

  # disable sctransforms internal parallelism 
  old_mc_cores <- getOption("mc.cores")
  old_future   <- future::plan()
  options(mc.cores = 1L)
  future::plan("sequential")
  on.exit({
    options(mc.cores = old_mc_cores)
    future::plan(old_future)
  }, add = TRUE)


  .subset_vst <- function(vst_base, sel_cells) {
    vst_local <- vst_base
    vst_local$cell_attr <- vst_local$cell_attr[sel_cells, , drop = FALSE]
    vst_local
  }

  if (is.null(vst_out)) {
    if (is.null(so@assays$SCT@misc$vst.out)) stop("No vst.out found at so@assays$SCT@misc$vst.out")
    vst_base <- so@assays$SCT@misc$vst.out
  } else {
    vst_base <- vst_out
  }

  # loop over clusters
  out_list <- vector("list", length(cl_levels))
  names(out_list) <- cl_levels

  for (cl in cl_levels) {
    message("Cluster ", cl)
    cells_filter <- cl_vec == cl
    cell_names   <- colnames(so)[cells_filter]
    n_cells      <- length(cell_names)

    if (n_cells < 2L) {
      message("  Skip: only ", n_cells, " cell(s)")
      next
    }

    umi_cl <- so@assays[[count_slot_name]]@counts[, cell_names, drop = FALSE]

    # pre-generate unique IDs (reused every iteration)
    unique_ids <- paste0("cell_", seq_len(n_cells))

    hvgs <- res_var_cl %>%
      dplyr::filter(cluster == cl & hvg %in% TRUE)
    keep_genes <- if (keep_all) rownames(umi_cl) else hvgs$gene

    # parallel over iterations
    iter_tbl <-
      foreach::foreach(iter = seq_len(n_cycles),
                       .combine = "rbind",
                       .packages = c("tibble", "dplyr", "sctransform"),
                       .export   = c()) %op% {

        # bootstrap sample of cells with replacement
        sel <- sample.int(n_cells, size = n_cells, replace = TRUE)
        sel_cells <- cell_names[sel]

        # slice counts
        umi_subset     <- umi_cl[, sel, drop = FALSE]
        vst_out_subset <- .subset_vst(vst_base, sel_cells)

        colnames(umi_subset) <- unique_ids
        rownames(vst_out_subset$cell_attr) <- unique_ids

        # residual variance per gene for this bootstrap
        ResVar <- sctransform::get_residual_var(
          vst_out = vst_out_subset,
          umi     = umi_subset)

        boolean_hvg <- ResVar >= hvg_cutoff

        tb <- tibble::tibble(
          gene        = names(ResVar),
          boolean_hvg = unname(boolean_hvg),
          ResVar      = unname(ResVar),
          iter        = iter)

        if (!keep_all) {
          tb <- dplyr::filter(tb, gene %in% keep_genes)
        }

        tb
      }

    if (!is.null(iter_tbl) && nrow(iter_tbl)) {
      iter_tbl$cluster <- cl
      out_list[[cl]] <- iter_tbl
    }

    rm(umi_cl); gc()
  }

  out <- dplyr::bind_rows(out_list)
  if (is.null(out) || nrow(out) == 0L) {
    tibble::tibble(gene = character(), boolean_hvg = logical(),
                   ResVar = numeric(), iter = integer(), cluster = character())
  } else {
    out
  }
}



# mar 9th 2026: crashed bc of nested loops in sctransform::get_residual_var(), fixed issue and made new function
res_var_bootstrap_og2 <- function(
  so,
  clusters_slot,
  n_cycles    = 1000,
  n_cores     = 10,
  hvg_cutoff  = 5,
  res_var_cl,
  keep_all    = FALSE,
  count_slot_name = "originalexp",  
  vst_out          = NULL,          
  seed             = 42
 # lvg_cutoff       = 1
){
#clusters_slot <- "seurat_clusters"
#n_cycles <- 1
#n_cores  <-  1
#res_var_cl <- "1"
#count_slot_name = "originalexp"
 #seed             = 42
 #keep_all = FALSE
# lvg_cutoff       = 1
 #vst_out <- so@assays$SCT@misc$vst.out



  # sanity: levels
  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # parallel backend + reproducible RNG
  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    `%op%` <- `%dopar%`
  } else {
    `%op%` <- `%do%`
  }
  set.seed(seed); doRNG::registerDoRNG(seed)

  # helper to fetch vst.out and restrict cell_attr to sampled cells
  .subset_vst <- function(vst_base, sel_cells) {
    vst_local <- vst_base
    vst_local$cell_attr <- vst_local$cell_attr[sel_cells, , drop = FALSE]
    vst_local
  }

  # base vst.out (from object or provided)
  if (is.null(vst_out)) {
    if (is.null(so@assays$SCT@misc$vst.out)) stop("No vst.out found at so@assays$SCT@misc$vst.out")
    vst_base <- so@assays$SCT@misc$vst.out
  } else {
    vst_base <- vst_out
  }

  # iterate clusters sequentially (each cluster fans iterations in parallel)
  out_list <- vector("list", length(cl_levels))
  names(out_list) <- cl_levels

  for (cl in cl_levels) {
    message("Cluster ", cl)
    cells_filter <- cl_vec == cl
    cell_names   <- colnames(so)[cells_filter]
    n_cells      <- length(cell_names)

    if (n_cells < 2L) {
      message("  Skip: only ", n_cells, " cell(s)")
      next
    }

    # counts for this cluster once
    umi_cl <- so@assays[[count_slot_name]]@counts[, cell_names, drop = FALSE]

    hvgs <- res_var_cl %>%
    dplyr::filter(cluster == cl & hvg %in% TRUE)
    #hvgs <- dplyr::filter(res_var_cl, cluster == cl, hvg)
    #lvgs <- dplyr::filter(res_var_cl, cluster == cl, lvg)
    keep_genes <- if (keep_all) rownames(umi_cl) else hvgs$gene #union(hvgs$gene, lvgs$gene)

    # parallel over iterations (reproducible)
    iter_tbl <-
      foreach::foreach(iter = seq_len(n_cycles),
                       .combine = "rbind",
                       .packages = c("tibble","dplyr","sctransform"),
                       .export   = c()) %op% {

        # bootstrap sample of cells with replacement
        sel <- sample.int(n_cells, size = n_cells, replace = TRUE)
        sel_cells <- cell_names[sel]

        # slice counts; set up matching vst.out (restrict cell_attr)
        umi_subset     <- umi_cl[, sel, drop = FALSE]
        vst_out_subset <- .subset_vst(vst_base, sel_cells)

        # residual variance per gene for this bootstrap
        ResVar <- sctransform::get_residual_var(
          vst_out = vst_out_subset,
          umi     = umi_subset)

        boolean_hvg <- ResVar >= hvg_cutoff
       # boolean_lvg <- ResVar <= lvg_cutoff  

        tb <- tibble::tibble(
          gene         = names(ResVar),
          boolean_hvg  = unname(boolean_hvg),
          #boolean_lvg  = unname(boolean_lvg),   
          ResVar       = unname(ResVar),
          iter         = iter)


        if (!keep_all) {
          tb <- dplyr::filter(tb, gene %in% keep_genes)
        }

        tb
      }

    if (!is.null(iter_tbl) && nrow(iter_tbl)) {
      iter_tbl$cluster <- cl
      out_list[[cl]] <- iter_tbl
    }
    # free memory per cluster
    rm(umi_cl); gc()
  }

  out <- dplyr::bind_rows(out_list)
  if (is.null(out)) {
    tibble::tibble(gene=character(), boolean_hvg=logical(), #boolean_lvg=logical(),
               ResVar=numeric(), iter=integer(), cluster=character())
  } else {
    out
  }
}




res_var_bootstrap_og <- function(so, clusters_slot, n_cycles = 1000, n_cores, hvg_cutoff = 5, res_var_cl, keep_all = FALSE){
  
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
      boolean_hvg <- features %in% names(ResVar[ResVar >= hvg_cutoff])
      
      cbind(boolean_hvg, ResVar, iter)
      
    }
    
    #stopCluster(cluster)
    
    # hvgs for cluster
    hvgs <- res_var_cl %>% filter(cluster == cl,
                                  hvg)
    
    # filter results for hvgs
    if (keep_all) {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene")
    } else {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene") %>%
        filter(gene %in% hvgs$gene)
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