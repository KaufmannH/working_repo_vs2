

gene <- "Lst1"
cl   <- "Lung_3m_male_17"

# TESTING THE RESVAR 0 THING


source('RegulatedNoiseHel/src/config_file.R')
source('RegulatedNoiseHel/src/FunctionLib.R')
source('RegulatedNoiseHel/src/FunctionLibHel.R')


# importend from config file

SelectedTissues <- "Spleen"

impose_clusters <- TRUE
testrun <- FALSE
n_cells <- 200 # only used if testrun == TRUE

min_cells <- 100 # minimum number of cells per subgroup/condition (in case of TMS age/sex group)
min_n_genes <- 500
min_n_counts <- 2500
max_n_counts <- 40000
max_percent_mt <- 5
npc_pca <- 100
dim_umap <- 30
dim_neighbors <- 30
cluster_res <- 1.3 # clustering resolution
num.de_cutoff <- 10 # number of violating genes in doublet detection
quantile_perc <- c(.025, .975) # bootstrapping, CI for ResVar
cycles_bootstrap <- 1000 # 1000
hvg_cutoff <- 5




# SCTransform_v2 is an adaptation of the original Seurat wrapper, fixing the issue 
# that vst.out is not saved in the misc slot of the Seurat output object
source('RegulatedNoiseHel/src/SCTransform_v2.R')
# bind functions of this external version of SCT to Seurat
e1 <- loadNamespace("Seurat")
environment(SCTransform_v2) <- e1


# raw data folder
so_folder <- file.path("RegulatedNoiseHel/data/rawdata")


  # 1.1 Filtering, QC, Normalization
  print("filtering, qc and normalization...")
  
  so_full[["percent.mt"]] <- PercentageFeatureSet(so_full, pattern = "^Mt")
  so_full[["age_sex"]] <- so_full@meta.data %>% mutate(age_sex = paste0(age, "_", sex)) %>% pull(age_sex)
  so_full <- subset(so_full, subset = nFeature_originalexp > min_n_genes & nCount_originalexp > min_n_counts & nCount_originalexp < max_n_counts & percent.mt < max_percent_mt)
  
  # only analyze subsets with sufficient number of cells
  AgeSexGroups <- unique(so_full@meta.data$age_sex)
  keep <- vector(mode = "logical", length = length(AgeSexGroups))

  for (i in 1:length(AgeSexGroups)) {
    keep[i] <- so_full@meta.data %>% filter(age_sex == AgeSexGroups[i]) %>% nrow() >= min_cells
  }
  AgeSexGroups <- factor(AgeSexGroups[keep])
  print(AgeSexGroups)


    age_group <- '3m_male'
    tissue <- 'Spleen'
    so <- subset(so_full, subset = age_sex == age_group) # select relevant age group

    print(paste("processing", tissue, age_group, sep = " "))

  
    
    so <- SCTransform_v2(so, assay = "originalexp", verbose = FALSE)


    # 1.2 Run Dimensionality Reduction 
    so <- RunPCA(so, verbose = FALSE, npcs = npc_pca)
    so <- RunUMAP(so, dims = 1:dim_umap, verbose = FALSE)
      
    
    # 2.0 Cluster Identification
    print("clustering ...")

    so <- FindNeighbors(so, reduction = "pca", dims = 1:dim_neighbors, verbose = FALSE) %>%
      FindClusters(resolution = cluster_res, algorithm = 4, verbose = FALSE)
    


    # 2.1 Outlier pruning
    
    # 2.1.1 Assign Unstable Cell Flag
    so <- AssignUnstableCellFlag(so)
    frac_UnstableCells <- round(sum(so$UnstableCells) / nrow(so[[]]), 3)
  
    new_row <- c(tissue = tissue, age_group = age_group, fractionUnstable = frac_UnstableCells)
    fraction_removed <- bind_rows(fraction_removed, new_row)
    
    print(paste("Fraction of unstable cells being removed:", 
                frac_UnstableCells, 
                sep = " "))
  

    # 2.1.2 Subset the Seurat Object excluding the unstable cells
    filter_mask <- !so@meta.data$UnstableCells
    so <- subset(so, subset = UnstableCells == FALSE)
    
    # vst_out is not automatically subset, so it needs to be done manually
    so@assays$SCT@misc$vst.out$cell_attr <- so@assays$SCT@misc$vst.out$cell_attr[filter_mask,]
    
    #saveRDS(so, file.path(sceval_folder, tissue, 'full_sample_preprocess', paste0('subset_', SelectedFile)))


  
    


    # 6.0 HVG Identification with Stability Estimation (Bootstrapping)
    print("HVG detection and bootstrapping ...")
    
    
    # 6.1 Get cluster specific residual variances
    res_var_cl <- res_var_by_group(so, 
                                   meta_col = "seurat_clusters", 
                                   hvg_cutoff = hvg_cutoff,
                                   count_slot_name = "originalexp")
    
    # test
 print("checking res_var == 0 cases ...")
    
    zero_cases <- res_var_cl %>%
      filter(res_var == 0) %>%
      arrange(desc(gmean))
    
    cat("total res_var == 0 entries:", nrow(zero_cases), "\n")
    cat("clusters where it happens:\n"); print(table(zero_cases$cluster))
    cat("gmean summary among zeros:\n"); print(summary(zero_cases$gmean))
    



    if (nrow(zero_cases) > 0) {
      
      vst_out  <- so@assays$SCT@misc$vst.out
      mp_all   <- vst_out$model_pars_fit
      ca_all   <- vst_out$cell_attr
      
      # check each (gene, cluster) pair: is every cell at the upper clip?
      check_clip <- zero_cases %>%
        rowwise() %>%
        mutate(
          n_cells       = sum(so@meta.data$seurat_clusters == cluster),
          clip_hi       = sqrt(n_cells),
          frac_above_clip = {
            cf  <- so@meta.data$seurat_clusters == cluster
            y   <- as.numeric(so@assays$originalexp@counts[gene, cf])
            mp  <- mp_all[gene, ]
            lu  <- ca_all[cf, "log_umi"]
            mu  <- exp(mp["(Intercept)"] + mp["log_umi"] * lu)
            zu  <- (y - mu) / sqrt(mu + mu^2 / mp["theta"])
            mean(zu > sqrt(sum(cf)))
          }
        ) %>%
        ungroup()
      
      cat("\nfraction of zero-resvar cases where ALL cells are clipped (frac == 1):\n")
      cat("  ", round(mean(check_clip$frac_above_clip == 1) * 100, 1), "%\n")
      
      cat("\ndistribution of frac_above_clip:\n")
      print(summary(check_clip$frac_above_clip))
      
      cat("\ntop 15 zero-resvar cases with clipping diagnosis:\n")
      print(head(check_clip, 15))
      

      # save for inspection
      write_tsv(check_clip, 
                paste0("ImmuneNoise/scripts/delete_chlipcheck.tsv"))
    }




  df <- 'droplet/droplet_raw_data/full_results_10_22.tsv'






# TESTING THE RESVAR 0 THING


source('RegulatedNoiseHel/src/config_file.R')
source('RegulatedNoiseHel/src/FunctionLib.R')
source('RegulatedNoiseHel/src/FunctionLibHel.R')


# importend from config file

SelectedTissues <- "Spleen"

impose_clusters <- TRUE
testrun <- FALSE
n_cells <- 200 # only used if testrun == TRUE

min_cells <- 100 # minimum number of cells per subgroup/condition (in case of TMS age/sex group)
min_n_genes <- 500
min_n_counts <- 2500
max_n_counts <- 40000
max_percent_mt <- 5
npc_pca <- 100
dim_umap <- 30
dim_neighbors <- 30
cluster_res <- 1.3 # clustering resolution
num.de_cutoff <- 10 # number of violating genes in doublet detection
quantile_perc <- c(.025, .975) # bootstrapping, CI for ResVar
cycles_bootstrap <- 1000 # 1000
hvg_cutoff <- 5




# SCTransform_v2 is an adaptation of the original Seurat wrapper, fixing the issue 
# that vst.out is not saved in the misc slot of the Seurat output object
source('RegulatedNoiseHel/src/SCTransform_v2.R')
# bind functions of this external version of SCT to Seurat
e1 <- loadNamespace("Seurat")
environment(SCTransform_v2) <- e1


# raw data folder
so_folder <- file.path("RegulatedNoiseHel/data/rawdata")


  # 1.1 Filtering, QC, Normalization
  print("filtering, qc and normalization...")
  
  so_full[["percent.mt"]] <- PercentageFeatureSet(so_full, pattern = "^Mt")
  so_full[["age_sex"]] <- so_full@meta.data %>% mutate(age_sex = paste0(age, "_", sex)) %>% pull(age_sex)
  so_full <- subset(so_full, subset = nFeature_originalexp > min_n_genes & nCount_originalexp > min_n_counts & nCount_originalexp < max_n_counts & percent.mt < max_percent_mt)
  
  # only analyze subsets with sufficient number of cells
  AgeSexGroups <- unique(so_full@meta.data$age_sex)
  keep <- vector(mode = "logical", length = length(AgeSexGroups))

  for (i in 1:length(AgeSexGroups)) {
    keep[i] <- so_full@meta.data %>% filter(age_sex == AgeSexGroups[i]) %>% nrow() >= min_cells
  }
  AgeSexGroups <- factor(AgeSexGroups[keep])
  print(AgeSexGroups)


    age_group <- '3m_male'
    tissue <- 'Spleen'
    so <- subset(so_full, subset = age_sex == age_group) # select relevant age group

    print(paste("processing", tissue, age_group, sep = " "))

  
    
    so <- SCTransform_v2(so, assay = "originalexp", verbose = FALSE)


    # 1.2 Run Dimensionality Reduction 
    so <- RunPCA(so, verbose = FALSE, npcs = npc_pca)
    so <- RunUMAP(so, dims = 1:dim_umap, verbose = FALSE)
      
    
    # 2.0 Cluster Identification
    print("clustering ...")

    so <- FindNeighbors(so, reduction = "pca", dims = 1:dim_neighbors, verbose = FALSE) %>%
      FindClusters(resolution = cluster_res, algorithm = 4, verbose = FALSE)
    


    # 2.1 Outlier pruning
    
    # 2.1.1 Assign Unstable Cell Flag
    so <- AssignUnstableCellFlag(so)
    frac_UnstableCells <- round(sum(so$UnstableCells) / nrow(so[[]]), 3)
  
    new_row <- c(tissue = tissue, age_group = age_group, fractionUnstable = frac_UnstableCells)
    fraction_removed <- bind_rows(fraction_removed, new_row)
    
    print(paste("Fraction of unstable cells being removed:", 
                frac_UnstableCells, 
                sep = " "))
  

    # 2.1.2 Subset the Seurat Object excluding the unstable cells
    filter_mask <- !so@meta.data$UnstableCells
    so <- subset(so, subset = UnstableCells == FALSE)
    
    # vst_out is not automatically subset, so it needs to be done manually
    so@assays$SCT@misc$vst.out$cell_attr <- so@assays$SCT@misc$vst.out$cell_attr[filter_mask,]
    
    #saveRDS(so, file.path(sceval_folder, tissue, 'full_sample_preprocess', paste0('subset_', SelectedFile)))


  
    


    # 6.0 HVG Identification with Stability Estimation (Bootstrapping)
    print("HVG detection and bootstrapping ...")
    
    
    # 6.1 Get cluster specific residual variances
    res_var_cl <- res_var_by_group(so, 
                                   meta_col = "seurat_clusters", 
                                   hvg_cutoff = hvg_cutoff,
                                   count_slot_name = "originalexp")
    
    # test
 print("checking res_var == 0 cases ...")
    
    zero_cases <- res_var_cl %>%
      filter(res_var == 0) %>%
      arrange(desc(gmean))
    
    cat("total res_var == 0 entries:", nrow(zero_cases), "\n")
    cat("clusters where it happens:\n"); print(table(zero_cases$cluster))
    cat("gmean summary among zeros:\n"); print(summary(zero_cases$gmean))
    



    if (nrow(zero_cases) > 0) {
      
      vst_out  <- so@assays$SCT@misc$vst.out
      mp_all   <- vst_out$model_pars_fit
      ca_all   <- vst_out$cell_attr
      
      # check each (gene, cluster) pair: is every cell at the upper clip?
      check_clip <- zero_cases %>%
        rowwise() %>%
        mutate(
          n_cells       = sum(so@meta.data$seurat_clusters == cluster),
          clip_hi       = sqrt(n_cells),
          frac_above_clip = {
            cf  <- so@meta.data$seurat_clusters == cluster
            y   <- as.numeric(so@assays$originalexp@counts[gene, cf])
            mp  <- mp_all[gene, ]
            lu  <- ca_all[cf, "log_umi"]
            mu  <- exp(mp["(Intercept)"] + mp["log_umi"] * lu)
            zu  <- (y - mu) / sqrt(mu + mu^2 / mp["theta"])
            mean(zu > sqrt(sum(cf)))
          }
        ) %>%
        ungroup()
      
      cat("\nfraction of zero-resvar cases where ALL cells are clipped (frac == 1):\n")
      cat("  ", round(mean(check_clip$frac_above_clip == 1) * 100, 1), "%\n")
      
      cat("\ndistribution of frac_above_clip:\n")
      print(summary(check_clip$frac_above_clip))
      
      cat("\ntop 15 zero-resvar cases with clipping diagnosis:\n")
      print(head(check_clip, 15))
      

      # save for inspection
      write_tsv(check_clip, 
                paste0("ImmuneNoise/scripts/delete_chlipcheck.tsv"))
    }




  df <- read_tsv('ImmuneNoise/droplet/droplet_raw_data/full_results_10_22.tsv')

head(df)

zero_resvar <- df |>
  filter(res_var == 0) |>
  mutate(tissue_cluster = paste(tissue, age, cluster, sep = "_")) |>
  select(tissue_cluster, gene, gmean, res_var) |>
  arrange(gene)

zero_resvar

write_tsv(zero_resvar, 'res_var_0_table.tsv')

unique(zero_resvar$gene)












