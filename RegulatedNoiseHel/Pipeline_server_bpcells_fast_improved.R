# Masterscript pipeline server version
# adaptation 1: seurat 5 instead of seurat 4

Sys.setenv(RETICULATE_PYTHON = "/home/hkaufm49/anaconda3/envs/bpcells_env/bin/python")
library(reticulate)
library(tidyverse)
library(Seurat)
library(patchwork)
library(pryr) # RAM monitoring
library(BPCells)



#source('src/config_file.R')
source('src/FunctionLib_bpcells_fast_improved.R')
source('src/FunctionLibHel.R')

# importend from config file

SelectedTissues <- "Spleen"

impose_clusters <- TRUE
testrun <- FALSE # change bootstrap
n_cells <- 200 # only used if testrun == TRUE

min_cells <- 100 # minimum number of cells per subgroup/condition (in case of TMS age/sex group)
min_n_genes <- 500
min_n_counts <- 2500
max_n_counts <- 40000
max_percent_mt <- 5
npc_pca <- 100 #acutally 100
dim_umap <- 30 # acutally 30
dim_neighbors <- 30 # actually 30
cluster_res <- 1.3 # clustering resolution
num.de_cutoff <- 10 # number of violating genes in doublet detection
quantile_perc <- c(.025, .975) # bootstrapping, CI for ResVar
cycles_bootstrap <- 1000 # 1000
hvg_cutoff <- 5



# time monitoring
start <- Sys.time()
ram_start <- mem_used()


# SCTransform_v2 is an adaptation of the original Seurat wrapper, fixing the issue 
# that vst.out is not saved in the misc slot of the Seurat output object
#source('src/SCTransform_v2_seurat5.R')
# bind functions of this external version of SCT to Seurat
#e1 <- loadNamespace("Seurat")
#environment(SCTransform_v2) <- e1


# raw data folder
so_folder <- file.path("data/rawdata")

# Setup Output folders
out_folder <- "data/processed/added_bpcells_fast_improved"          
if (!file.exists(out_folder)) dir.create(out_folder)


fig_out_folder <- file.path(out_folder, "figures")         
if (!file.exists(fig_out_folder)) dir.create(fig_out_folder)


fig_UMAP <- file.path(fig_out_folder, "fig_UMAP")         
if (!file.exists(fig_UMAP)) dir.create(fig_UMAP)


res_doublets <- file.path(fig_out_folder, "res_doublets")         
if (!file.exists(res_doublets)) dir.create(res_doublets)


boot_folder <- paste0(out_folder, "/bootstrapping_tables")
if (!file.exists(boot_folder)) dir.create(boot_folder)

resvar_folder <- paste0(out_folder, "/res_var_tables")
if (!file.exists(resvar_folder)) dir.create(resvar_folder)

log_path <- file.path(out_folder, "logs")
if (!file.exists(log_path)) dir.create(log_path)   

object_path <- file.path(out_folder, "seurat_objects")
if (!file.exists(object_path)) dir.create(object_path) 


list_res <- list()

# Tissues <- c()
frac_UnstableCells <- c()

ResVar <- NULL
nCells_Cl <- NULL
gmeanClusters <- NULL
fraction_removed <- tibble(
  tissue = character(),
  age_group = character(),
  fractionUnstable = character())



#for (tissue in SelectedTissues) { # generalize tissue to input file (with different possible input file formats?)
  
  tissue <- "Spleen"
  gc()
  
  # 1.0 Select Seurat Obj
  #so.path <- file.path(so_folder, paste0(tissue, "_rawdata.rds"))
  so.path <- "data/rawdata/Spleen_rawdata.rds"
  so_full <- readRDS(so.path)


  # 0. Resaving matrix with Bpcells to optimize memory usage for large data sets (H)

  output_dir <- "data/"
  so_full <- resave_with_bp_cells(so_full, output_dir)




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
  #print(AgeSexGroups)


  for (age_group in AgeSexGroups) {

    #age_group <- "18m_female"

    so <- subset(so_full, subset = age_sex == age_group) # select relevant age group

    print(paste("processing", tissue, age_group, sep = " "))

    if (testrun) {
      so <- so[,1:n_cells]
    }
  

    # running seurat SCT (other verison does not support bpcell matrix)
    class(so@assays$originalexp$counts)
    so <- SCTransform(so, assay = "originalexp",  
                          verbose = FALSE)
  


    # 1.2 Run Dimensionality Reduction 
    so <- RunPCA(so, verbose = FALSE, npcs = npc_pca)
    so <- RunUMAP(so, dims = 1:dim_umap, verbose = FALSE)
    

    # 2.0 Cluster Identification
    print("clustering ...")

    so <- FindNeighbors(so, reduction = "pca", dims = 1:dim_neighbors, verbose = FALSE) %>%
      FindClusters(resolution = cluster_res, algorithm = 4, verbose = FALSE)
    # outlier pruning requires clusters to be stored in 'seurat_clusters' metadata slot (maybe change that?)
    #so$leiden_sct <- so$seurat_clusters
    #so$seurat_clusters <- NULL




    # OPTIONAL: Impose original clusters (only for reproduction)
    if(impose_clusters){

      print("Imposing original clusters ...")

      so <- impose_og_clusters(so, age_group) 

       so@meta.data[['seurat_clusters']] <- factor(so@meta.data[['seurat_clusters']])
        
        # test
        dim(so@assays$SCT@misc$vst.out$cell_attr)
        common_cells <- intersect(colnames(so),rownames(so@assays$SCT@misc$vst.out$cell_attr))

        so@assays$SCT@misc$vst.out$cell_attr <-  so@assays$SCT@misc$vst.out$cell_attr[common_cells, , drop = FALSE]



        ##DefaultAssay(so) <- "originalexp"
        ##so[["SCT"]] <- NULL



        #rerun SCT
        # test, the two # mean it was in the og version
        ##counts <- GetAssayData(so, assay = "originalexp", slot = "counts")

        ## vst_out <- sctransform::vst(
          ##umi           = counts,
         ## method        = "glmGamPoi", 
          ##residual_type = "pearson",
          ##verbosity     = 1, 
          ##return_cell_attr = TRUE)

        
        ##so[["SCT"]] <- Seurat::CreateAssayObject(data = vst_out$y)
        ##DefaultAssay(so) <- "SCT"
        ##Idents(so) <- "seurat_clusters"

        # save in misc for later
        ##so@assays$SCT@misc$vst.out <- vst_out

        # get residuals and store in scale.data
        ##resid_mat <- sctransform::get_residuals(
        ##  vst_out        = vst_out,
         ## umi            = counts)
          
        ##so <- SetAssayData(so, assay = "SCT", slot = "scale.data", new.data = resid_mat)

        # get resvar and store as variable features
        ##res_var  <- sctransform::get_residual_var(
         ## vst_out = vst_out, umi = counts)
        ##VariableFeatures(so) <- names(sort(res_var, decreasing = TRUE))[1:3000]

    }


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
    # H: becomes NULL
    so@assays$SCT@misc$vst.out$cell_attr <- so@assays$SCT@misc$vst.out$cell_attr[filter_mask,]
    
    #saveRDS(so, file.path(sceval_folder, tissue, 'full_sample_preprocess', paste0('subset_', SelectedFile)))


    
    # 3.0 Visualization of 1.2 & 2.0, UMAPs
    p1 <- UMAPPlot(so, group.by = "seurat_clusters") 
    
    p2 <- UMAPPlot(so, group.by = "cell_ontology_class")


    UMAPPatchwork <- (p1 | p2) +   plot_annotation(title = str_replace_all(tissue, '_', ' '), 
                                  subtitle = 'Clusters vs Cell Ontology Annotation', 
                                  tag_levels = c('A', '1'), 
                                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
                                  plot.subtitle = element_text(hjust = 0.5))) + plot_layout(ncol = 2) 
    
    ggsave(paste0(fig_UMAP, "/ResClustersCellStability_TMS_", tissue, "_", age_group, ".pdf"), UMAPPatchwork, width = 12, height = 6)
    
  
    

    # 5.0 Celltype Annotation
    # to be added
    
 
    
    # 6.0 HVG Identification with Stability Estimation (Bootstrapping)
    print("HVG detection ...")
    
    
    # 6.1 Get cluster specific residual variances
    res_var_cl <- res_var_by_group(so, 
                                   meta_col = "seurat_clusters", 
                                   hvg_cutoff = hvg_cutoff,
                                   count_slot_name = "originalexp")

    print("Bootstrapping ...")
    # 6.2 Bootstrapping
    bootstrap_res <- res_var_bootstrap(so = so, 
                                       clusters_slot = "seurat_clusters", 
                                       n_cycles = cycles_bootstrap, 
                                       n_cores = 10,
                                       hvg_cutoff = hvg_cutoff, 
                                       res_var_cl = res_var_cl,
                                       keep_all = FALSE)
    
    
    # 6.2.1 Plotting bootstrapping results (proposal)
    # resvar_plot <- plot_resvar(test_gene = "Actb", housek_gene = "Cxcl14", plot_cluster = "2")
    # 
    # res_var_cl_small <- res_var_cl %>% select(gene, cluster, hvg)
    # 
    # bs_hvg_cl <- bootstrap_res %>% 
    #   group_by(gene, cluster) %>% 
    #   summarize(mean_resvar_bs = mean(ResVar),
    #             perc.hvg = sum(boolean_hvg) / n(),
    #             quant_low = quantile(ResVar, probs = quantile_perc[1]), 
    #             quant_high = quantile(ResVar, probs = quantile_perc[2])) %>%
    #   left_join(res_var_cl_small, by = c("gene", "cluster"))
    # 
    # hist <- ggplot(subset(bs_hvg_cl, cluster == "2" & hvg), aes(x = perc.hvg)) +
    #   geom_histogram(bins = 20, color = "darkgray", fill = "#5b8dde") +
    #   ylab("# genes") +
    #   geom_vline(xintercept = 0.9, linetype = "dashed", 
    #              color = "red", size = 1) +
    #   theme_bw()
    # 
    # 
    # ggsave(file.path(proposal, "resvar_plot.pdf"), resvar_plot, device = "pdf", height = 2, width = 6)
    # ggsave(file.path(proposal, "hist.pdf"), hist, device = "pdf", height = 3, width = 4)
    
    
    
    saveRDS(res_var_cl, file = paste0(resvar_folder, "/res_var_cl_", tissue, "_", age_group, ".rds"))
    saveRDS(bootstrap_res, file = paste0(boot_folder,"/bootstrap_res", tissue, "_", age_group, ".rds"))
    

    # 6.3 Results
    stats <- bootstrap_res %>% 
      group_by(gene, cluster) %>% 
      summarize(mean_resvar_bs = mean(ResVar),
                perc.hvg = sum(boolean_hvg) / n(),
                quant_low = quantile(ResVar, probs = quantile_perc[1]), 
                quant_high = quantile(ResVar, probs = quantile_perc[2]))
    
    results_tmp_tibble <- left_join(res_var_cl, stats, by = c("gene" = "gene", "cluster" = "cluster")) %>%
      mutate(tissue = tissue, age = age_group)
    
    # adapt according to data structure
    # if (tissue == SelectedTissues[1]) {
    #   results_tibble <- results_tmp_tibble
    # } else {
    #   results_tibble <- rbind(results_tibble, results_tmp_tibble)
    # }
    
    if (age_group == AgeSexGroups[1]) {
      results_tibble <- results_tmp_tibble
    } else {
      results_tibble <- rbind(results_tibble, results_tmp_tibble)
    }
    
    
    # 7.0 save processed seurat object
    so_save.path <- file.path(object_path, paste0(tissue, "_", age_group, "_", "processed.rds"))
    saveRDS(so, file = so_save.path)
    
  }
  
  write_tsv(fraction_removed, file = file.path(out_folder, "fraction_unstable.tsv"))
  
#}


#saveRDS(results_tibble, file = file.path(out_folder, "results_tibble.rds"))
#write_tsv(results_tibble, file = file.path(out_folder, "results_tibble.tsv"))


  end <- Sys.time()
  ram_end <- mem_used()
  duration <- end - start

# time info
  log_file <- file.path(log_path, "pipeline_log.txt")
  writeLines(
  c(
    paste("Start:", start),
    paste("End:", end),
    paste("Duration (seconds):", as.numeric(duration, units = "secs")),
    paste("RAM at start:", round(as.numeric(ram_start) / 1024^3, 3), "GB"),
    paste("RAM at end:  ", round(as.numeric(ram_end) / 1024^3, 3), "GB"),
    paste("RAM diff:    ", round(as.numeric(ram_end - ram_start) / 1024^3, 3), "GB")),
    con = log_file)

