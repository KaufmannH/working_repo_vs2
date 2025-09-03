
library(readr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tibble)
library(tibble)
library(DElegate)
library(edgeR)
library(future)

# Part I: differential gene expression and selection of marker genes per condition
# Part II:  calculate module score of each condition specific genes in control cells


data <- readRDS("data/RDS_objects/ref_data_Macrophage.RDS")
data
head(data@meta.data)


# get marker genes per condition
Idents(data) <- "sample"

conditions_list <- unique(data@meta.data$sample)
testing_list <- list()

for (condition in conditions_list){
    row <- list(c("PBS", condition))
    names(row) <- condition 
   testing_list <- append(testing_list, row)
}
testing_list <- testing_list[order(names(testing_list))]

results_list <- list()

plan("multiprocess", workers = 4)
for (test in testing_list){ 
  print(test)
  condition_name <- test[2]
  de_results <- DElegate::findDE(object = data, compare = test, replicate_column = 'rep')
  de_results_signif <- de_results |> 
    filter(pvalue < 0.05) |>
    select('feature', 'log_fc', 'pvalue')
  results_list[condition_name] <- de_results_signif
}
plan("sequential")


# number of degs
number_of_degs_df <- data.frame(condition = names(results_list),  
                                number_degs = lengths(results_list)) |>
  arrange(desc(condition))
number_of_degs_df$condition <- factor(number_of_degs_df$condition, levels = number_of_degs_df$condition)


# how many genes were DEGs in each condition 
plots <- ggplot(number_of_degs_df, aes(x = condition , y = number_degs, fill = condition)) +
  geom_col(width = 0.5) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") + 
  labs(
    x = "Condition",          
    y = "Number of DEGs"      )
  

pdf("plots/number_of_degs.pdf", width = 5, height = 10)
plot(plots)
dev.off()

# save
saveRDS(results_list, "data/marker_genes_condition.rds")





# show module score

data <- readRDS("data/RDS_objects/ref_data_Macrophage.RDS")
marker_genes_list <- readRDS('data/marker_genes_condition.rds')
head(marker_genes_list)


for (condition in names(marker_genes_list)) {
  print(condition)
  genes_per_condition <- list(c(marker_genes_list[condition]))
  data <- AddModuleScore(object = data, features = genes_per_condition, ctrl = 100, name = condition)
  score <- colnames(data@meta.data) == paste0(condition, 1)
  colnames(data@meta.data)[score] <- paste0("score_", condition)
}
head(data@meta.data)

#violin plots for control cells with gene sets
extracted_data <- data@meta.data |>
  select(grep("^score_", colnames(data@meta.data), value = TRUE)) |>
  pivot_longer(cols = -1,
                names_to = "condition",
                values_to = "score") |>
  mutate(condition = sub("^[^_]+_", "", condition)) |>
  arrange(desc(condition))
extracted_data <- extracted_data |> mutate(condition = factor(condition, levels = unique(condition)))


plots <- ggplot(extracted_data, aes(x = condition , y = score, fill = condition)) +
  geom_violin(scale = "width" ) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 0.005) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
   labs(
    x = "Condition",          
    y = "Module score" )
  
  

pdf("plots/module_score.pdf", width = 5, height = 15)
plot(plots)
dev.off()





