# according to Thron et al.: testing for introduced biases: Thron et al. 2025
# Is there a lot of gene expression correlation in the data because of technical artefacts? 

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(tibble)



# add and filter main df from TMS
df_main_all <- read.csv("data/combined_data.csv") 
colnames(df_main_all)


# correlation between two ranked variable lists
df_ranked <- df_main_all %>%
  mutate(age = str_extract(age, "\\d+m")) %>%  
  mutate(sample = paste(age, tissue, sep="_")) %>%
  mutate(age = as.numeric(str_remove(age, "m"))) %>% 
  filter(age == 3) %>%
  group_by(cluster_id) %>%
  select(gene, gmean, cluster_id)
df_ranked

wide_df <- df_ranked %>%
  pivot_wider (
    values_from = gmean,
    names_from = cluster_id
  )
wide_df

s_wide_df <- wide_df %>%
  group_by(gene) %>%
  column_to_rownames(var = "gene") 
  #summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) # get rid of multiple rows per gene
head(s_wide_df)


original_correlation <- cor(s_wide_df, method = "spearman",  use = "pairwise.complete.obs")

# shuffle to get the baseline correlation
shuffled_df <- s_wide_df
rownames(shuffled_df) <- sample(rownames(shuffled_df))
shuffled_correlation <- cor(shuffled_df, method = "pearson", use = "pairwise.complete.obs")
shuffled_correlation


cor_df <- data.frame(
  Original = as.vector(original_correlation),
  Shuffled = as.vector(shuffled_correlation)
)

shuffled <- ggplot(cor_df, aes(x = Shuffled, y = Original)) +
  geom_point(alpha = 0.5, color = "#70916B") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Original vs. Shuffled Correlation",
       x = "Shuffled Correlation",
       y = "Original Correlation") +
  theme_classic()
ggsave('plots/3_m/corr_QC.png', plot = shuffled,  width = 7, height = 7 )


# points above the diagonal means that the biological correlation is bigger than the technical one, on the diagonal means that its quite equal, bad sign. 
# takes super long