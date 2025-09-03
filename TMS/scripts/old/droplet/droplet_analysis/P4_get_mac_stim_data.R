# xue et al 2014 response genes to certain cytokines
# here filtered before reading it in P5 load gene sets


library(readxl)
library(dplyr)
library(tidyr)

df <- read_excel("data/Xue_2014.xlsx", skip = 2, sheet = 'Table S1D') 

condition_cols <- setdiff(colnames(df), c("Illumina Probe ID", "Gene.Symbol", "Baseline"))

df_l2fc <- df %>%
  mutate(across(all_of(condition_cols), ~ . - Baseline, .names = "L2FC_{.col}"))

# significance threshold
threshold <- 1 

df_significant <- df_l2fc %>%
  pivot_longer(cols = starts_with("L2FC_"), names_to = "Condition", values_to = "LFC") %>%
  mutate(Condition = sub("L2FC_", "", Condition)) %>%
  filter(abs(LFC) > threshold) %>%
  group_by(Condition) %>%
  summarise(
    Upregulated = list(Gene.Symbol[LFC > threshold]),
    Downregulated = list(Gene.Symbol[LFC < -threshold])
  ) %>%
  ungroup()

df_significant$Condition
df_selection <- df_significant  %>% filter(Condition %in% c("IL10", "IL4", "TNF", "IFNg"))

check <- df_selection %>% select(Condition, Upregulated) %>% unnest(Upregulated)  %>%
     group_by(Condition) %>% summarise(count = n())
check


saveRDS(df_selection, "data/xue_response_genes.rds")

