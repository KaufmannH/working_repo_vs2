# What are the distributions of LVG/HVG and marker genes in all genes vs in innate immune genes?
# only in Lung, 3m, monocytes, macrophages, DCs


library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(purrr)


iros <- "👑"
ei <- "🥚"
stephannes <- "👨‍❤️‍👨"


# colors
palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")

colors <- c('#111D4E', '#4472C4', '#5FB4E5', '#3CBFAE', '#70AD47', '#EE8339', '#8E0C72')

# load tagged df
df_main_with_markers_expr_level <- readRDS("droplet_analysis/data/df_with_markers_2.rds") # NEW  DF (01.03.!!!!!!)

check <- df_main_with_markers_expr_level %>%
  group_by(gene_set) %>%
  summarise(count = n(), .groups = 'drop')
check
colnames(df_main_with_markers_expr_level)
unique(df_main_with_markers_expr_level$cell_type)

df_main_with_markers_expr_level <- df_main_with_markers_expr_level %>%
 mutate(cluster_id = recode(cluster_id,
    "Heart_3m_female_7"        = "Heart (female, cluster 7)",
    "Heart_3m_male_5"          = "Heart (male, cluster 5)",
    "Heart_3m_male_6"          = "Heart (male, cluster 6)",
    "Heart_3m_male_13"         = "Heart (male, cluster 13)",
    "Aorta_3m_male_6"          = "Aorta (male, cluster 6)",
    "Diaphragm_3m_male_7"      = "Diaphragm (male, cluster 7)",
    "MAT_3m_female_6"          = "Adipose (female, cluster 6)",
    "MAT_3m_female_9"          = "Adipose (female, cluster 9)",
    "MAT_3m_male_7"            = "Adipose (male, cluster 7)",
    "Brain_Myeloid_3m_male_8"  = "Brain (male, cluster 8)",
    "Kidney_3m_female_5"       = "Kidney (female, cluster 5)",
    "Kidney_3m_male_4"         = "Kidney (male, cluster 4)",
    "Limb_Muscle_3m_female_7"  = "LimbMuscle (female, cluster 7)",
    "Liver_3m_male_3"          = "Liver (male, cluster 3)",
    "Lung_3m_female_3"         = "Lung (female, cluster 3)",
    "Lung_3m_male_3"           = "Lung (male, cluster 3)",
    "Marrow_3m_male_14"        = "Marrow (male, cluster 14)",
    "Marrow_3m_female_10"      = "Marrow (female, cluster 10)",
    "Trachea_3m_male_8"        = "Trachea (male, cluster 8)",
    "Trachea_3m_male_11"       = "Trachea (male, cluster 11)",
    "Trachea_3m_female_8"      = "Trachea (female, cluster 8)",
    "Trachea_3m_female_10"     = "Trachea (female, cluster 10)"
  ))



# TEST 1
# version 1: get only lps response genes vs other genes. exclude non expressed genes. from non blown up df

df_plot_1 <- df_main_with_markers_expr_level %>%
  filter(!is.na(cell_type), !is.na(cluster_id)) %>%
  group_by(cluster_id, cell_type, gene_set, lps_stim) %>%
  summarise(count = n(), .groups = 'drop')
head(df_plot_1)

my_colors <- c("TRUE" = "#556967","FALSE" = "#5AB4AC" )
df_plot_1$gene_set <- factor(df_plot_1$gene_set, levels = c("HVG", "Other", "LVG","Response gene not expressed"))

barplot_relative <- ggplot(df_plot_1, aes(x = gene_set, y = count, fill = lps_stim)) +
  geom_bar(stat = "identity", position = "fill" ) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type",
       title = 'LPS stimulated genes') +

  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)
  )

ggsave('plots/3_m/Test_1/relative_abundances1.png', plot = barplot_relative)
# There are more LPS response genes in HVG category than in LVG. 




# version 2: distribution of all tags in the data set, with genes being in multiple categories

long_df_prep <- df_main_with_markers_expr_level %>%
   rename(
    `Immune Response Gene` = response_expressed,
    `Housekeeping Gene Lin` = housekeeping_gene_lin,
    `Housekeeping Gene` = housekeeping_gene, 
    `Response gene not expressed` = response_not_expressed, 
    `LPS response early` = lps_stim_early,
    `LPS response late` = lps_stim_late,
    `IFNy response` = ifng_gene,
    `Il10 response` = il10_gene,
    `Il4 response` = il4_gene,
    `TNF response` = tnf_gene,
    `Meiosis` = meiosis_gene,
    `Sperm DNA condensation` = sperm_gene) 

long_df <- long_df_prep %>%
  pivot_longer(
    cols = c('Immune Response Gene', "Housekeeping Gene Lin", "Housekeeping Gene", "Response gene not expressed", "LPS response early", "LPS response late", "IFNy response", "Il10 response" , "Il4 response", "TNF response", "Meiosis","Sperm DNA condensation", "Other"),
    names_to = "category",
    values_to = "values"
  ) 
head(long_df )#[1:10, 10:13]


cleaned_long_df <- long_df %>%
   filter(values == TRUE) 


category_counts <- long_df %>%
  filter(values == TRUE) %>%
  filter(!is.na(cell_type), !is.na(cluster_id)) %>%
  group_by(cluster_id,cell_type, gene_set, category) %>%
  summarise(count = n(), .groups = "drop") 
category_counts

#my_colors <- c("Marker Gene" = "#B8C7C7", "Other" = "#556967", "Immune Response Gene" = "#E0D5C8", "Stable Housekeeping Gene" = "#AE9982", "Housekeeping Gene" = "#5AB4AC", "Response gene not expressed" = "#A4C089")
category_counts$gene_set <- factor(category_counts$gene_set, levels = c("HVG", "Other", "LVG", "Response gene not expressed"))

barplot_relative <- ggplot(category_counts, aes(x = gene_set, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "fill" ) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +

  #scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)
  )
ggsave('plots/3_m/Test_1/relative_abundances2_1.png', plot = barplot_relative)

# marker genes are more in HVG, house keeping more in LVG. 



# verison 3: like version 2 but turned

category_counts <- bright_future_df %>%
  filter(values == TRUE) %>%
  filter(!is.na(cell_type), !is.na(cluster_id)) %>%
  group_by(cluster_id,cell_type, gene_set, category) %>%
  summarise(count = n(), .groups = "drop") 
category_counts

my_colors <- c("HVG" = "#A4C089", "Other"= "#556967", "LVG"= "#5AB4AC", "Response gene not expressed"= "#B8C7C7")
category_counts$category <- factor(category_counts$category, levels = c("Marker Gene" , "Other" , "Immune Response Gene" , "Stable Housekeeping Gene" , "Housekeeping Gene", "Response gene not expressed" ) )


barplot_relative <- ggplot(category_counts, aes(x = category, y = count, fill = gene_set)) +
  geom_bar(stat = "identity", position = "fill" ) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
ggsave('plots/3_m/Test_1/relative_abundances3.png', plot = barplot_relative, width = 8, height = 11)


# 5: cytokine response (old version)

plot_5_df <- df_main_with_markers %>%
  mutate(category = case_when(
    cytokine_gene & is_marker ~ "Cytokine response gene and Marker Gene",
    cytokine_gene ~ "Cytokine response gene",
    is_marker ~ "Marker Gene",
    TRUE ~ "Other"
  ))
head(plot_5_df)


p <- plot_5_df %>%
  group_by(category) %>%
  summary(count = n(), drop.groups() )
p

my_colors <- c("Marker Gene" = "#5AB4AC",  "Other" = "#556967", "Cytokine response gene" = "#A4C089", "Cytokine response gene and Marker Gene" = "#B8C7C7")
plot_5_df$category <- factor(plot_5_df$category, levels = c("Marker Gene", "Other", "Cytokine response gene", "Cytokine response gene and Marker Gene"))


 boxplot_expression <- ggplot(plot_5_df, aes(x = category, y = gmean, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.colour = "black", outlier.alpha = 0.1, size = 0.1) + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(y = "Expression Level") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5)
  )

ggsave('plots/3_m/Test_1/relative_abundances5.png', plot = boxplot_expression,  width = 9, height = 10)



# 6: TLR pathway (old version)
sum(df_main_with_markers$tlr_gene)


plot_6_df <- df_main_with_markers %>%
  mutate(category = case_when(
    tlr_gene & is_marker ~ "TLR pathway and Marker Gene",
    tlr_gene ~ "TLR pathway gene",
    is_marker ~ "Marker Gene",
    TRUE ~ "Other"
  ))
head(plot_6_df)


p <- plot_6_df %>%
  group_by(category) %>%
  summary(count = n(), drop.groups() )
p

my_colors <- c("Marker Gene" = "#5AB4AC",  "Other" = "#556967", "TLR pathway gene" = "#A4C089", "TLR pathway and Marker Gene" = "#B8C7C7")
plot_6_df$category <- factor(plot_6_df$category, levels = c("Marker Gene", "Other", "TLR pathway gene", "TLR pathway and Marker Gene"))


 boxplot_expression <- ggplot(plot_6_df, aes(x = category, y = gmean, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.colour = "black", outlier.alpha = 0.1, size = 0.1) + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(y = "Expression Level") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5)
  )

ggsave('plots/3_m/Test_1/relative_abundances6.png', plot = boxplot_expression,  width = 10, height = 12)



# version 7: like v1 but classes of lps genes (get only lps response genes vs other genes. exclude non expressed genes. from non blown up df.)

ear <- unique(subset(df_main_with_markers_expr_level, lps_stim_early == TRUE)$gene)
late <- unique(subset(df_main_with_markers_expr_level, lps_stim_late == TRUE)$gene)
inter <- intersect(ear, late)
inter
length(late)

df_plot_1 <- df_main_with_markers_expr_level %>%
  filter(!is.na(cell_type), !is.na(cluster_id)) %>%
  filter(cell_type == "Macrophage") %>%
  mutate(
    lps_stim = case_when(
      lps_stim_early == TRUE ~ "Early LPS response genes",
      lps_stim_late == TRUE ~ "Late LPS response genes",
      TRUE ~ "Not LPS activated")) %>%
  group_by(cluster_id, cell_type, gene_set, lps_stim) %>%
  summarise(count = n(), .groups = 'drop')
head(df_plot_1)


colors <- c("Early LPS response genes" =  '#7F7F7F', "Late LPS response genes" = '#EE8339', "Not LPS activated" = '#3CBFAE')
df_plot_1$gene_set <- factor(df_plot_1$gene_set, levels = c("HVG", "Other", "LVG","Response gene not expressed"))
df_plot_1$lps_stim <- factor(df_plot_1$lps_stim, levels = c( "Not LPS activated", "Early LPS response genes", "Late LPS response genes"))

barplot_relative <- ggplot(df_plot_1, aes(x = gene_set, y = count, fill = lps_stim)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Variability class",
       y = "Fraction of genes",
       fill = "Gene type",
       title = 'LPS stimulated genes') +
  scale_fill_manual(values = colors) +
  theme_classic() +
  coord_flip() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 9)
  )

ggsave('plots/3_m/Test_1/relative_abundances1_2.png', plot = barplot_relative, width = 10, height = 5)

# flipped 
colors <- c("Other" =  '#7F7F7F', "LVG" = '#EE8339', "HVG" = '#3CBFAE')
barplot_relative <- ggplot(df_plot_1, aes(x = lps_stim, y = count, fill = gene_set)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene type",
       y = "Fraction of genes",
       fill = "Variability class") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  coord_flip() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 9)
  )

ggsave('plots/3_m/Test_1/relative_abundances1_3.png', plot = barplot_relative, width = 10, height = 5)

# There are more LPS response genes in HVG category than in LVG. 


colnames(df_main_with_markers)




# TEST 2
# version 4: expression comparison of marker, lps, housekeeping 

plot_4_df <- df_main_with_markers %>%
  mutate(category = case_when(
    lps_stim & is_marker ~ "LPS response gene and Marker Gene",
    lps_stim ~ "LPS response gene",
    is_marker ~ "Marker Gene",
    housekeeping_gene ~ "Housekeeping Gene",
    housekeeping_gene & lps_stim ~ "LPS & Housekeeping Gene",
    housekeeping_gene & is_marker ~ "Marker & Housekeeping Gene",
    TRUE ~ "Other"
  ))
head(plot_4_df)
sum(plot_4_df$housekeeping_gene)


p <- plot_4_df %>%
  group_by(cell_type, cluster_id, category) %>%
  summarise(count = n(), .groups = 'drop')
p

#my_colors <- c("Marker Gene" = "#5AB4AC",  "Other" = "#556967", "LPS response gene" = "#A4C089", "LPS response gene and Marker Gene" = "#B8C7C7")
plot_4_df$category <- factor(plot_4_df$category, levels = c("Marker Gene", "Other", "LPS response gene", "LPS response gene and Marker Gene", "Housekeeping Gene","LPS & Housekeeping Gene", "Marker & Housekeeping Gene"))


 boxplot_expression <- ggplot(plot_4_df, aes(x = category, y = gmean, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.colour = "black", outlier.alpha = 0.1, size = 0.1) + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(y = "Expression Level") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5)
  )
ggsave('plots/3_m/Test_2/boxplot4.png', plot = boxplot_expression,  width = 8, height = 8)


#version 4.1: boxplot but with longer format = non exclusive genes

cleaned_long_df <- long_df %>%
   filter(values == TRUE) 


#my_colors <- c("Marker Gene" = "#B8C7C7", "Other" = "#556967", "Immune Response Gene" = "#E0D5C8", "Stable Housekeeping Gene" = "#AE9982", "Housekeeping Gene" = "#5AB4AC", "Response gene not expressed" = "#A4C089")
#summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


 boxplot_expression <- ggplot(cleaned_long_df, aes(x = category, y = gmean, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.colour = "black", outlier.alpha = 0.1, size = 0.1) + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(y = "Expression Level") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5)
  )
ggsave('plots/3_m/Test_2/boxplot4_1.png', plot = boxplot_expression,  width = 8, height = 8)


# verison 4.2 longer format (non exclusive) but just lps

colors <- c('#111D4E', '#4472C4', '#5FB4E5', '#3CBFAE', '#70AD47', '#EE8339', '#8E0C72')

long_df <- long_df_prep %>%
  pivot_longer(
    cols = c('Marker Gene' , "Housekeeping Gene Lin", "Housekeeping Gene", "LPS response early", "LPS response late", "Other"),
    names_to = "category",
    values_to = "values"
  ) 

cleaned_long_df <- long_df %>%
   filter(values == TRUE) %>%
   filter(cell_type == "Macrophage")

cleaned_long_df$category <- factor(cleaned_long_df$category, 
                                   levels = c('Marker Gene', "Housekeeping Gene Lin", "Housekeeping Gene", 
                                              "LPS response early", "LPS response late", "Other"))

 boxplot_expression <- ggplot(cleaned_long_df, aes(x = category, y = gmean, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.colour = "black", outlier.alpha = 0.1, size = 0.1) +  #outlier dots
  #geom_boxplot(alpha = 0.6, outlier.shape = NA, coef = Inf, size = 0.3) + # whisters instead of dots
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(y = "Expression Level", fill = "Gene Category") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )
ggsave('plots/3_m/Test_2/boxplot4_2.png', plot = boxplot_expression, width = 8, height = 4 ) 






# TEST 3

# stacked barplot: same as verison 1 but non exclusive and more groups

cleaned_long_df <- long_df %>%
   filter(values == TRUE) 


summary_table <- cleaned_long_df %>%
  group_by(cell_type, cluster_id, gene_set, category) %>%
  summarise(count = n(), .groups = 'drop') 


#my_colors <- c("Marker Gene" = "#B8C7C7", "Other" = "#556967", "Immune Response Gene" = "#E0D5C8", "Stable Housekeeping Gene" = "#AE9982", "Housekeeping Gene" = "#5AB4AC", "Response gene not expressed" = "#A4C089")
summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


barplot_relative <- ggplot(summary_table, aes(x = gene_set, y = count, fill =  category)) +
  geom_bar(stat = "identity", position = "fill") +  # Using 'fill' for 100% stacking
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  #scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_3/stacked_barplot1_1.png', plot = barplot_relative,  width = 10, height = 14 )





# stacked barplot: same as verison 3 but with more tags

my_colors <- c("HVG" = "#A4C089", "Other"= "#556967", "LVG"= "#5AB4AC", "Response gene not expressed"= "#B8C7C7")
summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))

# Create the 100% stacked barplot
barplot_relative <- ggplot(summary_table, aes(x = category, y = count, fill =  gene_set)) +
  geom_bar(stat = "identity", position = "fill") +  # Using 'fill' for 100% stacking
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_3/stacked_barplot2.png', plot = barplot_relative,  width = 10, height = 14 )





# make stacked bar plot with longer format: expression level

long_df <- df_main_with_markers_expr_level %>%
  pivot_longer(
    cols = c('response_expressed', 'is_marker', 'housekeeping_gene_stable', 'housekeeping_gene', 'response_not_expressed', 'lps_stim', 'cytokine_gene', 'tlr_gene', 'Other'),
    values_to = 'values',
    names_to = 'category'
  )

cleaned_long_df <- long_df %>%
   filter(values == TRUE)  



summary_table <- cleaned_long_df %>%
  group_by(cell_type, cluster_id, expression_level, category) %>%
  summarise(count = n(), .groups = 'drop') 


#my_colors <- c("Marker Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
#summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))

# Create the 100% stacked barplot
barplot_relative <- ggplot(summary_table, aes(x = category, y = count, fill =  expression_level)) +
  geom_bar(stat = "identity", position = "fill") +  
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  #scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_3/stacked_barplot3.png', plot = barplot_relative,  width = 10, height = 14 )

# the response pathway groups are part of the lower quartile in expression






# TEST 4
# replaced by A8: stratification

# make stacked barplot for stratifying both: expression and gene set
stratified_df <- df_main_with_markers_expr_level %>%
  mutate(low_low = gene_set == 'LVG' & expression_level == 'Lowly Expressed') %>%
  mutate(high_high = gene_set == 'HVG' & expression_level == 'Highly Expressed') %>%
  mutate(high_low = gene_set == 'HVG' & expression_level == 'Lowly Expressed') %>%
  mutate(low_high = gene_set == 'LVG' & expression_level == 'Highly Expressed')
stratified_df

long_df1 <- stratified_df %>%
  pivot_longer(
    # no spec cytokines: cols = c('response_expressed', 'is_marker', 'housekeeping_gene_lin', 'housekeeping_gene', 'response_not_expressed', 'lps_stim', 'cytokine_gene', 'tlr_gene', 'Other'),
    c('is_marker', 'housekeeping_gene_lin', 'housekeeping_gene',   'Other')
    values_to = 'values',
    names_to = 'category')


cleaned_long_df1 <- long_df1 %>%
   filter(values == TRUE) 


long_df2 <- cleaned_long_df1 %>%
  pivot_longer(
    cols = c(
  "low_low", "high_low", "low_high", "high_high"),
   values_to = 'values2',
    names_to = 'groups')


cleaned_long_df2 <- long_df2 %>%
   filter(values2 == TRUE) 
cleaned_long_df2



summary_table <- cleaned_long_df2 %>%
  group_by(cell_type, cluster_id, category, groups) %>%
  summarise(count = n(), .groups = 'drop') 
print(summary_table, n= 100)


#my_colors <- c("Marker Gene" = "#A4C089", "Other" = "#556967", "Immune Response Gene" = "#5AB4AC")
#summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


barplot_relative <- ggplot(summary_table, aes(x = category, y = count, fill =  groups)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_4/stacked_barplot5.png', plot = barplot_relative,  width = 10, height = 14 )


# get the low_low genes
report <- cleaned_long_df2 %>%
  filter(groups == 'low_low') %>%
  filter(category %in% c('cytokine_gene', 'lps_stim', 'tlr_gene')) %>%
  #group_by(category,cell_type, cluster_id, ) %>%
  #summarise(count = n(), .groups = 'drop') 
  select(gene, cluster_id, cell_type) %>%
  unique() %>%
  group_by(gene) %>%
  summarise( count = n(), .groups = 'drop') %>%
  arrange(desc(count))

print(report, n= 100)


#write.csv(report, 'data/low_low_gene_top_100.csv')



#TEST 5

# same as test 3 but just with specific cytokines

long_df_specific <- df_main_with_markers %>%
   rename(
    `Immune Response Gene` = response_expressed,
    `Housekeeping Gene Lin` = housekeeping_gene_lin,
    `Housekeeping Gene` = housekeeping_gene, 
    `Response gene not expressed` = response_not_expressed, 
    `LPS response early` = lps_stim_early,
    `LPS response late` = lps_stim_late,
    `IFNy response` = ifng_gene,
    `Il10 response` = il10_gene,
    `Il4 response` = il4_gene,
    `TNF response` = tnf_gene,
    `Meiosis` = meiosis_gene,
    `Sperm DNA condensation` = sperm_gene, 
    `Autoimmune risk gene` = autoimmune_gene) %>%
  pivot_longer(
    cols = c('Marker Gene' , "Housekeeping Gene Lin", "Housekeeping Gene", "LPS response early", "LPS response late", "IFNy response", "Il10 response" , "Il4 response", "TNF response", "Autoimmune risk gene", "Other"),
    names_to = "category",
    values_to = "values"
  ) 
head(long_df_specific)



cleaned_long_df <- long_df_specific %>%
   filter(values == TRUE) 

summary_table <- cleaned_long_df %>%
  group_by(cell_type, cluster_id, gene_set, category) %>%
  summarise(count = n(), .groups = 'drop') 


#my_colors <- c("Marker Gene" = "#B8C7C7", "Other" = "#556967", "Immune Response Gene" = "#E0D5C8", "Stable Housekeeping Gene" = "#AE9982", "Housekeeping Gene" = "#5AB4AC", "Response gene not expressed" = "#A4C089")
summary_table$gene_set <- factor(summary_table$gene_set, levels = c("HVG", "Other", "LVG"))


barplot_relative <- ggplot(summary_table, aes(x = gene_set, y = count, fill =  category)) +
  geom_bar(stat = "identity", position = "fill") +  # Using 'fill' for 100% stacking
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  #scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_5/bar_specific_1.png', plot = barplot_relative,  width = 10, height = 14 )



barplot_relative <- ggplot(summary_table, aes(x = category, y = count, fill =  gene_set)) +
  geom_bar(stat = "identity", position = "fill") +  # Using 'fill' for 100% stacking
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  labs(x = "Gene Set",
       y = "Percentage of Genes (%)",
       fill = "Gene type") +
  #scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave('plots/3_m/Test_5/bar_specific_2.png', plot = barplot_relative,  width = 10, height = 14 )




# TEST 6
# iros code and reporting


# update 31.03.: get negative control (added meiosis and sperm signature before)




# aggregated expression

df_agg <- cleaned_long_df %>% 
  group_by(category, cluster_id, gene) 
df_agg


output_dir <- "droplet_analysis/plots/3_m/Test_6"
for (cat in unique(df_agg$category)) {
  df_subset <- df_agg %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cluster_id), y = gmean, fill = as.factor(cluster_id))) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.2,  aes(colour = as.factor(cluster_id))) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = paste("Expression in", cat),
         x = "Cluster ID", y = "Gene Expression") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  filename <- paste0(output_dir, "/violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}
colnames(cleaned_long_df)


# just macrophages
unique(cleaned_long_df$cell_type)

df_agg_macs <- cleaned_long_df %>% 
  filter(cell_type == "Macrophage") %>%
  group_by(category, cluster_id, gene) %>%
  summarise(gmean = mean(gmean, na.rm = TRUE), .groups = "drop")

output_dir <- "droplet_analysis/plots/3_m/Test_6/Macs"
for (cat in unique(df_agg_macs$category)) {
  df_subset <- df_agg_macs %>% filter(category == cat)
  p <- ggplot(df_subset, aes(x = as.factor(cluster_id), y = gmean, fill = as.factor(cluster_id))) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.2,  aes(colour = as.factor(cluster_id))) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = paste("Expression in", cat),
         x = "Cluster ID", y = "Gene Expression") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  filename <- paste0(output_dir, "/violin_", gsub(" ", "_", cat), ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}

#saveRDS(cleaned_long_df, "droplet_analysis/cleaned_long_df.rds")









# TEST 7
head(df_main_with_markers_expr_level)


# violin plots and sctrans figure reproduction

bins <- df_main_with_markers_expr_level %>%
  filter(!is.na(cluster_id)) %>%
  mutate(gmean_log = log10(gmean)) %>%
  mutate(gmean_bin = cut_number((gmean_log), n = 4))  %>%
  mutate(gmean_group = factor(gmean_bin, labels = c("Group 1", "Group 2", "Group 3", "Group 4")))
  #mutate(gmean_group = (gmean_bin))

bin_plot <- ggplot(bins, aes(x = (gmean_group), y = res_var, fill = gmean_group)) +
  geom_point(alpha = 0.2) +  
  #geom_boxplot(outlier.shape = NA,  alpha = 0.2) +  
  geom_violin(alpha = 0.6, scale = "width") + 
  #geom_jitter(width = 0.2, alpha = 0.6, color = "black") + 
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) 
  theme_classic() +
  labs(x = "Gene mean", 
       y = "Residual variance", 
       color = "Expression Bin",
       fill = "Gene expression bin") +
 # scale_fill_discrete(labels = levels(bins$gmean_bin)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('plots/3_m/Test_7/bin_plot1.png', plot = bin_plot, width = 8, height = 8, dpi = 300)




# reproduce the sctrans figure (scatterplot)
bin_plot <- ggplot(bins, aes(x = gmean, y = res_var)) +
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) 
  ylim(0, 60) +
  theme_classic() +
  labs(x = "Gene mean", 
       y = "Residual variance", 
       color = "Expression Bin",
       fill = "Gene expression bin") +

  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('plots/3_m/Test_7/bin_plot2.png', plot = bin_plot, width = 8, height = 8, dpi = 300)



# add the marker genes 
sum(bins$ifng_gene)


bin_plot <- ggplot(bins, aes(x = gmean, y = res_var)) +
  geom_point(aes(color = housekeeping_gene), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~paste(cell_type, cluster_id, sep = " - "), scales = "free", labeller = label_parsed) +
  ylim(0, 60) +
  theme_classic() +
  scale_color_manual(values = c("TRUE" = "#9A1463", "FALSE" = "black"), labels = c("TRUE" = "Marker Genes", "FALSE" = "Other Genes")) +  
  labs(x = "Gene mean", 
       y = "Residual variance", 
       color = "Expression Bin",
       fill = "Gene expression bin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('plots/3_m/Test_7/bin_plot4.png', plot = bin_plot, width = 8, height = 8, dpi = 300)

