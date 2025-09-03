
# complex heatmap

heatmap_data <- df_main_with_markers %>%
  select(gene, gmean, cluster_id) %>%
  pivot_wider(names_from = cluster_id, values_from = gmean) 
print(heatmap_data, n = 30)


heatmap_matrix <- heatmap_data %>% 
  select(c(-gene) %>%
  as.matrix()
rownames(heatmap_matrix) <- heatmap_data$gene
heatmap_matrix[is.na(heatmap_matrix)] <- 0
head(heatmap_matrix)

df_annotation <- df_main_with_markers %>%
  select(gene_set, gene)
df_annotation

colors <- colorRampPalette(c("navy", "white", "firebrick"))(100)

heatmap <- Heatmap(
  heatmap_matrix,
  col = colors,
  show_row_names = FALSE,
  show_column_names = TRUE
  bottom_annotation = bottom_annotation)
  #row_annotation = row_annotation,

png('plots/3_m/heatmap_lineage_expression_complex.png', width = 1000, height = 800, res = 300)  
draw(heatmap)
dev.off()



# reporting ggplot
heatmap_plot <- ggplot(df_reporting, aes(x = gene, y = mouse_id, fill = mean_expression)) +
  geom_tile() +
  geom_text(aes(label = gene_set), color = "black", size = 3, vjust = 1) + 
  scale_fill_gradientn(colors = colors, values = scales::rescale(breakpoints), name = "Mean Gene Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank()
  ) 

ggsave('plots/3_m/heatmap_lineage_expression_3m.png', plot = heatmap_plot, width = 40, height = 10, dpi = 300)

# to make it wider, add to ggsave
#width = 8, height = 4, dpi = 150      


