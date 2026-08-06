###############################################################################
# BUILD MASTER CELL TABLE
###############################################################################

cell_data <- as.data.frame(expr_all)

cell_data$Cluster <- factor(cell_meta20)

cell_data$Sample <- sample_id

###############################################################################
# Add sample metadata
###############################################################################

cell_data <- left_join(
         
         cell_data,
         
         metadata,
         
         by = c("Sample" = "sample_id")
         
)

str(cell_data)

dim(cell_data)


library(dplyr)

cluster_medians <-
         
         cell_data %>%
         
         group_by(Cluster) %>%
         
         summarise(
                  
                  across(
                           
                           all_of(cluster_channels),
                           
                           median
                           
                  )
                  
         )

cluster_medians


library(pheatmap)

heat <- cluster_medians %>%
         
         column_to_rownames("Cluster")

heat <- as.matrix(heat)

pheatmap(
         
         heat,
         
         scale = "row",
         
         clustering_method = "ward.D2",
         
         border_color = NA,
         
         fontsize = 10,
         
         angle_col = 45
         
)

cluster_freq <-
         
         cell_data %>%
         
         count(
                  
                  Sample,
                  
                  Cluster
                  
         ) %>%
         
         group_by(Sample) %>%
         
         mutate(
                  
                  Frequency = n / sum(n)
                  
         )

cluster_freq
cluster_freq <- left_join(
         
         cluster_freq,
         
         metadata,
         
         by = c("Sample" = "sample_id")
         
)


library(ggplot2)

ggplot(
         
         cluster_freq,
         
         aes(
                  
                  Cluster,
                  
                  Frequency,
                  
                  fill = tissue
                  
         )
         
) +
         
         geom_boxplot(
                  
                  outlier.shape = NA
                  
         ) +
         
         geom_jitter(
                  
                  width = 0.2,
                  
                  alpha = 0.6,
                  
                  size = 2
                  
         ) +
         
         theme_classic(base_size = 14)


saveRDS(
         
         cell_data,
         
         "Results/cell_data.rds"
         
)

saveRDS(
         
         cluster_freq,
         
         "Results/cluster_frequency.rds"
         
)


###############################################################################
# Rename fluorescence channels
###############################################################################

marker_lookup <- c(
         
         "FJComp-BV785-A-1"        = "CD3",
         "FJComp-APC-Cy7-A-1"      = "CD4",
         "FJComp-AF700-A-1"        = "CD8",
         "FJComp-BV750-A-1"        = "CD19",
         "FJComp-BV421-A-1"        = "CD14",
         "FJComp-BV650-A-1"        = "CD16",
         "FJComp-PE-A-1"           = "CD66b",
         "FJComp-APC-A-1"          = "CD56",
         "FJComp-BV605-A-1"        = "CD69",
         "FJComp-PE-Dazzle594-A-1" = "CD103",
         "FJComp-FITC-A-1"         = "CD45RO",
         "FJComp-PE-Cy7-A-1"       = "CCR7"
         
)

names(cell_data)[match(names(marker_lookup), names(cell_data))] <-
         unname(marker_lookup)

gc()
colnames(cell_data)


library(dplyr)

cluster_summary <-
         
         cell_data %>%
         
         group_by(Cluster) %>%
         
         summarise(
                  
                  Cells = n(),
                  
                  across(
                           
                           CD3:CCR7,
                           
                           list(
                                    
                                    Median = median,
                                    
                                    Mean = mean,
                                    
                                    Positive = ~100*mean(. > 1)
                                    
                           ),
                           
                           .names = "{.col}_{.fn}"
                           
                  )
                  
         )

cluster_summary
write.csv(
         cluster_summary,
         "Results/Cluster_Summary.csv",
         row.names = FALSE
)
library(ComplexHeatmap)
library(circlize)

heat_matrix <- cluster_summary |>
         dplyr::select(Cluster, ends_with("_Median"))

rownames(heat_matrix) <- paste0("MC", heat_matrix$Cluster)

heat_matrix <- as.matrix(heat_matrix[,-1])

Heatmap(
         
         heat_matrix,
         
         name = "Median\nExpression",
         
         cluster_rows = TRUE,
         
         cluster_columns = FALSE,
         
         row_names_gp = grid::gpar(fontsize = 9),
         
         column_names_gp = grid::gpar(fontsize = 10),
         
         border = FALSE
         
)

plot_df <- data.frame(
         
         UMAP1 = umap[,1],
         
         UMAP2 = umap[,2],
         
         cell_data[, c(
                  "CD3",
                  "CD4",
                  "CD8",
                  "CD19",
                  "CD14",
                  "CD16",
                  "CD66b",
                  "CD56",
                  "CD69",
                  "CD103",
                  "CD45RO",
                  "CCR7"
         )]
         
)

library(ggplot2)
library(patchwork)

feature_plot <- function(marker){
         
         ggplot(plot_df,
                aes(
                         UMAP1,
                         UMAP2,
                         colour = .data[[marker]]
                )) +
                  
                  geom_point(
                           size = 0.12,
                           alpha = 0.6
                  ) +
                  
                  scale_colour_viridis_c() +
                  
                  coord_equal() +
                  
                  theme_classic(base_size = 12) +
                  
                  ggtitle(marker)
         
}


plots <- lapply(
         
         c(
                  "CD3","CD4","CD8",
                  "CD19","CD14","CD16",
                  "CD66b","CD56",
                  "CD69","CD103",
                  "CD45RO","CCR7"
         ),
         
         feature_plot
         
)

wrap_plots(plots, ncol = 4)


###############################################################################
# Cluster sizes
###############################################################################

cluster_size <- cell_data %>%
         count(Cluster, name = "Cells") %>%
         mutate(
                  Percent = round(100 * Cells / sum(Cells), 2)
         ) %>%
         arrange(desc(Cells))

cluster_size

write.csv(
         cluster_size,
         "Results/Cluster_Size.csv",
         row.names = FALSE
)

cluster_annotation <- cluster_summary %>%
         left_join(cluster_size, by = "Cluster") %>%
         arrange(desc(Cells))

cluster_annotation

write.csv(
         cluster_annotation,
         "Results/Cluster_Annotation_Table.csv",
         row.names = FALSE
)


library(ComplexHeatmap)
library(circlize)

heat_matrix <- cluster_annotation %>%
         select(
                  Cluster,
                  CD3_Median,
                  CD4_Median,
                  CD8_Median,
                  CD19_Median,
                  CD14_Median,
                  CD16_Median,
                  CD66b_Median,
                  CD56_Median,
                  CD69_Median,
                  CD103_Median,
                  CD45RO_Median,
                  CCR7_Median
         )

rownames(heat_matrix) <-
         paste0("MC", heat_matrix$Cluster,
                " (", cluster_annotation$Percent, "%)")

heat_matrix <- as.matrix(heat_matrix[,-1])

colnames(heat_matrix) <-
         c("CD3","CD4","CD8","CD19",
           "CD14","CD16","CD66b","CD56",
           "CD69","CD103","CD45RO","CCR7")

Heatmap(
         
         heat_matrix,
         
         name = "Median",
         
         cluster_rows = TRUE,
         
         cluster_columns = FALSE,
         
         row_names_gp = grid::gpar(fontsize = 9),
         
         column_names_gp = grid::gpar(fontsize = 10)
         
)

library(patchwork)
library(ggplot2)

cluster_plots <- lapply(levels(cell_data$Cluster), function(cl){
         
         ggplot(plot_df,
                aes(UMAP1, UMAP2)) +
                  
                  geom_point(
                           colour = "grey90",
                           size = 0.08
                  ) +
                  
                  geom_point(
                           data = subset(plot_df,
                                         Cluster == cl),
                           aes(UMAP1, UMAP2),
                           colour = "#D81B60",
                           size = 0.15
                  ) +
                  
                  coord_equal() +
                  
                  theme_void() +
                  
                  ggtitle(paste("MC", cl))
         
})

wrap_plots(cluster_plots, ncol = 4)


markers <- c(
         "CD3","CD4","CD8",
         "CD19","CD14","CD16",
         "CD66b","CD56",
         "CD69","CD103",
         "CD45RO","CCR7"
)

feature_plot <- function(marker){
         
         ggplot(plot_df,
                aes(
                         UMAP1,
                         UMAP2,
                         colour = .data[[marker]]
                )) +
                  
                  geom_point(
                           size = 0.12,
                           alpha = 0.6
                  ) +
                  
                  scale_colour_viridis_c(
                           option = "C"
                  ) +
                  
                  coord_equal() +
                  
                  theme_classic(base_size = 12) +
                  
                  ggtitle(marker)
         
}

plots <- lapply(markers, feature_plot)

wrap_plots(plots, ncol = 4)


cluster_abundance <- cell_data %>%
         count(
                  Sample,
                  Cluster,
                  tissue
         ) %>%
         group_by(Sample) %>%
         mutate(
                  Frequency = n / sum(n)
         ) %>%
         ungroup()

write.csv(
         cluster_abundance,
         "Results/Cluster_Abundance.csv",
         row.names = FALSE
)


cluster_abundance <- cell_data %>%
  count(
    Sample,
    Cluster,
    tissue
  ) %>%
  group_by(Sample) %>%
  mutate(
    Frequency = n / sum(n)
  ) %>%
  ungroup()

write.csv(
  cluster_abundance,
  "Results/Cluster_Abundance.csv",
  row.names = FALSE
)
###############################################################################
# Create folders
###############################################################################

dir.create("Figures", showWarnings = FALSE)
dir.create("Figures/UMAP", showWarnings = FALSE)
dir.create("Figures/Heatmaps", showWarnings = FALSE)
dir.create("Figures/FeaturePlots", showWarnings = FALSE)
dir.create("Figures/Clusters", showWarnings = FALSE)
dir.create("Results", showWarnings = FALSE)
library(ComplexHeatmap)

ht <- Heatmap(
         
         heat_matrix,
         
         name = "Median",
         
         cluster_rows = TRUE,
         
         cluster_columns = FALSE,
         
         row_names_gp = grid::gpar(fontsize = 9),
         
         column_names_gp = grid::gpar(fontsize = 10)
         
)

pdf(
         "Figures/Heatmaps/Cluster_Median_Heatmap.pdf",
         width = 8,
         height = 10
)

draw(ht)

dev.off()

png(
         "Figures/Heatmaps/Cluster_Median_Heatmap.png",
         width = 2400,
         height = 3000,
         res = 300
)

draw(ht)

dev.off()


###############################################################################
# Save UMAP
###############################################################################

umap_plot <-
         
         ggplot(plot_df,
                aes(UMAP1,
                    UMAP2,
                    colour = Cluster)) +
         
         geom_point(
                  size = 0.12,
                  alpha = 0.6
         ) +
         
         coord_equal() +
         
         theme_classic(base_size = 14)

ggsave(
         
         filename = "Figures/UMAP/UMAP_Clusters.pdf",
         
         plot = umap_plot,
         
         width = 8,
         
         height = 6
         
)

ggsave(
         
         filename = "Figures/UMAP/UMAP_Clusters.png",
         
         plot = umap_plot,
         
         width = 8,
         
         height = 6,
         
         dpi = 600
         
)
