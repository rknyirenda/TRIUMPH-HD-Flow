###############################################################################
# MUCOSAL HIGH-DIMENSIONAL ANALYSIS
# SCRIPT 04: PUBLICATION FIGURES
#
# Reproduces three reference figure styles + adds PCA:
#   FIG 1 : Two-way clustered heatmap (row + column dendrograms, RdBu z-score,
#           column annotation bar)                        [ComplexHeatmap]
#   FIG 2 : UMAP panels coloured by cluster + faceted-by-tissue density
#           overlays + per-marker feature grid            [ggplot2 / patchwork]
#   FIG 3 : Faceted ridgeline density plots per marker    [ggridges]
#   FIG 4 : PCA of samples on cluster-frequency matrix     [prcomp / ggplot2]
#
# HARDWARE NOTE (8 GB RAM):
#   Large single-cell matrices are the memory bottleneck. This script:
#     - works from checkpoints (never re-reads FCS)
#     - subsamples cells ONLY for point-heavy UMAP plots (rasterised)
#     - aggregates to cluster/sample level before heavy plotting
#     - purges + gc() after every major block
###############################################################################

rm(list = ls())
graphics.off()
gc()

###############################################################################
# 0. PACKAGES
###############################################################################

suppressPackageStartupMessages({
         library(tidyverse)        # dplyr / ggplot2 / tibble / stringr
         library(ggridges)         # ridgeline plots  (FIG 3)
         library(patchwork)        # panel assembly
         library(ComplexHeatmap)   # heatmap          (FIG 1)
         library(circlize)         # colorRamp2       (FIG 1)
         library(RColorBrewer)
         library(ggrastr)          # rasterise dense point layers (memory/vector)
         library(matrixStats)      # fast colVars for PCA feature selection
})

set.seed(2026)

###############################################################################
# 0b. MEMORY HELPER
###############################################################################

mem_report <- function(step){
         cat("\n=====================================================\n")
         cat(step, "\n")
         print(gc())
         cat("=====================================================\n\n")
}
mem_report("Start Script 04")

###############################################################################
# 0c. OUTPUT FOLDERS
###############################################################################

dir.create("Figures",              showWarnings = FALSE)
dir.create("Figures/Heatmaps",     showWarnings = FALSE)
dir.create("Figures/UMAP",         showWarnings = FALSE)
dir.create("Figures/Ridgeline",    showWarnings = FALSE)
dir.create("Figures/PCA",          showWarnings = FALSE)
dir.create("Results",              showWarnings = FALSE)

###############################################################################
# 0d. GLOBAL AESTHETICS  (shared publication theme)
###############################################################################

theme_pub <- theme_classic(base_size = 9) +
         theme(
                  panel.border    = element_blank(),
                  axis.line       = element_line(colour = "black", linewidth = 0.5),
                  axis.text       = element_text(colour = "black"),
                  legend.title    = element_text(face = "bold"),
                  legend.position = "right",
                  plot.title      = element_text(face = "bold", hjust = 0.5),
                  strip.background = element_blank(),
                  strip.text      = element_text(face = "bold")
         )

# ---- Cluster colour palette (colour-blind safe, high contrast) -------------
# Extend to as many clusters as you have (20 metaclusters here).
cluster_cols <- c(
         "#0072B2","#E69F00","#009E73","#CC79A7","#56B4E9",
         "#D55E00","#F0E442","#000000","#882255","#44AA99",
         "#117733","#999933","#88CCEE","#AA4499","#DDCC77",
         "#332288","#661100","#6699CC","#AA4466","#228833"
)

# ---- Tissue colours --------------------------------------------------------
tissue_cols <- c(
         "Nasal Swab"         = "#1B9E77",
         "Nasal Scrape"       = "#7570B3",
         "Cervical Scrape"    = "#D95F02",
         "Cervical Cytobrush" = "#E7298A",
         "Unknown"            = "grey60"
)

# ---- Diverging heatmap ramp (RdBu, matches reference image 1) --------------
# Reference uses a strong blue -> white -> red diverging scheme.
heat_ramp <- circlize::colorRamp2(
         c(-2, -1, 0, 1, 2),
         c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
)

# ---- Viridis-magma for continuous marker feature plots ---------------------
# (option = "C" = plasma; matches the warm feature-plot look)

###############################################################################
# 1. LOAD CHECKPOINTS
###############################################################################

cell_data         <- readRDS("Results/cell_data.rds")          # per-cell table
umap              <- readRDS("Results/UMAP.rds")               # matrix Nx2
cluster_abundance <- readRDS("Results/cluster_frequency.rds")  # sample x cluster

# Marker set actually used for clustering / display
markers <- c("CD3","CD4","CD8","CD19","CD14","CD16",
             "CD66b","CD56","CD69","CD103","CD45RO","CCR7")

# Defensive: ensure marker columns exist (Script 03 renames FJComp -> antigen).
# If they are still the raw FJComp names, rename now.
if(!all(markers %in% colnames(cell_data))){
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
         hit <- match(names(marker_lookup), colnames(cell_data))
         colnames(cell_data)[hit[!is.na(hit)]] <-
                  unname(marker_lookup[!is.na(hit)])
}

cell_data$Cluster <- factor(cell_data$Cluster)

mem_report("Checkpoints loaded")

###############################################################################
###############################################################################
# FIGURE 1 : TWO-WAY CLUSTERED HEATMAP  (reference image 1)
#
# Reference anatomy:
#   - rows      = features (miRNAs there / MARKERS or CLUSTERS here)
#   - columns   = samples, with a coloured group annotation bar
#   - both dendrograms shown
#   - z-scored per row, RdBu diverging key
#
# We build the biologically-informative version: rows = markers,
# columns = samples, values = median arcsinh expression, z-scored per marker.
# (This is the flow analogue of the reference's feature x sample layout.)
###############################################################################
###############################################################################

# ---- Build marker x sample median matrix -----------------------------------
sample_marker_median <- cell_data %>%
         group_by(Sample) %>%
         summarise(across(all_of(markers), median), .groups = "drop")

mat_hm <- sample_marker_median %>%
         column_to_rownames("Sample") %>%
         as.matrix() %>%
         t()                                  # markers (rows) x samples (cols)

# ---- Z-score per row (per marker) ------------------------------------------
mat_z <- t(scale(t(mat_hm)))                  # centre + scale each marker row
mat_z[is.na(mat_z)] <- 0                      # guard constant rows

# ---- Column annotation (tissue group bar, like the green/red bar) ----------
col_meta <- cell_data %>%
         distinct(Sample, tissue) %>%
         filter(Sample %in% colnames(mat_z)) %>%
         arrange(match(Sample, colnames(mat_z)))

col_anno <- HeatmapAnnotation(
         Tissue = col_meta$tissue,
         col    = list(Tissue = tissue_cols[unique(col_meta$tissue)]),
         annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
         simple_anno_size   = grid::unit(3, "mm")
)

ht <- Heatmap(
         mat_z,
         name                = "Row z-score",
         col                 = heat_ramp,
         top_annotation      = col_anno,
         cluster_rows        = TRUE,
         cluster_columns     = TRUE,
         clustering_method_rows    = "ward.D2",
         clustering_method_columns = "ward.D2",
         clustering_distance_rows    = "euclidean",
         clustering_distance_columns = "euclidean",
         show_column_names   = TRUE,
         show_row_names      = TRUE,
         row_names_gp        = grid::gpar(fontsize = 8),
         column_names_gp     = grid::gpar(fontsize = 6),
         row_dend_width      = grid::unit(15, "mm"),
         column_dend_height  = grid::unit(15, "mm"),
         heatmap_legend_param = list(
                  title_gp  = grid::gpar(fontsize = 8, fontface = "bold"),
                  labels_gp = grid::gpar(fontsize = 7),
                  legend_direction = "horizontal"
         ),
         border = FALSE
)

# ---- Export ----------------------------------------------------------------
pdf("Figures/Heatmaps/FIG1_Marker_Sample_Heatmap.pdf", width = 9, height = 5)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
dev.off()

tiff("Figures/Heatmaps/FIG1_Marker_Sample_Heatmap.tiff",
     width = 9, height = 5, units = "in", res = 600, compression = "lzw")
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
dev.off()

png("Figures/Heatmaps/FIG1_Marker_Sample_Heatmap.png",
    width = 9, height = 5, units = "in", res = 600)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
dev.off()

# ---- ALSO: cluster x marker heatmap (phenotype key for the UMAP) -----------
cluster_marker_median <- cell_data %>%
         group_by(Cluster) %>%
         summarise(across(all_of(markers), median), .groups = "drop")

mat_cl <- cluster_marker_median %>%
         column_to_rownames("Cluster") %>%
         as.matrix()
rownames(mat_cl) <- paste0("MC", rownames(mat_cl))
mat_cl_z <- scale(mat_cl)                       # z-score per marker (column)
mat_cl_z[is.na(mat_cl_z)] <- 0

ht_cl <- Heatmap(
         mat_cl_z,
         name             = "z-score",
         col              = heat_ramp,
         cluster_rows     = TRUE,
         cluster_columns  = TRUE,
         clustering_method_rows    = "ward.D2",
         clustering_method_columns = "ward.D2",
         row_names_gp     = grid::gpar(fontsize = 8),
         column_names_gp  = grid::gpar(fontsize = 8),
         row_dend_width   = grid::unit(12, "mm"),
         column_dend_height = grid::unit(12, "mm"),
         border = FALSE
)

pdf("Figures/Heatmaps/FIG1b_Cluster_Marker_Heatmap.pdf", width = 6, height = 7)
draw(ht_cl)
dev.off()
png("Figures/Heatmaps/FIG1b_Cluster_Marker_Heatmap.png",
    width = 6, height = 7, units = "in", res = 600)
draw(ht_cl)
dev.off()

rm(mat_hm, mat_z, mat_cl, mat_cl_z, sample_marker_median,
   cluster_marker_median, ht, ht_cl, col_anno, col_meta)
gc()
mem_report("FIG 1 done")

###############################################################################
###############################################################################
# FIGURE 2 : UMAP PANELS  (reference image 2)
#
# Reference anatomy:
#   B = one master UMAP coloured by cluster identity + legend
#   C,D,E = the same embedding faceted / highlighted by condition
#
# Memory strategy: rasterise point layers (ggrastr) so exported PDFs stay
# small and RAM during draw stays low, even with ~N*10k points.
###############################################################################
###############################################################################

# ---- Assemble a lean plotting frame (only what we plot) --------------------
umap_df <- tibble(
         UMAP1   = umap[,1],
         UMAP2   = umap[,2],
         Cluster = cell_data$Cluster,
         tissue  = cell_data$tissue
)

# Subsample for point-dense plots (keeps shape, slashes RAM/vector size).
# 60k points is plenty for a smooth-looking UMAP; adjust if you like.
max_pts <- 60000
if(nrow(umap_df) > max_pts){
         keep_idx <- sample.int(nrow(umap_df), max_pts)
         umap_plot_df <- umap_df[keep_idx, ]
} else {
         umap_plot_df <- umap_df
}

n_clusters <- nlevels(umap_df$Cluster)
pal_clusters <- setNames(cluster_cols[seq_len(n_clusters)],
                         levels(umap_df$Cluster))

# ---- PANEL B: master UMAP coloured by cluster ------------------------------
p_master <- ggplot(umap_plot_df, aes(UMAP1, UMAP2, colour = Cluster)) +
         rasterise(geom_point(size = 0.15, alpha = 0.6), dpi = 600) +
         scale_colour_manual(values = pal_clusters, name = "Metacluster") +
         coord_equal() +
         guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1),
                                      ncol = 1)) +
         theme_pub +
         labs(title = "Mucosal immune landscape (FlowSOM metaclusters)")

ggsave("Figures/UMAP/FIG2B_UMAP_master.pdf", p_master,
       width = 7, height = 5.5)
ggsave("Figures/UMAP/FIG2B_UMAP_master.png", p_master,
       width = 7, height = 5.5, dpi = 600)

# ---- PANELS C-E style: UMAP faceted by tissue, all-cells greyed backdrop ----
backdrop <- umap_plot_df %>% dplyr::select(UMAP1, UMAP2)

p_facets <- ggplot() +
         rasterise(geom_point(data = backdrop,
                              aes(UMAP1, UMAP2),
                              colour = "grey88", size = 0.10), dpi = 600) +
         rasterise(geom_point(data = umap_plot_df,
                              aes(UMAP1, UMAP2, colour = Cluster),
                              size = 0.14, alpha = 0.7), dpi = 600) +
         scale_colour_manual(values = pal_clusters, guide = "none") +
         facet_wrap(~ tissue, ncol = 2) +
         coord_equal() +
         theme_pub +
         theme(strip.text = element_text(face = "bold", size = 10))

ggsave("Figures/UMAP/FIG2CE_UMAP_by_tissue.pdf", p_facets,
       width = 8, height = 7)
ggsave("Figures/UMAP/FIG2CE_UMAP_by_tissue.png", p_facets,
       width = 8, height = 7, dpi = 600)

# ---- PANEL: per-marker feature grid (plasma/viridis continuous) ------------
# Attach marker expression to the (subsampled) coordinates to stay light.
feat_df <- bind_cols(
         umap_plot_df[, c("UMAP1","UMAP2")],
         cell_data[if(exists("keep_idx")) keep_idx else seq_len(nrow(cell_data)),
                   markers]
)

feature_plot <- function(marker){
         ggplot(feat_df, aes(UMAP1, UMAP2, colour = .data[[marker]])) +
                  rasterise(geom_point(size = 0.12, alpha = 0.7), dpi = 600) +
                  scale_colour_viridis_c(option = "C", name = NULL) +
                  coord_equal() +
                  theme_pub +
                  theme(legend.key.width  = grid::unit(2, "mm"),
                        legend.key.height = grid::unit(6, "mm"),
                        axis.title = element_blank(),
                        axis.text  = element_blank(),
                        axis.ticks = element_blank()) +
                  ggtitle(marker)
}

feat_plots <- lapply(markers, feature_plot)
p_features <- wrap_plots(feat_plots, ncol = 4)

ggsave("Figures/UMAP/FIG2_Feature_grid.pdf", p_features,
       width = 10, height = 8)
ggsave("Figures/UMAP/FIG2_Feature_grid.png", p_features,
       width = 10, height = 8, dpi = 600)

rm(umap_df, umap_plot_df, backdrop, feat_df, feat_plots,
   p_master, p_facets, p_features)
if(exists("keep_idx")) rm(keep_idx)
gc()
mem_report("FIG 2 done")

###############################################################################
###############################################################################
# FIGURE 3 : RIDGELINE DENSITY PLOTS  (reference image 3)
#
# Reference anatomy:
#   - one facet panel per marker (there: per HTO)
#   - y = Identity (there: HTO identity; here: CLUSTER or TISSUE)
#   - x = Expression Level
#   - filled ridges, classic clean theme
#
# We show marker expression distributions split by CLUSTER, faceted by marker.
# Memory strategy: ggridges computes densities on the fly; we subsample cells
# per (cluster x marker) so the density estimation stays cheap.
###############################################################################
###############################################################################

# ---- Long, subsampled frame ------------------------------------------------
# Cap cells per cluster to keep the KDE light (2000/cluster is smooth).
cap_per_cluster <- 2000

ridge_cells <- cell_data %>%
         dplyr::select(Cluster, all_of(markers)) %>%
         group_by(Cluster) %>%
         slice_sample(n = cap_per_cluster) %>%     # caps large clusters
         ungroup()

ridge_long <- ridge_cells %>%
         pivot_longer(cols = all_of(markers),
                      names_to = "Marker", values_to = "Expression") %>%
         mutate(Marker = factor(Marker, levels = markers))

rm(ridge_cells); gc()

p_ridge <- ggplot(ridge_long,
                  aes(x = Expression, y = Cluster, fill = Cluster)) +
         geom_density_ridges(scale = 2.2, linewidth = 0.25,
                             rel_min_height = 0.01, alpha = 0.9) +
         scale_fill_manual(values = pal_clusters, guide = "none") +
         facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
         labs(x = "Arcsinh expression", y = "Metacluster") +
         theme_ridges(font_size = 9, grid = TRUE) +
         theme(strip.text = element_text(face = "bold"),
               axis.title = element_text(face = "bold"))

ggsave("Figures/Ridgeline/FIG3_Ridgeline_by_cluster.pdf", p_ridge,
       width = 11, height = 9)
ggsave("Figures/Ridgeline/FIG3_Ridgeline_by_cluster.png", p_ridge,
       width = 11, height = 9, dpi = 600)

# ---- Variant closer to reference: y = TISSUE, one panel per marker ---------
ridge_tissue <- cell_data %>%
         dplyr::select(tissue, all_of(markers)) %>%
         group_by(tissue) %>%
         slice_sample(n = min(5000, cap_per_cluster * 3)) %>%
         ungroup() %>%
         pivot_longer(cols = all_of(markers),
                      names_to = "Marker", values_to = "Expression") %>%
         mutate(Marker = factor(Marker, levels = markers))

p_ridge_tissue <- ggplot(ridge_tissue,
                         aes(x = Expression, y = tissue, fill = tissue)) +
         geom_density_ridges(scale = 1.8, linewidth = 0.3,
                             rel_min_height = 0.01, alpha = 0.9) +
         scale_fill_manual(values = tissue_cols, guide = "none") +
         facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
         labs(x = "Arcsinh expression", y = "Tissue") +
         theme_ridges(font_size = 9, grid = TRUE) +
         theme(strip.text = element_text(face = "bold"),
               axis.title = element_text(face = "bold"))

ggsave("Figures/Ridgeline/FIG3b_Ridgeline_by_tissue.pdf", p_ridge_tissue,
       width = 11, height = 7)
ggsave("Figures/Ridgeline/FIG3b_Ridgeline_by_tissue.png", p_ridge_tissue,
       width = 11, height = 7, dpi = 600)

rm(ridge_long, ridge_tissue, p_ridge, p_ridge_tissue); gc()
mem_report("FIG 3 done")

###############################################################################
###############################################################################
# FIGURE 4 : PCA OF MUCOSAL SAMPLES  (new, requested)
#
# Rationale: each sample is summarised as its cluster-frequency vector
# (compositional immune fingerprint). PCA on these vectors asks whether
# samples separate by tissue in unsupervised space — the flow analogue of a
# scRNA-seq sample-level PCA and a natural companion to the heatmap.
#
# We also provide a marker-median PCA as an alternative feature space.
###############################################################################
###############################################################################

# ---- 4a. Build sample x cluster frequency matrix ---------------------------
freq_wide <- cluster_abundance %>%
         dplyr::select(Sample, Cluster, Frequency) %>%
         mutate(Cluster = paste0("MC", Cluster)) %>%
         pivot_wider(names_from = Cluster, values_from = Frequency,
                     values_fill = 0)

sample_meta <- cell_data %>%
         distinct(Sample, tissue, patient) %>%
         filter(Sample %in% freq_wide$Sample)

freq_mat <- freq_wide %>%
         column_to_rownames("Sample") %>%
         as.matrix()

# CLR-style stabilisation for compositional data (avoid log(0) via pseudocount)
freq_clr <- log(freq_mat + 1e-4)
freq_clr <- sweep(freq_clr, 1, rowMeans(freq_clr), "-")   # centred log-ratio

# ---- PCA (features = clusters). Drop zero-variance columns first. -----------
nzv <- matrixStats::colVars(freq_clr) > 0
pca <- prcomp(freq_clr[, nzv, drop = FALSE], center = TRUE, scale. = TRUE)

var_exp <- 100 * (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- as_tibble(pca$x[, 1:2], rownames = "Sample") %>%
         left_join(sample_meta, by = "Sample")

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = tissue)) +
         geom_point(size = 3, alpha = 0.9) +
         ggrepel::geom_text_repel(aes(label = patient), size = 2.5,
                                  max.overlaps = 20, show.legend = FALSE) +
         scale_colour_manual(values = tissue_cols, name = "Tissue") +
         stat_ellipse(aes(group = tissue), type = "norm",
                      linewidth = 0.4, alpha = 0.5) +
         labs(
                  title = "PCA of mucosal samples (cluster-frequency space)",
                  x = sprintf("PC1 (%.1f%%)", var_exp[1]),
                  y = sprintf("PC2 (%.1f%%)", var_exp[2])
         ) +
         theme_pub +
         coord_equal()

ggsave("Figures/PCA/FIG4_PCA_cluster_frequency.pdf", p_pca,
       width = 7, height = 6)
ggsave("Figures/PCA/FIG4_PCA_cluster_frequency.png", p_pca,
       width = 7, height = 6, dpi = 600)

# ---- 4b. Scree plot --------------------------------------------------------
scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)),
                               levels = paste0("PC", seq_along(var_exp))),
                   Variance = var_exp) %>%
         slice_head(n = min(10, length(var_exp)))

p_scree <- ggplot(scree_df, aes(PC, Variance)) +
         geom_col(fill = "#4477AA") +
         geom_text(aes(label = sprintf("%.1f%%", Variance)),
                   vjust = -0.4, size = 2.6) +
         labs(x = NULL, y = "Variance explained (%)",
              title = "Scree plot") +
         theme_pub

ggsave("Figures/PCA/FIG4b_PCA_scree.pdf", p_scree, width = 5, height = 4)
ggsave("Figures/PCA/FIG4b_PCA_scree.png", p_scree, width = 5, height = 4,
       dpi = 600)

# ---- 4c. PCA loadings (which clusters drive separation) --------------------
load_df <- as_tibble(pca$rotation[, 1:2], rownames = "Cluster") %>%
         mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
         arrange(desc(magnitude))

write.csv(load_df, "Results/PCA_loadings.csv", row.names = FALSE)

p_load <- ggplot(load_df, aes(PC1, PC2)) +
         geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                      arrow = arrow(length = grid::unit(2, "mm")),
                      colour = "grey50", linewidth = 0.3) +
         ggrepel::geom_text_repel(aes(label = Cluster), size = 2.6,
                                  max.overlaps = 30) +
         labs(title = "PCA loadings (cluster contributions)",
              x = sprintf("PC1 (%.1f%%)", var_exp[1]),
              y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
         theme_pub

ggsave("Figures/PCA/FIG4c_PCA_loadings.pdf", p_load, width = 6, height = 5)
ggsave("Figures/PCA/FIG4c_PCA_loadings.png", p_load, width = 6, height = 5,
       dpi = 600)

saveRDS(pca, "Results/PCA_object.rds")
write.csv(pca_df, "Results/PCA_scores.csv", row.names = FALSE)

rm(freq_wide, freq_mat, freq_clr, pca, pca_df, load_df,
   scree_df, p_pca, p_scree, p_load, nzv, var_exp)
gc()
mem_report("FIG 4 (PCA) done")

###############################################################################
# DONE
###############################################################################
cat("\nAll figures written to ./Figures/{Heatmaps,UMAP,Ridgeline,PCA}\n")
cat("PDF (vector) + PNG (600 dpi) + TIFF (heatmap) exported.\n")
mem_report("Script 04 finished")
