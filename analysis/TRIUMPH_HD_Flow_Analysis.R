###############################################################################
#                                                                             #
#   MUCOSAL HIGH-DIMENSIONAL IMMUNE PROFILING — MASTER PIPELINE                #
#   Spectral flow cytometry (Sony ID7000) | Mucosal compartments              #
#                                                                             #
#   Author : Robert Nyirenda (MLW, Infection & Immunity Group)                #
#                                                                             #
#   STAGES                                                                     #
#     01  Import, QC, channel curation, panel annotation                      #
#     02  Filter (>=10k cells) & downsample to 10k / sample                    #
#     03  Arcsinh transform, FlowSOM + ConsensusClusterPlus, UMAP             #
#     04  Master cell table, cluster summaries, abundances                     #
#     05  Publication figures (heatmaps, UMAP, ridgelines, PCA)               #
#                                                                             #
#   DESIGN NOTES                                                               #
#     - Checkpoint-driven: each stage saves .rds and can resume alone.         #
#       To re-run only figures, set RUN_STAGES <- 5 (loads stage-04 output).  #
#     - Tuned for 8 GB RAM: aggressive rm() + gc(), downsampling, and         #
#       rasterised point layers for figures.                                   #
#     - One source of truth for panel / markers / palettes / theme (Sec 1).   #
#                                                                             #
###############################################################################


###############################################################################
# SECTION 0 : ENVIRONMENT
###############################################################################

setwd("C:/Users/LENOVO/OneDrive - Malawi-Liverpool Wellcome Research Programme/Desktop/LAB WORK/Exported FCS Files for Mucosal Analysis/High Dimensionality reductionl-2026/DownSample")

rm(list = ls())
graphics.off()
gc()

suppressPackageStartupMessages({
         ## Stage 01-04 (data)
         library(flowCore)
         library(FlowSOM)
         library(ConsensusClusterPlus)
         library(uwot)
         library(matrixStats)
         ## Stage 05 (figures)
         library(ComplexHeatmap)
         library(circlize)
         library(ggridges)
         library(ggrastr)
         library(ggrepel)
         library(RColorBrewer)
         ## Tidyverse last (namespace priority for dplyr verbs)
         library(tidyverse)
})

set.seed(2026)


###############################################################################
# SECTION 1 : GLOBAL CONFIGURATION  (the only block you routinely edit)
###############################################################################

## ---- 1.1  Which stages to run this session --------------------------------
#   c(1,2,3,4,5) = full pipeline from FCS.
#   c(5)         = figures only (requires Results/cell_data.rds etc.).
RUN_STAGES <- c(1, 2, 3, 4, 5)

## ---- 1.1b  Live plot display (VS Code Plots pane) --------------------------
#   TRUE  = each figure is also print()/draw()'n to the active graphics
#           device (shows up in VS Code as it's produced), on top of the
#           PDF/PNG/TIFF files always written to Figures/.
#   FALSE = disk-only (faster; no rendering overhead for dense point layers).
SHOW_PLOTS <- TRUE

## ---- 1.2  Working directory (raw exported FCS live here) ------------------
WORKDIR <- paste0(
         "C:/Users/LENOVO/OneDrive - Malawi-Liverpool Wellcome Research ",
         "Programme/Desktop/LAB WORK/Exported FCS Files for Mucosal Analysis/",
         "High Dimensionality reductionl-2026/DownSample"
)
if(dir.exists(WORKDIR)) setwd(WORKDIR)

## ---- 1.3  Output folders --------------------------------------------------
paths <- list(
         results   = "Results",
         fig_heat  = "Figures/Heatmaps",
         fig_umap  = "Figures/UMAP",
         fig_ridge = "Figures/Ridgeline",
         fig_pca   = "Figures/PCA",
         fig_qc    = "Figures/QC"
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

## ---- 1.4  Tunable analysis parameters -------------------------------------
params <- list(
         min_cells       = 10000,   # QC threshold: drop samples below this
         downsample_to   = 10000,   # cells retained per sample
         cofactor        = 150,     # arcsinh cofactor (fluorescence)
         som_xdim        = 10,      # FlowSOM grid x
         som_ydim        = 10,      # FlowSOM grid y  (-> 100 SOM nodes)
         som_rlen        = 20,      # FlowSOM training epochs
         consensus_maxK  = 30,      # ConsensusClusterPlus max metaclusters
         consensus_reps  = 100,
         meta_k          = 20,      # chosen metacluster resolution
         umap_neighbors  = 30,
         umap_min_dist   = 0.3,
         max_umap_pts    = 60000,   # cells rendered in UMAP scatter panels
         cap_per_cluster = 2000,    # cells / cluster for ridgeline KDE
         cap_per_tissue  = 5000,    # cells / tissue  for ridgeline KDE
         raster_dpi      = 600,
         export_dpi      = 600,
         pseudocount     = 1e-4     # CLR zero-handling for PCA
)

## ---- 1.5  Panel definition (single source of truth) -----------------------
#   fcs_colname : detector name in the FCS
#   antigen     : biological marker
#   role        : "cluster" markers drive FlowSOM/UMAP; "type" are kept but not
#                 clustered on; "none" are technical (CD45, LiveDead).
panel <- tibble::tribble(
         ~fcs_colname,                        ~antigen,   ~role,
         "FJComp-BV785-A-1",                  "CD3",      "cluster",
         "FJComp-APC-Cy7-A-1",                "CD4",      "cluster",
         "FJComp-AF700-A-1",                  "CD8",      "cluster",
         "FJComp-BV750-A-1",                  "CD19",     "cluster",
         "FJComp-BV421-A-1",                  "CD14",     "cluster",
         "FJComp-BV650-A-1",                  "CD16",     "cluster",
         "FJComp-PE-A-1",                     "CD66b",    "cluster",
         "FJComp-APC-A-1",                    "CD56",     "cluster",
         "FJComp-BV605-A-1",                  "CD69",     "cluster",
         "FJComp-PE-Dazzle594-A-1",           "CD103",    "cluster",
         "FJComp-FITC-A-1",                   "CD45RO",   "cluster",
         "FJComp-PE-Cy7-A-1",                 "CCR7",     "cluster",
         "FJComp-PerCP-Cy5.5-A-1",            "CD45",     "none",
         "FJComp-LiveDeadFixableAqua-A-1",    "LiveDead", "none"
)

# Scatter + time channels retained through QC (not clustered on).
scatter_channels <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","TIME-1")

# Convenience vectors derived from the panel.
cluster_channels <- panel$fcs_colname[panel$role == "cluster"]
markers          <- panel$antigen[panel$role == "cluster"]
marker_lookup    <- setNames(panel$antigen, panel$fcs_colname)  # channel->antigen
keep_channels    <- c(scatter_channels, panel$fcs_colname)

## ---- 1.6  Publication theme -----------------------------------------------
theme_pub <- theme_classic(base_size = 9) +
         theme(
                  panel.border     = element_blank(),
                  axis.line        = element_line(colour = "black", linewidth = 0.5),
                  axis.text        = element_text(colour = "black"),
                  legend.title     = element_text(face = "bold"),
                  legend.position  = "right",
                  plot.title       = element_text(face = "bold", hjust = 0.5),
                  strip.background = element_blank(),
                  strip.text       = element_text(face = "bold")
         )
theme_set(theme_pub)

## ---- 1.7  Palettes & colour ramp ------------------------------------------
cluster_cols <- colorRampPalette(c("#00BFC4", "#FF3E3E"))(20)

tissue_cols <- c(
         "Nasopharynx"         = "#1B9E77",
         "Inferior turbinate"  = "#7570B3",
         "Ectocervix"          = "#D95F02",
         "Endocervix"          = "#E7298A",
         "Unknown"             = "grey60"
)

# Inferno ramp for z-scored heatmaps (matches the FIG 2C feature-plot palette).
heat_ramp <- circlize::colorRamp2(
         seq(-2, 2, length.out = 100),
         viridisLite::inferno(100)
)

feature_colors <- c("#2166AC", "#B2182B")   # blue-red marker expression ramp


###############################################################################
# SECTION 2 : HELPER FUNCTIONS
###############################################################################

## ---- 2.1  Memory reporter -------------------------------------------------
mem_report <- function(step){
         cat("\n=====================================================\n")
         cat(step, "\n")
         print(gc())
         cat("=====================================================\n\n")
}

## ---- 2.2  Tissue assignment from filename ---------------------------------
assign_tissue <- function(fname){
         dplyr::case_when(
                  stringr::str_detect(fname, "Nasal Swab")         ~ "Nasopharynx",
                  stringr::str_detect(fname, "Nasal Scrape")       ~ "Inferior turbinate",
                  stringr::str_detect(fname, "Cervical Scrape")    ~ "Ectocervix",
                  stringr::str_detect(fname, "Cervical Cytobrush") ~ "Endocervix",
                  TRUE                                             ~ "Unknown"
         )
}

## ---- 2.3  Multi-format ggplot export --------------------------------------
save_fig <- function(plot, stem, width, height,
                     dpi = params$export_dpi, tiff = FALSE){
         if(isTRUE(SHOW_PLOTS)) print(plot)
         ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
         ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi)
         if(tiff)
                  ggsave(paste0(stem, ".tiff"), plot, width = width,
                         height = height, dpi = dpi, compression = "lzw")
         invisible(NULL)
}

## ---- 2.4  ComplexHeatmap export (draw() is not a ggplot) ------------------
save_heatmap <- function(ht, stem, width, height,
                         dpi = params$export_dpi, tiff = FALSE, ...){
         if(isTRUE(SHOW_PLOTS)) draw(ht, ...)
         pdf(paste0(stem, ".pdf"), width = width, height = height)
         draw(ht, ...); dev.off()
         png(paste0(stem, ".png"), width = width, height = height,
             units = "in", res = dpi)
         draw(ht, ...); dev.off()
         if(tiff){
                  tiff(paste0(stem, ".tiff"), width = width, height = height,
                       units = "in", res = dpi, compression = "lzw")
                  draw(ht, ...); dev.off()
         }
         invisible(NULL)
}

## ---- 2.5  Auto phenotype labels from z-scored marker MFI ------------------
#   Ranks markers within each cluster by z-scored median arcsinh MFI and
#   concatenates the top hits above z_min into a label, e.g. "CD3+CD4+CD69+".
#   Clusters with no marker above z_min are labelled "Unassigned".
label_clusters_by_mfi <- function(mat_z, top_n = 3, z_min = 0.5){
         labels <- apply(mat_z, 1, function(z_row){
                  hits <- sort(z_row[z_row > z_min], decreasing = TRUE)
                  if(length(hits) == 0) return("Unassigned")
                  paste0(names(hits)[seq_len(min(top_n, length(hits)))],
                         "+", collapse = "")
         })
         make.unique(labels, sep = "_")
}

mem_report("Configuration & helpers loaded")


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 01 : IMPORT, QC, CHANNEL CURATION, PANEL ANNOTATION                ##
##                                                                           ##
###############################################################################
###############################################################################

if(1 %in% RUN_STAGES){
         
         mem_report("START Stage 01")
         
         ## ---- 1.1  Locate FCS, drop PBMC (mucosal-only analysis) -----------
         fcs_files <- list.files(pattern = "\\.fcs$", full.names = TRUE,
                                 ignore.case = TRUE)
         cat(length(fcs_files), "FCS files detected.\n")
         
         fcs_files <- fcs_files[!grepl("PBMC", fcs_files, ignore.case = TRUE)]
         cat(length(fcs_files), "mucosal samples retained.\n")
         
         stopifnot(length(fcs_files) > 0)
         
         ## ---- 1.2  Read FlowSet --------------------------------------------
         fs <- read.flowSet(files = fcs_files, truncate_max_range = FALSE)
         mem_report("FlowSet imported")
         
         ## ---- 1.3  Event counts (pre-filter QC record) ---------------------
         event_counts <- fsApply(fs, nrow)
         cat("Event-count summary (pre-filter):\n"); print(summary(event_counts))
         
         ## ---- 1.4  Sample metadata -----------------------------------------
         metadata <- tibble(
                  file_name = basename(fcs_files),
                  sample_id = basename(fcs_files),
                  patient   = str_extract(basename(fcs_files), "DKP\\d+[A-Z]?"),
                  tissue    = assign_tissue(basename(fcs_files))
         )
         if(any(metadata$tissue == "Unknown"))
                  warning("Some files did not match a tissue pattern (tissue = 'Unknown').")
         
         ## ---- 1.5  Curate channels (keep ONE copy of each fluorochrome) ----
         missing_ch <- setdiff(keep_channels, colnames(fs))
         if(length(missing_ch))
                  stop("Channels absent from FCS: ", paste(missing_ch, collapse = ", "))
         fs <- fs[, keep_channels]
         mem_report("Channels curated")
         
         ## ---- 1.6  QC figure: cell yield per sample ------------------------
         qc_df <- tibble(sample_id = sampleNames(fs),
                         cells     = as.numeric(event_counts)) %>%
                  left_join(metadata, by = "sample_id")
         
         p_qc <- ggplot(qc_df, aes(reorder(sample_id, cells), cells,
                                   fill = tissue)) +
                  geom_col() +
                  geom_hline(yintercept = params$min_cells, linetype = 2,
                             colour = "red") +
                  scale_fill_manual(values = tissue_cols, name = "Tissue") +
                  coord_flip() +
                  labs(x = NULL, y = "Events (pre-downsampling)",
                       title = "Cell yield per mucosal sample")
         save_fig(p_qc, file.path(paths$fig_qc, "QC_cell_yield"),
                  width = 7, height = max(4, 0.18 * nrow(qc_df)))
         
         ## ---- 1.7  Checkpoint ----------------------------------------------
         saveRDS(fs,       file.path(paths$results, "fs_clean.rds"))
         saveRDS(metadata, file.path(paths$results, "metadata.rds"))
         saveRDS(panel,    file.path(paths$results, "panel.rds"))
         write.csv(qc_df, file.path(paths$results, "Cell_Counts_Before_Filtering.csv"),
                   row.names = FALSE)
         
         rm(event_counts, qc_df, p_qc, missing_ch); gc()
         mem_report("END Stage 01")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 02 : FILTER (>= min_cells) & DOWNSAMPLE                            ##
##                                                                           ##
###############################################################################
###############################################################################

if(2 %in% RUN_STAGES){
         
         mem_report("START Stage 02")
         
         if(!exists("fs")) fs <- readRDS(file.path(paths$results, "fs_clean.rds"))
         if(!exists("metadata"))
                  metadata <- readRDS(file.path(paths$results, "metadata.rds"))
         
         ## ---- 2.1  Keep samples with >= min_cells --------------------------
         event_counts <- fsApply(fs, nrow)
         keep <- event_counts >= params$min_cells
         cat(sum(keep),  "samples retained (>=", params$min_cells, "cells)\n")
         cat(sum(!keep), "samples removed\n")
         
         fs       <- fs[keep]
         metadata <- metadata %>% filter(sample_id %in% sampleNames(fs))
         rm(event_counts, keep); gc()
         
         ## ---- 2.2  Downsample to exactly downsample_to cells ---------------
         set.seed(2026)
         for(i in seq_len(length(fs))){
                  ff   <- fs[[i]]
                  expr <- exprs(ff)
                  idx  <- sample.int(nrow(expr), params$downsample_to, replace = FALSE)
                  exprs(ff) <- expr[idx, , drop = FALSE]
                  fs[[i]]   <- ff
                  rm(ff, expr, idx)
                  if(i %% 5 == 0) gc(verbose = FALSE)
         }
         mem_report("Downsampling complete")
         
         ## ---- 2.3  Verify & align metadata to FlowSet order ----------------
         new_counts <- fsApply(fs, nrow)
         stopifnot(all(new_counts == params$downsample_to))
         
         metadata <- metadata %>%
                  filter(sample_id %in% sampleNames(fs)) %>%
                  arrange(match(sample_id, sampleNames(fs)))
         stopifnot(all(metadata$sample_id == sampleNames(fs)))
         metadata$cells_downsampled <- params$downsample_to
         
         ## ---- 2.4  Checkpoint ----------------------------------------------
         saveRDS(fs,       file.path(paths$results, "fs_10k.rds"))
         saveRDS(metadata, file.path(paths$results, "metadata_10k.rds"))
         write.csv(metadata, file.path(paths$results, "Metadata_10k.csv"),
                   row.names = FALSE)
         
         rm(new_counts); gc()
         mem_report("END Stage 02")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 03 : TRANSFORM, FlowSOM + CONSENSUS METACLUSTERING, UMAP           ##
##                                                                           ##
###############################################################################
###############################################################################

if(3 %in% RUN_STAGES){
         
         mem_report("START Stage 03")
         
         if(!exists("fs")) fs <- readRDS(file.path(paths$results, "fs_10k.rds"))
         if(!exists("metadata"))
                  metadata <- readRDS(file.path(paths$results, "metadata_10k.rds"))
         
         ## ---- 3.1  Concatenate cluster channels + per-cell sample id -------
         expr_list <- vector("list", length(fs))
         sid_list  <- vector("list", length(fs))
         for(i in seq_along(fs)){
                  expr_list[[i]] <- exprs(fs[[i]])[, cluster_channels, drop = FALSE]
                  sid_list[[i]]  <- rep(sampleNames(fs)[i], nrow(expr_list[[i]]))
         }
         expr_all  <- do.call(rbind, expr_list)
         sample_id <- unlist(sid_list)
         rm(expr_list, sid_list); gc()
         cat("Concatenated matrix: ", nrow(expr_all), "cells x ",
             ncol(expr_all), "markers\n")
         
         ## ---- 3.2  Arcsinh transform ---------------------------------------
         expr_all <- asinh(expr_all / params$cofactor)
         gc()
         
         ## ---- 3.3  FlowSOM (100-node SOM) ----------------------------------
         set.seed(2026)
         fsom <- ReadInput(expr_all, transform = FALSE, scale = FALSE)
         fsom <- BuildSOM(fsom, colsToUse = seq_len(ncol(expr_all)),
                          xdim = params$som_xdim, ydim = params$som_ydim,
                          rlen = params$som_rlen)
         fsom <- BuildMST(fsom)
         saveRDS(fsom, file.path(paths$results, "fsom_100nodes.rds"))
         gc()
         
         ## ---- 3.4  Consensus metaclustering of SOM codes -------------------
         codes <- fsom$map$codes
         set.seed(2026)
         cc <- ConsensusClusterPlus(
                  t(codes),
                  maxK        = params$consensus_maxK,
                  reps        = params$consensus_reps,
                  pItem       = 0.9,
                  pFeature    = 1,
                  clusterAlg  = "hc",
                  distance    = "euclidean",
                  seed        = 2026,
                  plot        = "png",
                  title       = file.path(paths$results, "Consensus")
         )
         saveRDS(cc, file.path(paths$results, "ConsensusCluster.rds"))
         
         ## ---- 3.5  Map chosen k back to cells ------------------------------
         meta_assign  <- cc[[params$meta_k]]$consensusClass
         cell_cluster <- fsom$map$mapping[, 1]
         cell_meta    <- meta_assign[cell_cluster]
         saveRDS(cell_meta, file.path(paths$results, "Cell_Metaclusters.rds"))
         cat("Metacluster sizes (k =", params$meta_k, "):\n")
         print(table(cell_meta))
         
         ## ---- 3.6  UMAP embedding ------------------------------------------
         set.seed(2026)
         umap_coords <- uwot::umap(
                  expr_all,
                  n_neighbors = params$umap_neighbors,
                  min_dist    = params$umap_min_dist,
                  metric      = "euclidean",
                  verbose     = TRUE
         )
         saveRDS(umap_coords, file.path(paths$results, "UMAP.rds"))
         
         ## ---- 3.7  Persist transformed matrix pieces for Stage 04 ----------
         saveRDS(expr_all,  file.path(paths$results, "expr_all_asinh.rds"))
         saveRDS(sample_id, file.path(paths$results, "sample_id.rds"))
         
         rm(fsom, cc, codes, meta_assign, cell_cluster); gc()
         mem_report("END Stage 03")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 04 : MASTER CELL TABLE, CLUSTER SUMMARIES, ABUNDANCES              ##
##                                                                           ##
###############################################################################
###############################################################################

if(4 %in% RUN_STAGES){
         
         mem_report("START Stage 04")
         
         if(!exists("expr_all"))
                  expr_all  <- readRDS(file.path(paths$results, "expr_all_asinh.rds"))
         if(!exists("sample_id"))
                  sample_id <- readRDS(file.path(paths$results, "sample_id.rds"))
         if(!exists("cell_meta"))
                  cell_meta <- readRDS(file.path(paths$results, "Cell_Metaclusters.rds"))
         if(!exists("metadata"))
                  metadata  <- readRDS(file.path(paths$results, "metadata_10k.rds"))
         
         ## ---- 4.1  Master per-cell table (rename channels -> antigens) -----
         cell_data <- as.data.frame(expr_all)
         colnames(cell_data) <- marker_lookup[colnames(cell_data)]   # channel->antigen
         cell_data$Cluster <- factor(cell_meta)
         cell_data$Sample  <- sample_id
         cell_data <- left_join(cell_data, metadata,
                                by = c("Sample" = "sample_id"))
         saveRDS(cell_data, file.path(paths$results, "cell_data.rds"))
         rm(expr_all); gc()
         
         ## ---- 4.2  Cluster phenotype summary (median/mean/%pos) ------------
         cluster_summary <- cell_data %>%
                  group_by(Cluster) %>%
                  summarise(
                           Cells = n(),
                           across(all_of(markers),
                                  list(Median   = median,
                                       Mean     = mean,
                                       Positive = ~100 * mean(.x > 1)),
                                  .names = "{.col}_{.fn}"),
                           .groups = "drop"
                  )
         write.csv(cluster_summary,
                   file.path(paths$results, "Cluster_Summary.csv"),
                   row.names = FALSE)

         ## ---- 4.2b  Auto phenotype labels (z-scored median MFI per marker) -
         mat_median <- cluster_summary %>%
                  dplyr::select(Cluster, ends_with("_Median")) %>%
                  column_to_rownames("Cluster") %>%
                  rename_with(~ str_remove(.x, "_Median$")) %>%
                  as.matrix()
         mat_median_z <- scale(mat_median)
         mat_median_z[is.na(mat_median_z)] <- 0

         cluster_labels_df <- tibble(
                  Cluster   = factor(rownames(mat_median_z),
                                     levels = levels(cell_data$Cluster)),
                  Phenotype = label_clusters_by_mfi(mat_median_z)
         )
         write.csv(cluster_labels_df,
                   file.path(paths$results, "Cluster_Phenotype_Labels.csv"),
                   row.names = FALSE)

         cell_data <- cell_data %>% left_join(cluster_labels_df, by = "Cluster")
         saveRDS(cell_data, file.path(paths$results, "cell_data.rds"))

         rm(mat_median, mat_median_z); gc()

         ## ---- 4.3  Cluster sizes -------------------------------------------
         cluster_size <- cell_data %>%
                  count(Cluster, name = "Cells") %>%
                  mutate(Percent = round(100 * Cells / sum(Cells), 2)) %>%
                  arrange(desc(Cells))
         write.csv(cluster_size,
                   file.path(paths$results, "Cluster_Size.csv"),
                   row.names = FALSE)
         
         ## ---- 4.4  Per-sample cluster abundance (frequencies) --------------
         cluster_abundance <- cell_data %>%
                  count(Sample, Cluster, tissue, name = "n") %>%
                  group_by(Sample) %>%
                  mutate(Frequency = n / sum(n)) %>%
                  ungroup()
         saveRDS(cluster_abundance,
                 file.path(paths$results, "cluster_frequency.rds"))
         write.csv(cluster_abundance,
                   file.path(paths$results, "Cluster_Abundance.csv"),
                   row.names = FALSE)
         
         rm(cluster_summary, cluster_size); gc()
         mem_report("END Stage 04")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 05 : PUBLICATION FIGURES                                           ##
##     FIG 1  Clustered heatmap (cluster x marker)                            ##
##     FIG 2  UMAP panels (master; tissue facets; marker feature grid)        ##
##     FIG 3  Ridgeline densities (by cluster; by tissue)                     ##
##                                                                           ##
###############################################################################
###############################################################################

if(5 %in% RUN_STAGES){
         
         mem_report("START Stage 05")
         
         ## ---- Load figure inputs -------------------------------------------
         if(!exists("cell_data"))
                  cell_data <- readRDS(file.path(paths$results, "cell_data.rds"))
         if(!exists("umap_coords"))
                  umap_coords <- readRDS(file.path(paths$results, "UMAP.rds"))
         if(!exists("cluster_abundance"))
                  cluster_abundance <- readRDS(file.path(paths$results,
                                                         "cluster_frequency.rds"))
         
         # Defensive rename guard (in case cell_data still holds channel names).
         if(!all(markers %in% colnames(cell_data))){
                  hit <- match(names(marker_lookup), colnames(cell_data))
                  ok  <- !is.na(hit)
                  colnames(cell_data)[hit[ok]] <- unname(marker_lookup[ok])
         }
         stopifnot(all(markers %in% colnames(cell_data)))
         cell_data$Cluster <- factor(cell_data$Cluster)

         # Defensive guard: attach auto phenotype labels if this cell_data.rds
         # predates the Stage 04 MFI-labelling step.
         if(!"Phenotype" %in% colnames(cell_data)){
                  cluster_labels_df <- read.csv(
                           file.path(paths$results, "Cluster_Phenotype_Labels.csv"),
                           stringsAsFactors = FALSE) %>%
                           mutate(Cluster = factor(Cluster, levels = levels(cell_data$Cluster)))
                  cell_data <- cell_data %>% left_join(cluster_labels_df, by = "Cluster")
         }

         clus_levels  <- levels(cell_data$Cluster)
         pal_clusters <- setNames(cluster_cols[seq_along(clus_levels)], clus_levels)

         # "MC3: CD3+CD4+CD69+" style combined label, keyed by cluster id, for
         # UMAP text/legend and heatmap row names.
         pheno_by_cluster <- cell_data %>%
                  distinct(Cluster, Phenotype) %>%
                  arrange(Cluster) %>%
                  mutate(Label = paste0("MC", Cluster, ": ", Phenotype))
         combo_labels <- setNames(pheno_by_cluster$Label, pheno_by_cluster$Cluster)

         mem_report("Figure inputs ready")
         
         #####################################################################
         # FIG 1 : CLUSTERED HEATMAPS
         #####################################################################
         
         ## ---- 1B  cluster x marker phenotype key ---------------------------
         mat_cm <- cell_data %>%
                  group_by(Cluster) %>%
                  summarise(across(all_of(markers), median), .groups = "drop") %>%
                  column_to_rownames("Cluster") %>% as.matrix()
         rownames(mat_cm) <- combo_labels[rownames(mat_cm)]
         mat_cm_z <- scale(mat_cm); mat_cm_z[is.na(mat_cm_z)] <- 0
         
         ht_cm <- Heatmap(
                  mat_cm_z, name = "z-score", col = heat_ramp,
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_gp = grid::gpar(fontsize = 8),
                  column_names_gp = grid::gpar(fontsize = 8),
                  row_dend_width = grid::unit(12, "mm"),
                  column_dend_height = grid::unit(12, "mm"),
                  border = FALSE
         )
         save_heatmap(ht_cm, file.path(paths$fig_heat, "FIG1B_Cluster_Marker_Heatmap"),
                      width = 6, height = 7)
         
         rm(mat_cm, mat_cm_z, ht_cm)
         gc(); mem_report("FIG 1 complete")

         #####################################################################
         # FIG 2 : UMAP PANELS
         #####################################################################
         
         umap_df <- tibble(UMAP1 = umap_coords[,1], UMAP2 = umap_coords[,2],
                           Cluster = cell_data$Cluster, tissue = cell_data$tissue)
         
         keep_idx <- if(nrow(umap_df) > params$max_umap_pts)
                  sample.int(nrow(umap_df), params$max_umap_pts) else seq_len(nrow(umap_df))
         umap_plot_df <- umap_df[keep_idx, ]
         
         cluster_labels <- umap_plot_df %>%
                  group_by(Cluster) %>%
                  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
                            .groups = "drop") %>%
                  mutate(Phenotype = pheno_by_cluster$Phenotype[
                           match(Cluster, pheno_by_cluster$Cluster)])

         ## ---- 2A  master UMAP -----------------------------------------------
         #   Cluster number anchors each point cloud; the auto phenotype label
         #   (top z-scored markers by median arcsinh MFI, Cluster_Phenotype_
         #   Labels.csv) is repelled alongside it, and the same "MC#: markers"
         #   string drives the legend so colour keys read as cell types.
         p_master <- ggplot(umap_plot_df, aes(UMAP1, UMAP2, colour = Cluster)) +
                  rasterise(geom_point(size = 0.15, alpha = 0.7),
                            dpi = params$raster_dpi) +
                  geom_text(data = cluster_labels,
                            aes(UMAP1, UMAP2, label = Cluster),
                            colour = "black", fontface = "bold",
                            size = 3.8, show.legend = FALSE) +
                  ggrepel::geom_text_repel(
                           data = cluster_labels,
                           aes(UMAP1, UMAP2, label = Phenotype),
                           colour = "black", size = 2.6, fontface = "italic",
                           bg.color = "white", bg.r = 0.12,
                           max.overlaps = Inf, show.legend = FALSE) +
                  scale_colour_viridis_d(option = "magma",
                                        name = "Metacluster",
                                        labels = combo_labels[clus_levels],
                                        guide = guide_legend(override.aes = list(size = 4),
                                                             title.position = "top",
                                                             ncol = 1)) +
                  coord_equal() +
                  labs(title = "Mucosal immune landscape (FlowSOM metaclusters)") +
                  theme(legend.position = "right",
                        legend.title = element_text(face = "bold"),
                        legend.text = element_text(size = 6))
         save_fig(p_master, file.path(paths$fig_umap, "FIG2A_UMAP_master"),
                  width = 8.5, height = 5.5)
 
         ## ---- 2B  faceted by tissue ----------------------------------------
         backdrop <- umap_plot_df %>% dplyr::select(UMAP1, UMAP2)
         p_facets <- ggplot() +
                  rasterise(geom_point(data = backdrop, aes(UMAP1, UMAP2),
                                       colour = "grey88", size = 0.10),
                            dpi = params$raster_dpi) +
                  rasterise(geom_point(data = umap_plot_df,
                                       aes(UMAP1, UMAP2, colour = Cluster),
                                       size = 0.14, alpha = 0.7),
                            dpi = params$raster_dpi) +
                  scale_colour_viridis_d(option = "magma", guide = "none") +
                  facet_wrap(~ tissue, ncol = 2) + coord_equal() +
                  theme(strip.text = element_text(face = "bold", size = 10))
         save_fig(p_facets, file.path(paths$fig_umap, "FIG2B_UMAP_by_tissue"),
                  width = 8, height = 7)
 
         ## ---- 2C  per-marker feature grid ------------------------------------
         #   Inferno scale: black = no/negative expression (arcsinh <= 0),
         #   through purple/orange to yellow = strongly positive. Expression is
         #   scaled per marker to that marker's own 1st-99th percentile range
         #   (so dim and bright markers are both legible side by side, and a
         #   few extreme events can't wash out the colour scale). Points are
         #   drawn in ascending expression order so positive cells layer on
         #   top of the negative background.
         feat_long <- bind_cols(umap_plot_df[, c("UMAP1", "UMAP2")],
                                cell_data[keep_idx, markers]) %>%
                  pivot_longer(all_of(markers), names_to = "Marker",
                               values_to = "Expression") %>%
                  mutate(Marker = factor(Marker, levels = markers)) %>%
                  group_by(Marker) %>%
                  mutate(
                           hi     = quantile(Expression, 0.99),
                           Clip   = pmin(pmax(Expression, 0), hi),
                           Scaled = if_else(hi > 0, Clip / hi, 0)
                  ) %>%
                  ungroup() %>%
                  arrange(Marker, Scaled)

         p_features <- ggplot(feat_long, aes(UMAP1, UMAP2, colour = Scaled)) +
                  rasterise(geom_point(size = 0.14), dpi = params$raster_dpi) +
                  scale_colour_viridis_c(option = "inferno",
                                        limits  = c(0, 1),
                                        breaks  = c(0, 1),
                                        labels  = c("Negative", "Positive"),
                                        name    = NULL) +
                  coord_equal() +
                  facet_wrap(~ Marker, ncol = 4) +
                  theme(legend.position = "bottom",
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(face = "bold"))
         save_fig(p_features, file.path(paths$fig_umap, "FIG2C_Feature_grid"),
                  width = 10, height = 8)

         rm(umap_df, umap_plot_df, backdrop, feat_long,
            p_master, p_facets, p_features, keep_idx)
         gc(); mem_report("FIG 2 complete")
         
         #####################################################################
         # FIG 3 : RIDGELINE DENSITIES
         #####################################################################
         
         ## ---- 3A  by metacluster -------------------------------------------
         ridge_cluster <- cell_data %>%
                  dplyr::select(Cluster, all_of(markers)) %>%
                  group_by(Cluster) %>%
                  slice_sample(n = params$cap_per_cluster) %>% ungroup() %>%
                  pivot_longer(all_of(markers), names_to = "Marker",
                               values_to = "Expression") %>%
                  mutate(Marker = factor(Marker, levels = markers))
         
         p_ridge_cl <- ggplot(ridge_cluster,
                              aes(Expression, Cluster, fill = Cluster)) +
                  geom_density_ridges(scale = 2.2, linewidth = 0.25,
                                      rel_min_height = 0.01, alpha = 0.9) +
                  scale_fill_manual(values = pal_clusters, guide = "none") +
                  facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
                  labs(x = "Arcsinh expression", y = "Metacluster") +
                  theme_ridges(font_size = 9, grid = TRUE) +
                  theme(strip.text = element_text(face = "bold"),
                        axis.title = element_text(face = "bold"))
         save_fig(p_ridge_cl, file.path(paths$fig_ridge, "FIG3A_Ridgeline_by_cluster"),
                  width = 11, height = 9)
         rm(ridge_cluster, p_ridge_cl); gc()
         
         ## ---- 3B  by tissue ------------------------------------------------
         ridge_tissue <- cell_data %>%
                  dplyr::select(tissue, all_of(markers)) %>%
                  group_by(tissue) %>%
                  slice_sample(n = params$cap_per_tissue) %>% ungroup() %>%
                  pivot_longer(all_of(markers), names_to = "Marker",
                               values_to = "Expression") %>%
                  mutate(Marker = factor(Marker, levels = markers))
         
         p_ridge_ti <- ggplot(ridge_tissue,
                              aes(Expression, tissue, fill = tissue)) +
                  geom_density_ridges(scale = 1.8, linewidth = 0.3,
                                      rel_min_height = 0.01, alpha = 0.9) +
                  scale_fill_manual(values = tissue_cols, guide = "none") +
                  facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
                  labs(x = "Arcsinh expression", y = "Tissue") +
                  theme_ridges(font_size = 9, grid = TRUE) +
                  theme(strip.text = element_text(face = "bold"),
                        axis.title = element_text(face = "bold"))
         save_fig(p_ridge_ti, file.path(paths$fig_ridge, "FIG3B_Ridgeline_by_tissue"),
                  width = 11, height = 7)
         rm(ridge_tissue, p_ridge_ti); gc(); mem_report("FIG 3 complete")
}

# NOTE: FIG 1A (marker x sample heatmap), FIG 4 (cluster-frequency PCA), and
# Stage 05b (feature-driven PCA) have been removed per request.


###############################################################################
# SESSION RECORD
###############################################################################

writeLines(capture.output(sessionInfo()),
           file.path(paths$results, "sessionInfo_MasterPipeline.txt"))

cat("\n---------------------------------------------------------------\n")
cat("MASTER PIPELINE finished.  Stages run:", paste(RUN_STAGES, collapse = ", "), "\n")
cat("Checkpoints & tables : ./Results/\n")
cat("Figures              : ./Figures/{QC,Heatmaps,UMAP,Ridgeline,PCA}/\n")
cat("---------------------------------------------------------------\n")
mem_report("END MASTER PIPELINE")
