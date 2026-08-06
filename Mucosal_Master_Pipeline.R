###############################################################################
#                                                                             #
#   MUCOSAL HIGH-DIMENSIONAL IMMUNE PROFILING - MASTER PIPELINE                #
#   Spectral flow cytometry (Sony ID7000) | Mucosal compartments              #
#                                                                             #
#   Author : Robert Nyirenda (MLW, Infection & Immunity Group)                #
#                                                                             #
#   This script goes from raw FCS files all the way to publication figures.   #
#   It is written to be READABLE: each step has a comment saying WHY it is     #
#   done, and the code avoids dense one-liners so you can follow and edit it.  #
#                                                                             #
#   STAGES (run in order the first time):                                     #
#     01  Import FCS, quality control, choose channels, build metadata        #
#     02  Remove low-cell samples, downsample everyone to the same count      #
#     03  Transform data, cluster (FlowSOM), embed (UMAP)                      #
#     03b Quick batch-check UMAPs (colour by sample, then by donor)           #
#     04  Build the master per-cell table + summary tables                    #
#     05  Publication figures (heatmaps, UMAP, ridgelines, PCA)               #
#     05b Feature-driven PCA (samples described by marker medians)            #
#                                                                             #
#   HOW IT SAVES MEMORY (important on an 8 GB laptop):                         #
#     - Each stage saves its result to an .rds file (a "checkpoint").          #
#     - Later stages reload those files, so you never redo slow steps.         #
#     - We delete big objects with rm() and call gc() to free RAM.             #
#     - Figures use only a subset of cells so plotting stays light.            #
#                                                                             #
###############################################################################


###############################################################################
# SECTION 0 : START FRESH & LOAD PACKAGES
###############################################################################

# Clear the workspace so nothing from a previous run interferes.
rm(list = ls())
graphics.off()
gc()

# ---- ONE-TIME INSTALL (uncomment and run once on a fresh machine) ---------
# CRAN packages:
#   install.packages(c("tidyverse","uwot","matrixStats","circlize",
#                      "ggridges","patchwork","ggrastr","ggrepel",
#                      "RColorBrewer","BiocManager"))
# Bioconductor packages:
#   BiocManager::install(c("flowCore","FlowSOM","ConsensusClusterPlus",
#                          "ComplexHeatmap"))
# ---------------------------------------------------------------------------

# Load every package the whole script needs, once, at the top.
# (suppressPackageStartupMessages just hides the noisy loading text.)
suppressPackageStartupMessages({
         library(flowCore)             # read/handle FCS flow cytometry files
         library(FlowSOM)              # clustering of cytometry data
         library(ConsensusClusterPlus) # decide/stabilise the number of metaclusters
         library(uwot)                 # UMAP embedding
         library(matrixStats)          # fast column variance (used in PCA)
         library(ComplexHeatmap)       # heatmaps
         library(circlize)             # colour ramp for heatmaps
         library(ggridges)             # ridgeline (stacked density) plots
         library(patchwork)            # combine ggplots into panels
         library(ggrastr)              # rasterise dense points (small, fast PDFs)
         library(ggrepel)              # non-overlapping text labels
         library(RColorBrewer)         # colour palettes
         library(tidyverse)            # dplyr, ggplot2, tidyr, etc. (load LAST)
})

# Fix the random seed so downsampling, clustering and UMAP are reproducible.
set.seed(2026)


###############################################################################
# SECTION 1 : SETTINGS YOU MIGHT CHANGE
#   Everything you would normally adjust lives here, as plain variables.
###############################################################################

## ---- 1.1  Which stages to run this session --------------------------------
# Put the stage numbers you want to run in this vector.
#   Full run from raw FCS : c("1", "2", "3", "3b", "4", "5", "5b")
#   Only redo the figures : c("5")
# NOTE: we use text ("1", "2", ...) so that "3b" fits in naturally.
run_stages <- c("1", "2", "3", "3b", "4", "5", "5b")

## ---- 1.2  Where the raw FCS files live ------------------------------------
work_dir <- paste0(
         "C:/Users/LENOVO/OneDrive - Malawi-Liverpool Wellcome Research ",
         "Programme/Desktop/LAB WORK/Exported FCS Files for Mucosal Analysis/",
         "High Dimensionality reductionl-2026/DownSample"
)
if (dir.exists(work_dir)) setwd(work_dir)

## ---- 1.3  Output folders (created if they do not exist) -------------------
dir_results <- "Results"
dir_heat    <- "Figures/Heatmaps"
dir_umap    <- "Figures/UMAP"
dir_ridge   <- "Figures/Ridgeline"
dir_pca     <- "Figures/PCA"
dir_qc      <- "Figures/QC"

# Make each folder. recursive = TRUE also creates the parent "Figures" folder.
for (folder in c(dir_results, dir_heat, dir_umap, dir_ridge, dir_pca, dir_qc)) {
         dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

## ---- 1.4  Analysis settings (numbers you may tune) ------------------------
min_cells       <- 5000    # drop any sample with fewer cells than this
downsample_to   <- 3000    # keep exactly this many cells per sample
cofactor        <- 150     # arcsinh cofactor for fluorescence data
som_xdim        <- 10      # FlowSOM grid width
som_ydim        <- 10      # FlowSOM grid height (10 x 10 = 100 SOM nodes)
som_rlen        <- 20      # FlowSOM training rounds
consensus_maxk  <- 30      # try up to this many metaclusters
consensus_reps  <- 100     # resampling repeats for stability
meta_k          <- 15      # metaclusters used. Consensus matrices were cleanest
#   at k=11 (last k with uniformly crisp blocks); k=15
#   is chosen here for finer phenotypic resolution.
#   Treat clusters beyond ~11 as exploratory
#   sub-structure and rely on the k-independent
#   lineage roll-up (Stage 06) for formal abundance
#   inference. A k=11 sensitivity check is advised.
umap_neighbors  <- 30      # UMAP: neighbourhood size
umap_min_dist   <- 0.3     # UMAP: how tightly points pack
max_umap_points <- 60000   # cap on cells drawn in a UMAP figure
cap_per_cluster <- 2000    # cells per cluster used for ridgeline density
cap_per_tissue  <- 5000    # cells per tissue used for ridgeline density
raster_dpi      <- 600     # resolution when rasterising points
export_dpi      <- 600     # resolution for saved PNG/TIFF
pseudocount     <- 1e-4    # small value so log() never sees a zero (PCA)

## ---- 1.5  The antibody panel (single place that defines the markers) ------
# One row per detector. "role" tells the script how to use each channel:
#   "cluster" -> used to build clusters and the UMAP
#   "none"    -> kept in the file but NOT used for clustering (CD45, LiveDead)
panel <- tibble::tribble(
         ~fcs_colname,                     ~antigen,   ~role,
         "FJComp-BV785-A-1",               "CD3",      "cluster",
         "FJComp-APC-Cy7-A-1",             "CD4",      "cluster",
         "FJComp-AF700-A-1",               "CD8",      "cluster",
         "FJComp-BV750-A-1",               "CD19",     "cluster",
         "FJComp-BV421-A-1",               "CD14",     "cluster",
         "FJComp-BV650-A-1",               "CD16",     "cluster",
         "FJComp-PE-A-1",                  "CD66b",    "cluster",
         "FJComp-APC-A-1",                 "CD56",     "cluster",
         "FJComp-BV605-A-1",               "CD69",     "cluster",
         "FJComp-PE-Dazzle594-A-1",        "CD103",    "cluster",
         "FJComp-FITC-A-1",                "CD45RO",   "cluster",
         "FJComp-PE-Cy7-A-1",              "CCR7",     "cluster",
         "FJComp-PerCP-Cy5.5-A-1",         "CD45",     "none",
         "FJComp-LiveDeadFixableAqua-A-1", "LiveDead", "none"
)

# Scatter and time channels we keep through QC but never cluster on.
scatter_channels <- c("FSC-A", "FSC-H", "FSC-W",
                      "SSC-A", "SSC-H", "SSC-W", "TIME-1")

# Pull a few handy vectors out of the panel table.
# cluster_channels : detector names used for clustering
# markers          : the human-readable marker names (CD3, CD4, ...)
# keep_channels    : everything we keep when we read the FCS files
cluster_channels <- panel$fcs_colname[panel$role == "cluster"]
markers          <- panel$antigen[panel$role == "cluster"]
keep_channels    <- c(scatter_channels, panel$fcs_colname)

# A named lookup so we can turn detector names into marker names later.
# e.g. marker_lookup["FJComp-BV785-A-1"] gives "CD3".
marker_lookup <- panel$antigen
names(marker_lookup) <- panel$fcs_colname

## ---- 1.6  A clean, consistent plot theme ----------------------------------
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
theme_set(theme_pub)   # apply this theme to every ggplot automatically

## ---- 1.7  Colours ---------------------------------------------------------
# A colour-blind-friendly palette for clusters (extend if you have >20).
cluster_cols <- c(
         "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
         "#D55E00", "#F0E442", "#000000", "#882255", "#44AA99",
         "#117733", "#999933", "#88CCEE", "#AA4499", "#DDCC77",
         "#332288", "#661100", "#6699CC", "#AA4466", "#228833"
)

# Fixed colour for each tissue, so every figure uses the same colours.
# Fixed colour for each tissue, so every figure uses the same colours.
# Ordered anatomically: nasal (mucosa -> turbinate) then cervical (ecto -> endo).
tissue_cols <- c(
         "Nasal Mucosa"       = "#1B9E77",
         "Inferior Turbinate" = "#7570B3",
         "Ectocervix"         = "#D95F02",
         "Endocervix"         = "#E7298A",
         "Unknown"            = "grey60"
)

# Expression colour scale for z-scored heatmaps.
# Yellow-black-blue diverging map (the classic microarray/TreeView look used
# in the reference figure): HIGH = bright yellow, MID = black, LOW = blue.
# Symmetric around 0; saturates at +/-4 z-scores (matches the reference legend).
heat_ramp <- circlize::colorRamp2(
         c(-4, -2, 0, 2, 4),
         c("#0099FF", "#0033CC", "#000000", "#CCCC00", "#FFFF00")
)

# viridis "plasma" option for continuous marker UMAPs.
feature_option <- "C"


###############################################################################
# SECTION 2 : SMALL HELPER FUNCTIONS
#   Little reusable tools so we do not repeat ourselves.
###############################################################################

## ---- 2.1  Print a memory / progress report --------------------------------
mem_report <- function(step) {
         cat("\n=====================================================\n")
         cat(step, "\n")
         print(gc())   # gc() both frees memory and prints how much is in use
         cat("=====================================================\n\n")
}

## ---- 2.2  Work out the tissue from a file name ----------------------------
# Returns one tissue label per file name. If nothing matches -> "Unknown".
assign_tissue <- function(file_names) {
         tissue <- rep("Unknown", length(file_names))
         # LEFT of grepl() = the string in the FCS FILE NAME (unchanged).
         # RIGHT of <-     = the anatomical label we assign (renamed).
         # If your file names ever change, update the grepl() patterns, not the labels.
         tissue[grepl("Nasal Swab",         file_names)] <- "Nasal Mucosa"       # was Nasal Swab
         tissue[grepl("Nasal Scrape",       file_names)] <- "Inferior Turbinate" # was Nasal Scrape
         tissue[grepl("Cervical Scrape",    file_names)] <- "Ectocervix"         # was Cervical Scrape
         tissue[grepl("Cervical Cytobrush", file_names)] <- "Endocervix"         # was Cervical Cytobrush
         tissue
}

## ---- 2.3  Save a ggplot as a PNG -----------------------------------------
# "stem" is the file path WITHOUT extension; we add .png ourselves.
# PNG only (no PDF/TIFF) - keeps outputs light and consistent.
save_fig <- function(plot, stem, width, height) {
         ggsave(paste0(stem, ".png"), plot, width = width, height = height,
                dpi = export_dpi)
}

## ---- 2.4  Save a ComplexHeatmap as a PNG (it needs draw(), not ggsave) -----
# PNG only (no PDF/TIFF).
save_heatmap <- function(ht, stem, width, height, ...) {
         png(paste0(stem, ".png"), width = width, height = height,
             units = "in", res = export_dpi)
         draw(ht, ...)
         dev.off()
}

## ---- 2.5  Z-score each ROW of a matrix ------------------------------------
# Z-scoring puts every marker on the same scale so colours are comparable.
# If a row is constant, scale() gives NaN; we set those back to 0.
zscore_rows <- function(m) {
         z <- t(scale(t(m)))
         z[is.na(z)] <- 0
         z
}

## ---- 2.6  Sample up to N rows PER GROUP (safe if a group is smaller) -------
# slice_sample(n = ...) errors when a group has fewer than n rows, and
# dplyr::n() cannot be used inside its n = argument. This helper samples
# min(n, group size) rows from each group of an already-grouped tibble.
sample_up_to <- function(df_grouped, n) {
         df_grouped %>%
                  group_modify(function(g, key) {
                           take <- min(n, nrow(g))
                           g[sample.int(nrow(g), take), , drop = FALSE]
                  }) %>%
                  ungroup()
}

mem_report("Settings and helper functions loaded")


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 01 : IMPORT, QUALITY CONTROL, CHOOSE CHANNELS, BUILD METADATA      ##
##                                                                           ##
###############################################################################
###############################################################################

if ("1" %in% run_stages) {
         
         mem_report("START Stage 01")
         
         ## ---- Step 1: find the FCS files, remove blood (PBMC) samples ------------
         # We only want mucosal samples in this analysis, so drop anything with
         # "PBMC" in the file name.
         fcs_files <- list.files(pattern = "\\.fcs$", full.names = TRUE,
                                 ignore.case = TRUE)
         cat(length(fcs_files), "FCS files found.\n")
         
         is_pbmc   <- grepl("PBMC", fcs_files, ignore.case = TRUE)
         fcs_files <- fcs_files[!is_pbmc]
         cat(length(fcs_files), "mucosal samples kept (PBMC removed).\n")
         
         # Stop early with a clear message if no files were found.
         if (length(fcs_files) == 0) {
                  stop("No mucosal FCS files found. Check the working directory in Section 1.")
         }
         
         ## ---- Step 2: read all files into one flowSet ----------------------------
         # truncate_max_range = FALSE keeps the full signal range (no clipping).
         fs <- read.flowSet(files = fcs_files, truncate_max_range = FALSE)
         mem_report("FCS files read")
         
         ## ---- Step 3: count cells per sample (for the QC record) -----------------
         # fsApply runs a function on each sample; here we count rows (= cells).
         event_counts <- fsApply(fs, nrow)
         cat("Cell counts per sample (before filtering):\n")
         print(summary(event_counts))
         
         ## ---- Step 4: build a small metadata table -------------------------------
         # We pull the donor ID (e.g. DKP104) and the tissue out of each file name.
         file_names <- basename(fcs_files)
         metadata <- tibble(
                  file_name = file_names,
                  sample_id = file_names,
                  patient   = str_extract(file_names, "DKP\\d+[A-Z]?"),
                  tissue    = assign_tissue(file_names)
         )
         
         # Warn if any file did not match a tissue pattern (typo in file name?).
         if (any(metadata$tissue == "Unknown")) {
                  warning("Some files did not match a tissue name (tissue = 'Unknown'). ",
                          "Check those file names.")
         }
         
         ## ---- Step 5: keep only the channels we need -----------------------------
         # First check every wanted channel is actually present; stop if not.
         channels_present <- colnames(fs)
         channels_missing <- setdiff(keep_channels, channels_present)
         if (length(channels_missing) > 0) {
                  stop("These channels are missing from the FCS files: ",
                       paste(channels_missing, collapse = ", "))
         }
         fs <- fs[, keep_channels]
         mem_report("Unneeded channels dropped")
         
         ## ---- Step 6: QC figure - how many cells did each sample give? -----------
         qc_df <- tibble(
                  sample_id = sampleNames(fs),
                  cells     = as.numeric(event_counts)
         )
         qc_df <- left_join(qc_df, metadata, by = "sample_id")
         
         # A horizontal bar chart, one bar per sample, with the QC cut-off in red.
         p_qc <- ggplot(qc_df, aes(x = reorder(sample_id, cells), y = cells,
                                   fill = tissue)) +
                  geom_col() +
                  geom_hline(yintercept = min_cells, linetype = 2, colour = "red") +
                  scale_fill_manual(values = tissue_cols, name = "Tissue") +
                  coord_flip() +
                  labs(x = NULL, y = "Cells (before downsampling)",
                       title = "Cell yield per mucosal sample")
         
         # Height grows with the number of samples so bars do not get squashed.
         save_fig(p_qc, file.path(dir_qc, "QC_cell_yield"),
                  width = 7, height = max(4, 0.18 * nrow(qc_df)))
         
         ## ---- Step 7: save checkpoints -------------------------------------------
         saveRDS(fs,       file.path(dir_results, "fs_clean.rds"))
         saveRDS(metadata, file.path(dir_results, "metadata.rds"))
         saveRDS(panel,    file.path(dir_results, "panel.rds"))
         write.csv(qc_df, file.path(dir_results, "Cell_Counts_Before_Filtering.csv"),
                   row.names = FALSE)
         
         # Free memory we no longer need.
         rm(event_counts, qc_df, p_qc, is_pbmc, file_names,
            channels_present, channels_missing)
         gc()
         mem_report("END Stage 01")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 02 : REMOVE LOW-CELL SAMPLES, THEN DOWNSAMPLE EVERYONE EQUALLY     ##
##                                                                           ##
###############################################################################
###############################################################################

if ("2" %in% run_stages) {
         
         mem_report("START Stage 02")
         
         # Reload the checkpoints if they are not already in memory.
         if (!exists("fs"))       fs       <- readRDS(file.path(dir_results, "fs_clean.rds"))
         if (!exists("metadata")) metadata <- readRDS(file.path(dir_results, "metadata.rds"))
         
         ## ---- Step 1: keep samples with at least min_cells cells -----------------
         event_counts <- fsApply(fs, nrow)
         keep <- event_counts >= min_cells
         cat(sum(keep),  "samples kept (>=", min_cells, "cells).\n")
         cat(sum(!keep), "samples removed (too few cells).\n")
         
         fs       <- fs[keep]
         metadata <- metadata[metadata$sample_id %in% sampleNames(fs), ]
         
         rm(event_counts, keep)
         gc()
         
         ## ---- Step 2: downsample every sample to the SAME number of cells --------
         # Reason: samples with more cells would otherwise dominate the clustering.
         # We loop over samples, randomly pick 'downsample_to' cells from each,
         # and clean up as we go to stay within 8 GB of RAM.
         set.seed(2026)
         for (i in seq_len(length(fs))) {
                  
                  one_sample <- fs[[i]]                    # a single flowFrame
                  expr       <- exprs(one_sample)          # its numeric matrix (cells x channels)
                  
                  # Randomly choose which rows (cells) to keep.
                  chosen_rows <- sample.int(n = nrow(expr), size = downsample_to,
                                            replace = FALSE)
                  exprs(one_sample) <- expr[chosen_rows, , drop = FALSE]
                  
                  fs[[i]] <- one_sample                    # put the trimmed sample back
                  
                  rm(one_sample, expr, chosen_rows)
                  if (i %% 5 == 0) gc(verbose = FALSE)     # tidy memory every 5 samples
         }
         mem_report("Downsampling done")
         
         ## ---- Step 3: check it worked, and line metadata up with the samples -----
         new_counts <- fsApply(fs, nrow)
         # Every sample must now have exactly downsample_to cells; stop if not.
         stopifnot(all(new_counts == downsample_to))
         
         # Reorder metadata rows to match the order of samples in fs.
         metadata <- metadata[metadata$sample_id %in% sampleNames(fs), ]
         metadata <- metadata[match(sampleNames(fs), metadata$sample_id), ]
         stopifnot(all(metadata$sample_id == sampleNames(fs)))
         
         metadata$cells_downsampled <- downsample_to
         
         ## ---- Step 4: save checkpoints -------------------------------------------
         saveRDS(fs,       file.path(dir_results, "fs_downsampled.rds"))
         saveRDS(metadata, file.path(dir_results, "metadata_downsampled.rds"))
         write.csv(metadata, file.path(dir_results, "Metadata_Downsampled.csv"),
                   row.names = FALSE)
         
         rm(new_counts)
         gc()
         mem_report("END Stage 02")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 03 : TRANSFORM, CLUSTER (FlowSOM), EMBED (UMAP)                    ##
##                                                                           ##
###############################################################################
###############################################################################

if ("3" %in% run_stages) {
         
         mem_report("START Stage 03")
         
         if (!exists("fs"))       fs       <- readRDS(file.path(dir_results, "fs_downsampled.rds"))
         if (!exists("metadata")) metadata <- readRDS(file.path(dir_results, "metadata_downsampled.rds"))
         
         ## ---- Step 1: stack all samples into one big matrix ----------------------
         # We build two things in a loop:
         #   expr_list  : the marker values for each sample
         #   sid_list   : which sample each cell came from (same length)
         n_samples <- length(fs)
         expr_list <- vector("list", n_samples)
         sid_list  <- vector("list", n_samples)
         
         for (i in seq_len(n_samples)) {
                  this_expr <- exprs(fs[[i]])[, cluster_channels, drop = FALSE]
                  expr_list[[i]] <- this_expr
                  sid_list[[i]]  <- rep(sampleNames(fs)[i], nrow(this_expr))
         }
         
         # rbind stacks the per-sample matrices; unlist flattens the sample IDs.
         expr_all  <- do.call(rbind, expr_list)
         sample_id <- unlist(sid_list)
         
         rm(expr_list, sid_list)
         gc()
         cat("Combined matrix:", nrow(expr_all), "cells x",
             ncol(expr_all), "markers.\n")
         
         ## ---- Step 2: arcsinh transform ------------------------------------------
         # Fluorescence values span a huge range. asinh(x / cofactor) compresses
         # them into a scale where clusters separate nicely (standard for cytometry).
         expr_all <- asinh(expr_all / cofactor)
         gc()
         
         ## ---- Step 3: FlowSOM - group cells into 100 small nodes -----------------
         # ReadInput wraps the matrix; BuildSOM does the self-organising map;
         # BuildMST connects the nodes (used later for plotting).
         set.seed(2026)
         fsom <- ReadInput(expr_all, transform = FALSE, scale = FALSE)
         fsom <- BuildSOM(fsom,
                          colsToUse = seq_len(ncol(expr_all)),
                          xdim = som_xdim, ydim = som_ydim, rlen = som_rlen)
         fsom <- BuildMST(fsom)
         saveRDS(fsom, file.path(dir_results, "fsom_100nodes.rds"))
         gc()
         
         ## ---- Step 4: merge the 100 nodes into metaclusters ----------------------
         # The 100 SOM nodes are too many to interpret. ConsensusClusterPlus groups
         # them into a stable, smaller set (we later use meta_k metaclusters).
         som_codes <- fsom$map$codes
         set.seed(2026)
         cc <- ConsensusClusterPlus(
                  t(som_codes),
                  maxK       = consensus_maxk,
                  reps       = consensus_reps,
                  pItem      = 0.9,
                  pFeature   = 1,
                  clusterAlg = "hc",
                  distance   = "euclidean",
                  seed       = 2026,
                  plot       = "png",
                  title      = file.path(dir_results, "Consensus")
         )
         saveRDS(cc, file.path(dir_results, "ConsensusCluster.rds"))
         
         ## ---- Step 5: give every CELL a metacluster label ------------------------
         # cc[[meta_k]]$consensusClass says which metacluster each NODE belongs to.
         # fsom$map$mapping[, 1] says which NODE each CELL belongs to.
         # Combining them gives each cell its metacluster.
         node_to_meta <- cc[[meta_k]]$consensusClass
         cell_to_node <- fsom$map$mapping[, 1]
         cell_meta    <- node_to_meta[cell_to_node]
         
         saveRDS(cell_meta, file.path(dir_results, "Cell_Metaclusters.rds"))
         cat("Metacluster sizes (k =", meta_k, "):\n")
         print(table(cell_meta))
         
         ## ---- Step 6: UMAP - a 2D map for visualising the cells ------------------
         set.seed(2026)
         umap <- uwot::umap(
                  expr_all,
                  n_neighbors = umap_neighbors,
                  min_dist    = umap_min_dist,
                  metric      = "euclidean",
                  verbose     = TRUE
         )
         saveRDS(umap, file.path(dir_results, "UMAP.rds"))
         
         ## ---- Step 7: save the transformed matrix + sample IDs for later ---------
         saveRDS(expr_all,  file.path(dir_results, "expr_all_asinh.rds"))
         saveRDS(sample_id, file.path(dir_results, "sample_id.rds"))
         
         rm(fsom, cc, som_codes, node_to_meta, cell_to_node)
         gc()
         mem_report("END Stage 03")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 03b : QUICK BATCH-CHECK UMAPs  (diagnostic only)                   ##
##                                                                           ##
##   Two simple UMAPs to eyeball whether a batch effect is a problem:         ##
##     - coloured by SAMPLE : if each sample forms its own separate island    ##
##                            (little overlap), suspect a technical effect.   ##
##     - coloured by DONOR  : tells a donor effect apart from a sample effect.##
##   This stage only draws pictures; it changes no analysis object.           ##
###############################################################################
###############################################################################

if ("3b" %in% run_stages) {
         
         mem_report("START Stage 03b (batch-check UMAPs)")
         
         if (!exists("umap"))      umap      <- readRDS(file.path(dir_results, "UMAP.rds"))
         if (!exists("sample_id")) sample_id <- readRDS(file.path(dir_results, "sample_id.rds"))
         if (!exists("metadata"))  metadata  <- readRDS(file.path(dir_results, "metadata_downsampled.rds"))
         
         ## ---- Step 1: put the UMAP coordinates + labels into one table -----------
         diag_df <- tibble(
                  UMAP1  = umap[, 1],
                  UMAP2  = umap[, 2],
                  Sample = sample_id
         )
         # Add the donor (patient) column by matching on sample_id.
         meta_small <- metadata[, c("sample_id", "patient")]
         diag_df <- left_join(diag_df, meta_small,
                              by = c("Sample" = "sample_id"))
         
         ## ---- Step 2: subsample so plotting stays light --------------------------
         if (nrow(diag_df) > max_umap_points) {
                  rows_to_plot <- sample.int(nrow(diag_df), max_umap_points)
         } else {
                  rows_to_plot <- seq_len(nrow(diag_df))
         }
         diag_plot <- diag_df[rows_to_plot, ]
         
         ## ---- Step 3: UMAP coloured by SAMPLE ------------------------------------
         p_sample <- ggplot(diag_plot, aes(x = UMAP1, y = UMAP2, colour = Sample)) +
                  rasterise(geom_point(size = 0.12, alpha = 0.5), dpi = raster_dpi) +
                  coord_equal() +
                  guides(colour = guide_legend(
                           override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
                  labs(title = "UMAP coloured by sample")
         save_fig(p_sample, file.path(dir_umap, "DIAG_UMAP_by_sample"),
                  width = 8, height = 6)
         
         ## ---- Step 4: UMAP coloured by DONOR -------------------------------------
         p_donor <- ggplot(diag_plot, aes(x = UMAP1, y = UMAP2, colour = patient)) +
                  rasterise(geom_point(size = 0.12, alpha = 0.5), dpi = raster_dpi) +
                  coord_equal() +
                  guides(colour = guide_legend(
                           override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
                  labs(title = "UMAP coloured by donor")
         save_fig(p_donor, file.path(dir_umap, "DIAG_UMAP_by_donor"),
                  width = 7.5, height = 6)
         
         rm(diag_df, diag_plot, meta_small, rows_to_plot, p_sample, p_donor)
         gc()
         mem_report("END Stage 03b - look at Figures/UMAP/DIAG_UMAP_by_sample and _by_donor")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 04 : MASTER PER-CELL TABLE + SUMMARY TABLES                        ##
##                                                                           ##
###############################################################################
###############################################################################

if ("4" %in% run_stages) {
         
         mem_report("START Stage 04")
         
         if (!exists("expr_all"))  expr_all  <- readRDS(file.path(dir_results, "expr_all_asinh.rds"))
         if (!exists("sample_id")) sample_id <- readRDS(file.path(dir_results, "sample_id.rds"))
         if (!exists("cell_meta")) cell_meta <- readRDS(file.path(dir_results, "Cell_Metaclusters.rds"))
         if (!exists("metadata"))  metadata  <- readRDS(file.path(dir_results, "metadata_downsampled.rds"))
         
         ## ---- Step 1: build the master table, one row per cell -------------------
         # Start from the marker matrix, rename detector columns to marker names,
         # then add the cluster, the sample, and the sample's metadata.
         cell_data <- as.data.frame(expr_all)
         
         # Rename columns: for each column, look up its marker name.
         colnames(cell_data) <- marker_lookup[colnames(cell_data)]
         
         cell_data$Cluster <- factor(cell_meta)
         cell_data$Sample  <- sample_id
         cell_data <- left_join(cell_data, metadata,
                                by = c("Sample" = "sample_id"))
         
         saveRDS(cell_data, file.path(dir_results, "cell_data.rds"))
         rm(expr_all)
         gc()
         
         ## ---- Step 2: per-cluster phenotype summary ------------------------------
         # For each cluster we report: how many cells, and for every marker its
         # median, mean, and the % of cells that are "positive" (value > 1).
         cluster_summary <- cell_data %>%
                  group_by(Cluster) %>%
                  summarise(
                           Cells = n(),
                           across(all_of(markers),
                                  list(
                                           Median   = ~ median(.x),
                                           Mean     = ~ mean(.x),
                                           Positive = ~ 100 * mean(.x > 1)
                                  ),
                                  .names = "{.col}_{.fn}"),
                           .groups = "drop"
                  )
         write.csv(cluster_summary,
                   file.path(dir_results, "Cluster_Summary.csv"),
                   row.names = FALSE)
         
         ## ---- Step 3: cluster sizes (counts and % of all cells) ------------------
         cluster_size <- cell_data %>%
                  count(Cluster, name = "Cells") %>%
                  mutate(Percent = round(100 * Cells / sum(Cells), 2)) %>%
                  arrange(desc(Cells))
         write.csv(cluster_size,
                   file.path(dir_results, "Cluster_Size.csv"),
                   row.names = FALSE)
         
         ## ---- Step 4: per-sample cluster frequencies -----------------------------
         # For each sample, what fraction of its cells fall in each cluster?
         # These frequencies feed the PCA and any abundance statistics.
         cluster_abundance <- cell_data %>%
                  count(Sample, Cluster, tissue, name = "n") %>%
                  group_by(Sample) %>%
                  mutate(Frequency = n / sum(n)) %>%
                  ungroup()
         
         saveRDS(cluster_abundance,
                 file.path(dir_results, "cluster_frequency.rds"))
         write.csv(cluster_abundance,
                   file.path(dir_results, "Cluster_Abundance.csv"),
                   row.names = FALSE)
         
         rm(cluster_summary, cluster_size)
         gc()
         mem_report("END Stage 04")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 05 : PUBLICATION FIGURES                                           ##
##     FIG 1  Heatmaps    (marker x sample ; cluster x marker)                ##
##     FIG 2  UMAP panels (all cells ; split by tissue ; per-marker)          ##
##     FIG 3  Ridgelines  (marker distributions by cluster ; by tissue)       ##
##     FIG 4  PCA of samples (based on cluster frequencies)                   ##
##                                                                           ##
###############################################################################
###############################################################################

if ("5" %in% run_stages) {
         
         mem_report("START Stage 05")
         
         ## ---- Load what the figures need -----------------------------------------
         if (!exists("cell_data"))         cell_data         <- readRDS(file.path(dir_results, "cell_data.rds"))
         if (!exists("umap"))              umap              <- readRDS(file.path(dir_results, "UMAP.rds"))
         if (!exists("cluster_abundance")) cluster_abundance <- readRDS(file.path(dir_results, "cluster_frequency.rds"))
         
         # Safety net: if cell_data still has detector names, rename them to markers.
         if (!all(markers %in% colnames(cell_data))) {
                  pos <- match(names(marker_lookup), colnames(cell_data))
                  ok  <- !is.na(pos)
                  colnames(cell_data)[pos[ok]] <- marker_lookup[ok]
         }
         stopifnot(all(markers %in% colnames(cell_data)))
         cell_data$Cluster <- factor(cell_data$Cluster)
         
         # Give each cluster a fixed colour.
         cluster_levels <- levels(cell_data$Cluster)
         cluster_pal    <- cluster_cols[seq_along(cluster_levels)]
         names(cluster_pal) <- cluster_levels
         
         mem_report("Figure inputs ready")
         
         #######################################################################
         # FIG 1 : HEATMAPS
         #######################################################################
         
         ## ---- FIG 1A: marker (rows) x sample (columns) ---------------------------
         # For each sample, the median of each marker; then z-score each marker row.
         sample_medians <- cell_data %>%
                  group_by(Sample) %>%
                  summarise(across(all_of(markers), ~ median(.x)), .groups = "drop")
         
         # Turn into a matrix with samples as columns and markers as rows.
         mat_ms <- as.matrix(sample_medians[, markers])
         rownames(mat_ms) <- sample_medians$Sample
         mat_ms <- t(mat_ms)                 # now markers = rows, samples = columns
         mat_ms_z <- zscore_rows(mat_ms)     # z-score each marker
         
         # Aesthetic follows the reference figure: yellow-black-blue map, symmetric
         # Gene Z-score legend, thin dendrograms, no column labels, no annotation bar.
         ht_ms <- Heatmap(
                  mat_ms_z,
                  name = "Gene\nZ-score",
                  col  = heat_ramp,
                  cluster_rows    = TRUE,
                  cluster_columns = TRUE,
                  clustering_method_rows      = "ward.D2",
                  clustering_method_columns   = "ward.D2",
                  clustering_distance_rows    = "euclidean",
                  clustering_distance_columns = "euclidean",
                  show_column_names  = FALSE,                 # too many samples to label
                  row_names_gp       = grid::gpar(fontsize = 8),
                  row_dend_width     = grid::unit(12, "mm"),
                  column_dend_height = grid::unit(15, "mm"),
                  heatmap_legend_param = list(
                           title_gp    = grid::gpar(fontsize = 8, fontface = "bold"),
                           labels_gp   = grid::gpar(fontsize = 7),
                           at          = c(-4, -2, 0, 2, 4),
                           legend_height = grid::unit(30, "mm")
                  ),
                  border = FALSE
         )
         
         save_heatmap(ht_ms, file.path(dir_heat, "FIG1A_Marker_Sample_Heatmap"),
                      width = 8, height = 5)
         
         ## ---- FIG 1B: cluster (rows) x marker (columns) --------------------------
         # This is the "phenotype key" you read to name each cluster.
         cluster_medians <- cell_data %>%
                  group_by(Cluster) %>%
                  summarise(across(all_of(markers), ~ median(.x)), .groups = "drop")
         
         mat_cm <- as.matrix(cluster_medians[, markers])
         rownames(mat_cm) <- paste0("MC", cluster_medians$Cluster)
         
         # Z-score each marker (each column) so colours compare markers fairly.
         mat_cm_z <- scale(mat_cm)
         mat_cm_z[is.na(mat_cm_z)] <- 0
         
         ht_cm <- Heatmap(
                  mat_cm_z,
                  name = "Gene\nZ-score",
                  col  = heat_ramp,
                  cluster_rows    = TRUE,
                  cluster_columns = TRUE,
                  clustering_method_rows    = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_gp       = grid::gpar(fontsize = 8),
                  column_names_gp    = grid::gpar(fontsize = 8),
                  row_dend_width     = grid::unit(12, "mm"),
                  column_dend_height = grid::unit(12, "mm"),
                  heatmap_legend_param = list(
                           title_gp  = grid::gpar(fontsize = 8, fontface = "bold"),
                           labels_gp = grid::gpar(fontsize = 7),
                           at        = c(-4, -2, 0, 2, 4)
                  ),
                  border = FALSE
         )
         
         save_heatmap(ht_cm, file.path(dir_heat, "FIG1B_Cluster_Marker_Heatmap"),
                      width = 6, height = 7)
         
         rm(sample_medians, cluster_medians, mat_ms, mat_ms_z,
            mat_cm, mat_cm_z, ht_ms, ht_cm)
         gc()
         mem_report("FIG 1 done")
         
         #######################################################################
         # FIG 2 : UMAP PANELS
         #######################################################################
         
         # Put UMAP coordinates, cluster and tissue into one plotting table.
         umap_df <- tibble(
                  UMAP1   = umap[, 1],
                  UMAP2   = umap[, 2],
                  Cluster = cell_data$Cluster,
                  tissue  = cell_data$tissue
         )
         
         # Subsample cells so the plots render quickly and stay small.
         if (nrow(umap_df) > max_umap_points) {
                  keep_rows <- sample.int(nrow(umap_df), max_umap_points)
         } else {
                  keep_rows <- seq_len(nrow(umap_df))
         }
         umap_plot_df <- umap_df[keep_rows, ]
         
         ## ---- FIG 2A: one UMAP, coloured by cluster ------------------------------
         p_master <- ggplot(umap_plot_df,
                            aes(x = UMAP1, y = UMAP2, colour = Cluster)) +
                  rasterise(geom_point(size = 0.15, alpha = 0.6), dpi = raster_dpi) +
                  scale_colour_manual(values = cluster_pal, name = "Metacluster") +
                  coord_equal() +
                  guides(colour = guide_legend(
                           override.aes = list(size = 3, alpha = 1), ncol = 1)) +
                  labs(title = "Mucosal immune landscape (FlowSOM metaclusters)")
         save_fig(p_master, file.path(dir_umap, "FIG2A_UMAP_master"),
                  width = 7, height = 5.5)
         
         ## ---- FIG 2B: same UMAP, one panel per tissue ----------------------------
         # Grey backdrop = all cells; coloured = cells of that tissue panel.
         backdrop <- umap_plot_df[, c("UMAP1", "UMAP2")]
         p_facets <- ggplot() +
                  rasterise(geom_point(data = backdrop,
                                       aes(x = UMAP1, y = UMAP2),
                                       colour = "grey88", size = 0.10), dpi = raster_dpi) +
                  rasterise(geom_point(data = umap_plot_df,
                                       aes(x = UMAP1, y = UMAP2, colour = Cluster),
                                       size = 0.14, alpha = 0.7), dpi = raster_dpi) +
                  scale_colour_manual(values = cluster_pal, guide = "none") +
                  facet_wrap(~ tissue, ncol = 2) +
                  coord_equal() +
                  theme(strip.text = element_text(face = "bold", size = 10))
         save_fig(p_facets, file.path(dir_umap, "FIG2B_UMAP_by_tissue"),
                  width = 8, height = 7)
         
         ## ---- FIG 2C: one small UMAP per marker, coloured by expression ----------
         feat_df <- bind_cols(umap_plot_df[, c("UMAP1", "UMAP2")],
                              cell_data[keep_rows, markers])
         
         # A helper that draws one marker's feature plot.
         feature_plot <- function(marker_name) {
                  ggplot(feat_df, aes(x = UMAP1, y = UMAP2,
                                      colour = .data[[marker_name]])) +
                           rasterise(geom_point(size = 0.12, alpha = 0.7), dpi = raster_dpi) +
                           scale_colour_viridis_c(option = feature_option, name = NULL) +
                           coord_equal() +
                           theme(legend.key.width  = grid::unit(2, "mm"),
                                 legend.key.height = grid::unit(6, "mm"),
                                 axis.title = element_blank(),
                                 axis.text  = element_blank(),
                                 axis.ticks = element_blank()) +
                           ggtitle(marker_name)
         }
         
         # Make one plot per marker, then arrange them in a 4-column grid.
         feature_plots <- lapply(markers, feature_plot)
         p_features <- wrap_plots(feature_plots, ncol = 4)
         save_fig(p_features, file.path(dir_umap, "FIG2C_Feature_grid"),
                  width = 10, height = 8)
         
         rm(umap_df, umap_plot_df, backdrop, feat_df, feature_plots,
            p_master, p_facets, p_features, keep_rows)
         gc()
         mem_report("FIG 2 done")
         
         #######################################################################
         # FIG 3 : RIDGELINE DENSITY PLOTS
         #######################################################################
         
         ## ---- FIG 3A: marker distributions split by cluster ----------------------
         # We cap the cells per cluster so the density estimation is fast.
         # slice_sample() errors if a group has fewer rows than n, so we take the
         # smaller of the cap and the group size (small clusters keep all their cells).
         ridge_cluster <- cell_data %>%
                  dplyr::select(Cluster, all_of(markers)) %>%
                  group_by(Cluster) %>%
                  sample_up_to(cap_per_cluster)
         
         # Convert to long format: one row per (cell, marker) so we can facet.
         ridge_cluster_long <- ridge_cluster %>%
                  pivot_longer(cols = all_of(markers),
                               names_to = "Marker", values_to = "Expression") %>%
                  mutate(Marker = factor(Marker, levels = markers))
         
         p_ridge_cluster <- ggplot(ridge_cluster_long,
                                   aes(x = Expression, y = Cluster, fill = Cluster)) +
                  geom_density_ridges(scale = 2.2, linewidth = 0.25,
                                      rel_min_height = 0.01, alpha = 0.9) +
                  scale_fill_manual(values = cluster_pal, guide = "none") +
                  facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
                  labs(x = "Arcsinh expression", y = "Metacluster") +
                  theme_ridges(font_size = 9, grid = TRUE) +
                  theme(strip.text = element_text(face = "bold"),
                        axis.title = element_text(face = "bold"))
         save_fig(p_ridge_cluster,
                  file.path(dir_ridge, "FIG3A_Ridgeline_by_cluster"),
                  width = 11, height = 9)
         
         rm(ridge_cluster, ridge_cluster_long, p_ridge_cluster)
         gc()
         
         ## ---- FIG 3B: marker distributions split by tissue -----------------------
         # Same safe-cap trick as FIG 3A (a small tissue keeps all its cells).
         ridge_tissue <- cell_data %>%
                  dplyr::select(tissue, all_of(markers)) %>%
                  group_by(tissue) %>%
                  sample_up_to(cap_per_tissue)
         
         ridge_tissue_long <- ridge_tissue %>%
                  pivot_longer(cols = all_of(markers),
                               names_to = "Marker", values_to = "Expression") %>%
                  mutate(Marker = factor(Marker, levels = markers))
         
         p_ridge_tissue <- ggplot(ridge_tissue_long,
                                  aes(x = Expression, y = tissue, fill = tissue)) +
                  geom_density_ridges(scale = 1.8, linewidth = 0.3,
                                      rel_min_height = 0.01, alpha = 0.9) +
                  scale_fill_manual(values = tissue_cols, guide = "none") +
                  facet_wrap(~ Marker, ncol = 4, scales = "free_x") +
                  labs(x = "Arcsinh expression", y = "Tissue") +
                  theme_ridges(font_size = 9, grid = TRUE) +
                  theme(strip.text = element_text(face = "bold"),
                        axis.title = element_text(face = "bold"))
         save_fig(p_ridge_tissue,
                  file.path(dir_ridge, "FIG3B_Ridgeline_by_tissue"),
                  width = 11, height = 7)
         
         rm(ridge_tissue, ridge_tissue_long, p_ridge_tissue)
         gc()
         mem_report("FIG 3 done")
         
         #######################################################################
         # FIG 4 : PCA OF SAMPLES (based on cluster frequencies)
         #######################################################################
         
         ## ---- Step 1: build a samples x clusters frequency table -----------------
         # pivot_wider turns the long abundance table into a wide matrix:
         # one row per sample, one column per cluster, cells = frequency.
         freq_wide <- cluster_abundance %>%
                  dplyr::select(Sample, Cluster, Frequency) %>%
                  mutate(Cluster = paste0("MC", Cluster)) %>%
                  pivot_wider(names_from = Cluster, values_from = Frequency,
                              values_fill = 0)
         
         # Sample metadata (tissue + donor) to colour and label the plot.
         sample_meta <- cell_data %>%
                  distinct(Sample, tissue, patient) %>%
                  filter(Sample %in% freq_wide$Sample)
         
         # Turn the frequency table into a plain numeric matrix.
         freq_mat <- as.matrix(freq_wide[, -1])          # drop the Sample column
         rownames(freq_mat) <- freq_wide$Sample
         
         ## ---- Step 2: CLR transform (proper for compositional data) --------------
         # Frequencies always add up to 1, which distorts ordinary PCA. The
         # centred-log-ratio (CLR) removes that constraint. We add a tiny
         # pseudocount first so log() never sees a zero.
         freq_log <- log(freq_mat + pseudocount)
         row_means <- rowMeans(freq_log)
         freq_clr  <- freq_log - row_means               # subtract each row's mean
         
         ## ---- Step 3: run PCA ----------------------------------------------------
         # Drop any cluster whose values never vary (would break scaling).
         col_var  <- matrixStats::colVars(freq_clr)
         freq_clr <- freq_clr[, col_var > 0, drop = FALSE]
         
         pca <- prcomp(freq_clr, center = TRUE, scale. = TRUE)
         
         # Percent of variance each PC explains (for the axis labels).
         var_explained <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
         
         # Scores = where each sample sits in PC space.
         pca_scores <- as.data.frame(pca$x[, 1:2])
         pca_scores$Sample <- rownames(pca_scores)
         pca_scores <- left_join(pca_scores, sample_meta, by = "Sample")
         
         ## ---- FIG 4A: score plot -------------------------------------------------
         p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = tissue)) +
                  stat_ellipse(aes(group = tissue), type = "norm",
                               linewidth = 0.4, alpha = 0.5) +
                  geom_point(size = 3, alpha = 0.9) +
                  geom_text_repel(aes(label = patient), size = 2.5,
                                  max.overlaps = 20, show.legend = FALSE) +
                  scale_colour_manual(values = tissue_cols, name = "Tissue") +
                  coord_equal() +
                  labs(title = "PCA of mucosal samples (cluster-frequency space)",
                       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                       y = sprintf("PC2 (%.1f%%)", var_explained[2]))
         save_fig(p_pca, file.path(dir_pca, "FIG4A_PCA_scores"),
                  width = 7, height = 6)
         
         ## ---- FIG 4B: scree plot (how much variance per PC) ----------------------
         n_show <- min(10, length(var_explained))
         scree_df <- tibble(
                  PC       = factor(paste0("PC", seq_len(n_show)),
                                    levels = paste0("PC", seq_len(n_show))),
                  Variance = var_explained[seq_len(n_show)]
         )
         p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
                  geom_col(fill = "#4477AA") +
                  geom_text(aes(label = sprintf("%.1f%%", Variance)),
                            vjust = -0.4, size = 2.6) +
                  labs(x = NULL, y = "Variance explained (%)", title = "Scree plot")
         save_fig(p_scree, file.path(dir_pca, "FIG4B_PCA_scree"),
                  width = 5, height = 4)
         
         ## ---- FIG 4C: loadings (which clusters drive the separation) -------------
         loadings_df <- as.data.frame(pca$rotation[, 1:2])
         loadings_df$Cluster <- rownames(loadings_df)
         
         p_load <- ggplot(loadings_df, aes(x = PC1, y = PC2)) +
                  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                               arrow = arrow(length = grid::unit(2, "mm")),
                               colour = "grey50", linewidth = 0.3) +
                  geom_text_repel(aes(label = Cluster), size = 2.6, max.overlaps = 30) +
                  labs(title = "PCA loadings (cluster contributions)",
                       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                       y = sprintf("PC2 (%.1f%%)", var_explained[2]))
         save_fig(p_load, file.path(dir_pca, "FIG4C_PCA_loadings"),
                  width = 6, height = 5)
         
         ## ---- Step 4: save PCA outputs -------------------------------------------
         saveRDS(pca, file.path(dir_results, "PCA_object.rds"))
         write.csv(pca_scores,  file.path(dir_results, "PCA_scores.csv"),   row.names = FALSE)
         write.csv(loadings_df, file.path(dir_results, "PCA_loadings.csv"), row.names = FALSE)
         
         rm(freq_wide, freq_mat, freq_log, freq_clr, sample_meta, col_var,
            pca, var_explained, pca_scores, scree_df, loadings_df,
            p_pca, p_scree, p_load, row_means)
         gc()
         mem_report("FIG 4 done")
}


###############################################################################
###############################################################################
##                                                                           ##
##   STAGE 05b : FEATURE-DRIVEN PCA  (samples described by marker medians)    ##
##                                                                           ##
##   Instead of cluster frequencies, here each sample is summarised by the    ##
##   MEDIAN of each marker. The PCA loadings then point directly at markers   ##
##   (e.g. CD103, CD69), which is easy to interpret biologically.             ##
##                                                                           ##
##   We run it twice:                                                         ##
##     "Full12"    : all 12 markers  -> overall phenotypic separation         ##
##     "Residency" : CD69/CD103/CD45RO/CCR7 -> targeted tissue-residency test ##
##                                                                           ##
##   Why scale. = TRUE: markers with a big numeric range (like CD3) would     ##
##   otherwise dominate. Scaling puts every marker on an equal footing        ##
##   (this matches CATALYST pseudobulk-MDS; Nowicka et al., F1000Res 2019).   ##
###############################################################################
###############################################################################

if ("5b" %in% run_stages) {
         
         mem_report("START Stage 05b (feature-driven PCA)")
         
         if (!exists("cell_data")) cell_data <- readRDS(file.path(dir_results, "cell_data.rds"))
         cell_data$Cluster <- factor(cell_data$Cluster)
         
         # The residency/memory markers used for the targeted PCA.
         residency_markers <- c("CD69", "CD103", "CD45RO", "CCR7")
         stopifnot(all(residency_markers %in% markers))
         
         ## -------------------------------------------------------------------------
         ## One function that does the whole feature-PCA for a chosen marker set.
         ## marker_set : which markers to use ; tag : short name for the files.
         ## -------------------------------------------------------------------------
         run_feature_pca <- function(marker_set, tag) {
                  
                  ## Step 1: sample x marker median matrix ---------------------------------
                  med_tbl <- cell_data %>%
                           group_by(Sample) %>%
                           summarise(across(all_of(marker_set), ~ median(.x)), .groups = "drop")
                  
                  mat <- as.matrix(med_tbl[, marker_set])
                  rownames(mat) <- med_tbl$Sample
                  
                  # Drop any marker that does not vary (would break scaling).
                  keep_marker <- matrixStats::colVars(mat) > 0
                  if (any(!keep_marker)) {
                           message("Dropping constant marker(s): ",
                                   paste(colnames(mat)[!keep_marker], collapse = ", "))
                  }
                  mat <- mat[, keep_marker, drop = FALSE]
                  
                  ## Step 2: PCA (centre + scale each marker) ------------------------------
                  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
                  var_explained <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
                  
                  ## Step 3: sample metadata for colouring/labels --------------------------
                  smeta <- cell_data %>%
                           distinct(Sample, tissue, patient) %>%
                           filter(Sample %in% rownames(pca$x))
                  
                  scores <- as.data.frame(pca$x[, 1:2])
                  scores$Sample <- rownames(scores)
                  scores <- left_join(scores, smeta, by = "Sample")
                  
                  ## Step 4: loadings, scaled so the arrows fit over the points ------------
                  loadings <- as.data.frame(pca$rotation[, 1:2])
                  loadings$Marker <- rownames(loadings)
                  
                  # Scale factor: make the longest arrow about as long as the point cloud.
                  score_span   <- max(abs(scores$PC1), abs(scores$PC2))
                  loading_span <- max(abs(as.matrix(loadings[, c("PC1", "PC2")])))
                  arrow_scale  <- 0.9 * score_span / loading_span
                  
                  loadings$PC1s <- loadings$PC1 * arrow_scale
                  loadings$PC2s <- loadings$PC2 * arrow_scale
                  
                  ## Step 5A: score plot ---------------------------------------------------
                  p_scores <- ggplot(scores, aes(x = PC1, y = PC2, colour = tissue)) +
                           stat_ellipse(aes(group = tissue), type = "norm",
                                        linewidth = 0.4, alpha = 0.5) +
                           geom_point(size = 3, alpha = 0.9) +
                           geom_text_repel(aes(label = patient), size = 2.5,
                                           max.overlaps = 20, show.legend = FALSE) +
                           scale_colour_manual(values = tissue_cols, name = "Tissue") +
                           coord_equal() +
                           labs(title = sprintf("PCA scores - %s markers", tag),
                                x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                                y = sprintf("PC2 (%.1f%%)", var_explained[2]))
                  
                  ## Step 5B: biplot (points + marker arrows) - the interpretable figure ---
                  p_biplot <- ggplot() +
                           stat_ellipse(data = scores,
                                        aes(x = PC1, y = PC2, group = tissue, colour = tissue),
                                        type = "norm", linewidth = 0.4, alpha = 0.4,
                                        show.legend = FALSE) +
                           geom_point(data = scores,
                                      aes(x = PC1, y = PC2, colour = tissue),
                                      size = 3, alpha = 0.9) +
                           geom_segment(data = loadings,
                                        aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
                                        arrow = arrow(length = grid::unit(2.2, "mm")),
                                        colour = "grey25", linewidth = 0.4) +
                           geom_text_repel(data = loadings,
                                           aes(x = PC1s, y = PC2s, label = Marker),
                                           size = 3, fontface = "bold",
                                           colour = "grey15", max.overlaps = 30) +
                           scale_colour_manual(values = tissue_cols, name = "Tissue") +
                           coord_equal() +
                           labs(title = sprintf("PCA biplot - %s markers", tag),
                                subtitle = "Arrows = markers driving the separation",
                                x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                                y = sprintf("PC2 (%.1f%%)", var_explained[2]))
                  
                  ## Step 5C: bar chart of marker contributions to PC1 and PC2 -------------
                  contrib <- as.data.frame(pca$rotation[, 1:2])
                  contrib$Marker <- rownames(contrib)
                  contrib_long <- pivot_longer(contrib, cols = c(PC1, PC2),
                                               names_to = "PC", values_to = "Loading")
                  contrib_long$Marker <- fct_reorder(contrib_long$Marker,
                                                     abs(contrib_long$Loading))
                  
                  p_contrib <- ggplot(contrib_long,
                                      aes(x = Loading, y = Marker, fill = Loading > 0)) +
                           geom_col() +
                           geom_vline(xintercept = 0, linewidth = 0.3) +
                           facet_wrap(~ PC, nrow = 1) +
                           scale_fill_manual(values = c("TRUE" = "#D7191C", "FALSE" = "#2C7BB6"),
                                             guide = "none") +
                           labs(x = "Loading", y = NULL,
                                title = sprintf("Marker contributions - %s", tag))
                  
                  ## Step 6: save the figures ----------------------------------------------
                  stem <- file.path(dir_pca, paste0("FIG5_FeaturePCA_", tag))
                  save_fig(p_scores,  paste0(stem, "_scores"),  width = 7,   height = 6)
                  save_fig(p_biplot,  paste0(stem, "_biplot"),  width = 7.5, height = 6.5)
                  save_fig(p_contrib, paste0(stem, "_contrib"), width = 8,   height = 4)
                  
                  # A combined panel (biplot next to the contribution bars).
                  p_combo <- (p_biplot | p_contrib) + plot_annotation(tag_levels = "A")
                  save_fig(p_combo, paste0(stem, "_composite"), width = 13, height = 6)
                  
                  ## Step 7: save the numbers behind the plots -----------------------------
                  write.csv(scores,
                            file.path(dir_results, paste0("FeaturePCA_", tag, "_scores.csv")),
                            row.names = FALSE)
                  write.csv(loadings[, c("Marker", "PC1", "PC2")],
                            file.path(dir_results, paste0("FeaturePCA_", tag, "_loadings.csv")),
                            row.names = FALSE)
                  
                  # Return the key pieces in case you want them at the console.
                  invisible(list(pca = pca, scores = scores,
                                 loadings = loadings, var_explained = var_explained))
         }
         
         ## ---- Run the function twice ---------------------------------------------
         result_full      <- run_feature_pca(markers,           tag = "Full12")
         gc()
         result_residency <- run_feature_pca(residency_markers, tag = "Residency")
         gc()
         
         ## ---- Print the strongest markers on PC1 for each run --------------------
         cat("\nTop PC1 markers (full 12-marker panel):\n")
         full_load <- result_full$loadings
         full_load <- full_load[order(-abs(full_load$PC1)), c("Marker", "PC1")]
         print(head(full_load, 5))
         
         cat("\nTop PC1 markers (residency panel):\n")
         res_load <- result_residency$loadings
         res_load <- res_load[order(-abs(res_load$PC1)), c("Marker", "PC1")]
         print(head(res_load, 4))
         
         rm(result_full, result_residency, full_load, res_load)
         gc()
         mem_report("END Stage 05b (feature-driven PCA)")
}


###############################################################################
# RECORD THE SESSION (which package versions were used) - good for methods
###############################################################################

writeLines(capture.output(sessionInfo()),
           file.path(dir_results, "sessionInfo_MasterPipeline.txt"))

cat("\n---------------------------------------------------------------\n")
cat("PIPELINE FINISHED.  Stages run:", paste(run_stages, collapse = ", "), "\n")
cat("Checkpoints & tables : ./Results/\n")
cat("Figures              : ./Figures/{QC,Heatmaps,UMAP,Ridgeline,PCA}/\n")
cat("---------------------------------------------------------------\n")
mem_report("END OF PIPELINE")
