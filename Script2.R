###############################################################################
# MUCOSAL HIGH-DIMENSIONAL ANALYSIS
# SCRIPT 03
# Create SingleCellExperiment
###############################################################################

rm(list = ls())
graphics.off()
gc()

library(flowCore)
library(tidyverse)

fs <- readRDS("Results/fs_10k.rds")
metadata <- readRDS("Results/metadata.rds")
panel <- readRDS("Results/panel.rds")

metadata <- metadata %>%
         filter(sample_id %in% sampleNames(fs)) %>%
         arrange(match(sample_id, sampleNames(fs)))

stopifnot(all(metadata$sample_id == sampleNames(fs)))

saveRDS(metadata, "Results/metadata_10k.rds")

dim(exprs(fs[[1]]))
pData(parameters(fs[[1]]))[, c("name", "desc")]


library(FlowSOM)
packageVersion("FlowSOM")

###############################################################################
# Packages
###############################################################################

library(flowCore)
library(FlowSOM)
library(uwot)
library(tidyverse)
library(RColorBrewer)
library(ConsensusClusterPlus)

set.seed(2026)

cluster_channels <- c(
         
         "FJComp-BV785-A-1",          # CD3
         "FJComp-APC-Cy7-A-1",        # CD4
         "FJComp-AF700-A-1",          # CD8
         "FJComp-BV750-A-1",          # CD19
         "FJComp-BV421-A-1",          # CD14
         "FJComp-BV650-A-1",          # CD16
         "FJComp-PE-A-1",             # CD66b
         "FJComp-APC-A-1",            # CD56
         "FJComp-BV605-A-1",          # CD69
         "FJComp-PE-Dazzle594-A-1",   # CD103
         "FJComp-FITC-A-1",           # CD45RO
         "FJComp-PE-Cy7-A-1"          # CCR7
         
)

###############################################################################
# Combine all samples
###############################################################################

expr_list <- vector("list", length(fs))

sample_id <- vector("list", length(fs))

for(i in seq_along(fs)){
         
         expr_list[[i]] <- exprs(fs[[i]])[, cluster_channels]
         
         sample_id[[i]] <- rep(sampleNames(fs)[i],
                               nrow(expr_list[[i]]))
         
}

expr_all <- do.call(rbind, expr_list)

sample_id <- unlist(sample_id)

rm(expr_list)

gc()
dim(expr_all)

summary(expr_all[,1])

range(expr_all)

###############################################################################
# Arcsinh transformation
###############################################################################

cofactor <- 150

expr_all <- asinh(expr_all / cofactor)

gc()

range(expr_all)

summary(expr_all[,1])

###############################################################################
# FlowSOM clustering
###############################################################################

set.seed(2026)

fsom <- ReadInput(
         expr_all,
         transform = FALSE,
         scale = FALSE
)

gc()

fsom <- BuildSOM(
         fsom,
         colsToUse = 1:ncol(expr_all),
         xdim = 10,
         ydim = 10,
         rlen = 20
)

gc()

fsom <- BuildMST(fsom)

gc()
PlotStars(fsom)
head(fsom$map$mapping)
table(fsom$map$mapping[,1])
saveRDS(
         fsom,
         "Results/fsom_100nodes.rds"
)

gc()


###############################################################################
# Extract SOM codes
###############################################################################

codes <- fsom$map$codes

dim(codes)

head(codes)
###############################################################################
# Consensus metaclustering
###############################################################################

library(ConsensusClusterPlus)

set.seed(2026)

cc <- ConsensusClusterPlus(
         
         t(codes),
         
         maxK = 30,
         
         reps = 100,
         
         pItem = 0.9,
         
         pFeature = 1,
         
         clusterAlg = "hc",
         
         distance = "euclidean",
         
         seed = 2026,
         
         plot = "png",
         
         title = "Results/Consensus"
         
)


meta10 <- cc[[10]]$consensusClass
meta15 <- cc[[15]]$consensusClass
meta20 <- cc[[20]]$consensusClass
meta25 <- cc[[25]]$consensusClass


cell_cluster <- fsom$map$mapping[,1]

cell_meta20 <- meta20[cell_cluster]

table(cell_meta20)

saveRDS(cc, "Results/ConsensusCluster.rds")

saveRDS(cell_meta20, "Results/Cell_Metaclusters.rds")

gc()


library(uwot)

set.seed(2026)

umap <- uwot::umap(
         
         expr_all,
         
         n_neighbors = 30,
         
         min_dist = 0.3,
         
         metric = "euclidean",
         
         verbose = TRUE
         
)

saveRDS(umap, "Results/UMAP.rds")

plot_df <- data.frame(
         
         UMAP1 = umap[,1],
         
         UMAP2 = umap[,2],
         
         Cluster = factor(cell_meta20),
         
         Sample = sample_id
         
)


library(ggplot2)

ggplot(plot_df,
       aes(UMAP1, UMAP2, colour = Cluster)) +
         
         geom_point(
                  size = 0.15,
                  alpha = 0.5
         ) +
         
         theme_classic(base_size = 14) +
         
         coord_equal()

###############################################################################
# Median expression per metacluster
###############################################################################

library(dplyr)

marker_medians <-
         aggregate(
                  expr_all,
                  by = list(Cluster = cell_meta20),
                  median
         )

head(marker_medians)

library(pheatmap)

mat <- as.matrix(marker_medians[,-1])

rownames(mat) <- paste0("MC", marker_medians$Cluster)

pheatmap(
         
         mat,
         
         scale = "row",
         
         clustering_method = "ward.D2",
         
         fontsize = 11,
         
         border_color = NA
         
)
