###############################################################################
# MUCOSAL HIGH-DIMENSIONAL ANALYSIS
# SCRIPT 01: IMPORT, QC & MEMORY-EFFICIENT PREPARATION
###############################################################################

rm(list = ls())
graphics.off()
gc()

###############################################################################
# Load packages
###############################################################################

library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(tidyverse)
library(stringr)
library(tibble)

###############################################################################
# Working directory
###############################################################################

setwd("C:/Users/LENOVO/OneDrive - Malawi-Liverpool Wellcome Research Programme/Desktop/LAB WORK/Exported FCS Files for Mucosal Analysis/High Dimensionality reductionl-2026/DownSample")

###############################################################################
# Memory helper
###############################################################################

mem_report <- function(step){
         
         cat("\n=================================================\n")
         cat(step,"\n")
         print(gc())
         cat("=================================================\n\n")
         
}

mem_report("Start")

###############################################################################
# Locate FCS files
###############################################################################

fcs_files <- list.files(
         pattern = "\\.fcs$",
         full.names = TRUE,
         ignore.case = TRUE
)

cat(length(fcs_files),"FCS files detected.\n")

###############################################################################
# Remove PBMC
###############################################################################

fcs_files <- fcs_files[
         !grepl("PBMC", fcs_files, ignore.case = TRUE)
]

cat(length(fcs_files),"Mucosal samples retained.\n")

###############################################################################
# Read FlowSet
###############################################################################

fs <- read.flowSet(
         files = fcs_files,
         truncate_max_range = FALSE
)

mem_report("FlowSet imported")

###############################################################################
# Event counts
###############################################################################

event_counts <- fsApply(fs, nrow)

summary(event_counts)

sort(event_counts)

###############################################################################
# Sample metadata
###############################################################################

metadata <- tibble(
         
         file_name = basename(fcs_files),
         
         sample_id = basename(fcs_files),
         
         patient = str_extract(file_name,"DKP\\d+[A-Z]?"),
         
         tissue = case_when(
                  
                  str_detect(file_name,"Nasal Swab") ~ "Nasal Swab",
                  
                  str_detect(file_name,"Nasal Scrape") ~ "Nasal Scrape",
                  
                  str_detect(file_name,"Cervical Scrape") ~ "Cervical Scrape",
                  
                  str_detect(file_name,"Cervical Cytobrush") ~ "Cervical Cytobrush",
                  
                  TRUE ~ "Unknown"
                  
         )
         
)

metadata

###############################################################################
# Keep ONE copy of fluorescence channels
###############################################################################

keep_channels <- c(
         
         "FSC-A",
         "FSC-H",
         "FSC-W",
         
         "SSC-A",
         "SSC-H",
         "SSC-W",
         
         "FJComp-BV785-A-1",               # CD3
         "FJComp-APC-Cy7-A-1",             # CD4
         "FJComp-AF700-A-1",               # CD8
         "FJComp-BV750-A-1",               # CD19
         "FJComp-BV421-A-1",               # CD14
         "FJComp-BV650-A-1",               # CD16
         "FJComp-PE-A-1",                  # CD66b
         "FJComp-APC-A-1",                 # CD56
         "FJComp-BV605-A-1",               # CD69
         "FJComp-PE-Dazzle594-A-1",        # CD103
         "FJComp-FITC-A-1",                # CD45RO
         "FJComp-PE-Cy7-A-1",              # CCR7
         "FJComp-PerCP-Cy5.5-A-1",         # CD45
         "FJComp-LiveDeadFixableAqua-A-1", # Live/Dead
         
         "TIME-1"
         
)

fs <- fs[, keep_channels]

mem_report("Duplicate channels removed")

###############################################################################
# Panel annotation
###############################################################################

panel <- tibble(
         
         fcs_colname = c(
                  
                  "FJComp-BV785-A-1",
                  "FJComp-APC-Cy7-A-1",
                  "FJComp-AF700-A-1",
                  "FJComp-BV750-A-1",
                  "FJComp-BV421-A-1",
                  "FJComp-BV650-A-1",
                  "FJComp-PE-A-1",
                  "FJComp-APC-A-1",
                  "FJComp-BV605-A-1",
                  "FJComp-PE-Dazzle594-A-1",
                  "FJComp-FITC-A-1",
                  "FJComp-PE-Cy7-A-1",
                  "FJComp-PerCP-Cy5.5-A-1",
                  "FJComp-LiveDeadFixableAqua-A-1"
                  
         ),
         
         antigen = c(
                  
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
                  "CCR7",
                  "CD45",
                  "LiveDead"
                  
         ),
         
         marker_class = c(
                  
                  rep("type",12),
                  "none",
                  "none"
                  
         )
         
)

###############################################################################
# Save checkpoint
###############################################################################

dir.create("Results",showWarnings = FALSE)

saveRDS(fs,"Results/fs_clean.rds")

saveRDS(metadata,"Results/metadata.rds")

saveRDS(panel,"Results/panel.rds")

###############################################################################
# Release RAM
###############################################################################

rm(event_counts)

gc()

mem_report("Script 01 finished")

#Number of cells
cell_summary %>%
         dplyr::group_by(Compartment) %>%
         dplyr::summarise(
                  `>=10k` = sum(cells >= 10000),
                  `>=25k` = sum(cells >= 25000),
                  `>=50k` = sum(cells >= 50000),
                  Total = dplyr::n()
         )




###############################################################################
# MUCOSAL HIGH-DIMENSIONAL ANALYSIS
# SCRIPT 02: FILTER & DOWNSAMPLE
#
# PURPOSE
#   1. Load cleaned FlowSet
#   2. Remove samples with <10,000 cells
#   3. Downsample remaining samples to exactly 10,000 cells
#   4. Save checkpoint for CATALYST
###############################################################################

###############################################################################
# Clean workspace
###############################################################################

rm(list = ls())
graphics.off()
gc()

###############################################################################
# Load packages
###############################################################################

library(flowCore)
library(tidyverse)

###############################################################################
# Memory helper
###############################################################################

mem_report <- function(step){
         
         cat("\n=====================================================\n")
         cat(step, "\n")
         print(gc())
         cat("=====================================================\n\n")
         
}

mem_report("Start Script 02")

###############################################################################
# Load checkpoint
###############################################################################

fs <- readRDS("Results/fs_clean.rds")
metadata <- readRDS("Results/metadata.rds")
panel <- readRDS("Results/panel.rds")

mem_report("Objects loaded")

###############################################################################
# Count cells
###############################################################################

event_counts <- fsApply(fs, nrow)

cell_summary <- tibble(
         
         sample = sampleNames(fs),
         
         cells = event_counts
         
)

###############################################################################
# Add metadata
###############################################################################

cell_summary <- bind_cols(
         
         metadata,
         
         cells = event_counts
         
)

###############################################################################
# Save original cell counts
###############################################################################

write.csv(
         
         cell_summary,
         
         "Results/Cell_Counts_Before_Filtering.csv",
         
         row.names = FALSE
         
)

###############################################################################
# Keep samples >=10,000 cells
###############################################################################

keep <- event_counts >= 10000

cat("\n")
cat(sum(keep), "samples retained\n")
cat(sum(!keep), "samples removed\n")
cat("\n")

###############################################################################
# Filter FlowSet
###############################################################################

fs <- fs[keep]

###############################################################################
# Filter metadata
###############################################################################

metadata <- metadata[keep, ]

###############################################################################
# Release RAM
###############################################################################

rm(event_counts)
rm(keep)
rm(cell_summary)

gc()

mem_report("Filtering complete")

###############################################################################
# Downsample
###############################################################################

set.seed(2026)

for(i in seq_len(length(fs))){
         
         ff <- fs[[i]]
         
         expr <- exprs(ff)
         
         idx <- sample.int(
                  
                  n = nrow(expr),
                  
                  size = 10000,
                  
                  replace = FALSE
                  
         )
         
         exprs(ff) <- expr[idx, , drop = FALSE]
         
         fs[[i]] <- ff
         
         rm(ff)
         rm(expr)
         rm(idx)
         
         if(i %% 5 == 0){
                  
                  gc(verbose = FALSE)
                  
         }
         
}

mem_report("Downsampling complete")

###############################################################################
# Verify
###############################################################################

new_counts <- fsApply(fs, nrow)

print(new_counts)

summary(new_counts)

stopifnot(all(new_counts == 10000))

###############################################################################
# Save retained sample information
###############################################################################

metadata$Cells_After_Downsampling <- 10000

write.csv(
         
         metadata,
         
         "Results/Metadata_10k.csv",
         
         row.names = FALSE
         
)

###############################################################################
# Save FlowSet
###############################################################################

saveRDS(
         
         fs,
         
         "Results/fs_10k.rds"
         
)

###############################################################################
# Save panel
###############################################################################

saveRDS(
         
         panel,
         
         "Results/panel.rds"
         
)

###############################################################################
# Release RAM
###############################################################################

rm(new_counts)

gc()

mem_report("Script 02 Finished")
