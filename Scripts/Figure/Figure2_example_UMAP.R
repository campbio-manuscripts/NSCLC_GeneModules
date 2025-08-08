#### generate example UMAP for figure 2
library(celda)
library(singleCellTK)
library(lme4)
library(dplyr)
library(insight)
library(ggplot2)

setwd(here::here("./Scripts/Figure"))


### 1. Load dataset
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))


### 2. subset to following samples
samples <- c("Sinjab_et-P2-1-Epcam-pos", "Sinjab_et-P3-1-Epcam-pos", 
             "Sinjab_et-P5-1-Epcam-pos", "Qian_et-LC_3", "Qian_et-LC_4", "Qian_et-LC_6")

sub_sce <- data_sce[, data_sce$sample %in% samples]


### 3. Re-generated the example UMAP
sub_sce <- runNormalization(sub_sce, normalizationMethod = "LogNormalize", outAssayName = 'LogNormalized')
sub_sce <- runUMAP(sub_sce, useAssay = 'LogNormalized', useReducedDim = NULL)


### 4. Visualize
l <- 'sample'

fig_path <- '../../Figure/Figure2'
fn <- file.path(fig_path, 'Figure2_sample_UMAP.pdf')
pdf(fn, width = 3, height = 3)
plotSCEDimReduceColData(sub_sce, 
                        colorBy = l, reducedDimName = 'UMAP', 
                        clusterLabelSize = 2, dotSize = 0.1, labelClusters = T, 
                        colorScale = NULL) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  ) 
dev.off()