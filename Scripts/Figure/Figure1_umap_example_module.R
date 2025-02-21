library(celda)
library(singleCellTK)
library(lme4)
library(dplyr)
library(insight)
library(ggplot2)
library(here)

setwd(here::here("./Scripts/Figure"))

### 1. Load dataset
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))


red_umap <- readRDS('../../Mixed_model_output/corrected_module_UMAP.rds')
reducedDim(data_sce, 'corrected_umap') <- red_umap



### 2. visualize umap with study label
#### 2.1 sample
l <- 'sample'

fig_path <- '../..//Figure/Figure1'
fn <- file.path(fig_path, 'Figure1_sample_UMAP.pdf')
pdf(fn, width = 4, height = 4)
plotSCEDimReduceColData(data_sce, 
                        colorBy = l, reducedDimName = 'corrected_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = F, 
                        colorScale = NULL, title = l, titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  ) 
dev.off()

#### 2.2 stage
l <- 'Stage'
cols <- c("GGO" = "#F0F9E8", "I" = "#BAE4BC", "II" = "#7BCCC4", "III/IV" = "#2B8CBE")

fn <- file.path(fig_path, 'Figure1_Stage_UMAP.pdf')
pdf(fn, width = 5.5, height = 4)
plotSCEDimReduceColData(data_sce, 
                        colorBy = l, reducedDimName = 'corrected_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = F, 
                        colorScale = NULL, title = l, titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank()
    #legend.position="none"
  ) + 
  scale_colour_manual(values = cols)
dev.off()

#### 2.3 Histology
l <- 'Histology'
cols <- c("LUAD" = "#1B9E77", "LUSC" = "#D95F02", "NSCLC" = "#7570B3")

fn <- file.path(fig_path, 'Figure1_Histology_UMAP.pdf')
pdf(fn, width = 5.5, height = 4)
plotSCEDimReduceColData(data_sce, 
                        colorBy = l, reducedDimName = 'corrected_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = F, 
                        colorScale = NULL, title = l, titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank()
    #legend.position="none"
  ) + 
  scale_colour_manual(values = cols)
dev.off()


### 3. Visualize the corrected module score
histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]

#### load corrected module score
cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')

#### create sce object containing corrected module score instead of gene expression
sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                              colnames(data_sce)]),
                              rowData = NULL,
                              colData = colData(data_sce))
reducedDim(sce_m, 'corrected_umap') <- red_umap

#### 3.1 L10: NAPSA
l <- 'L10'
t <- 'L12 (NAPSA)'

#### 3.1 L12: NKX2-1
l <- 'L12'
t <- 'L12 (NKX2-1)'

#### 3.2 L16 KRT5
l <- 'L16'
t <- 'L16 (KRT5)'

#### 3.3 L64 MHC-II_1
l <- 'L64'
t <- 'L64 (MHC-II)'

#### 3.3 L65 MHC-II_2
l <- 'L65'
t <- 'L65 (MHC-II)'

#### 3.4 L35: Detox_4(AKR1C1)
l <- 'L35'
t <- 'L35 (Detox_4)'

#### 3.5 L55: VIM
l <- 'L55'
t <- 'L55 (VIM)'

#### 3.3 L59 ECM_4
l <- 'L59'
t <- 'L59 (ECM_4)'

#### 3.6 L43: TOP2A
l <- 'L43'
t <- 'L43 (TOP2A)'

l <- 'L15'
t <- 'L15 (TP63;NFE2L2)'

l <- 'L14'
t <- 'L14 (SOX2)'


#### 3.7 L36: Detox5(AKR1C2)
l <- 'L36'
t <- 'L36 (Detox_5)'


#### 3.7 L33: Detox2
l <- 'L33'
t <- 'L33 (Detox_2)'

#### 3.7 L39: DNA_replication_initiation
l <- 'L39'
t <- 'L39 (DNA_replication)'

#### 3.7 L41: DNA_replication_initiation
l <- 'L41'
t <- 'L41 (G1/S_Transition)'

col_scale <- c(0, 0.3, 0.6, 0.8, 0.995)
colors <- c('white', 'snow','orange','red', 'red3')

fig_path <- '../../Figure/Figure1'
fn <- file.path(fig_path, paste0('Figure1_module_',l,'_UMAP.tiff'))
tiff(fn, width = 400, height = 350)

plotSCEDimReduceFeatures(sce_m, feature=l, reducedDimName = 'corrected_umap', useAssay = 'module',
                         legendTitle = 'Corrected\nScore',
                         dotSize = 0.1, title = t, titleSize = 28, legendSize = 10, legendTitleSize = 15) + 
  scale_color_gradientn(values = col_scale, 
                       colours = colors,
                        na.value = "grey50") + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(5, 'mm'), legend.margin=margin(0,0,0,0)
  )
dev.off()


### 4. Example module heatmap
set.seed(2022)
data_sce_sub <- data_sce[, sample(colnames(data_sce), 30000)]


# moduleHeatmap(data_sce_sub, featureModule = l, useAssay = 'counts', altExpName = 'featureSubset', 
#               topCells = NULL, topFeatures = NULL, showFeatureNames = F)

for(l in c(65, 59, 10, 33, 39, 41)) { #c(15, 16, 35, 36, 43, 55, 64)
  fig_path <- '../../Figure/Figure1'
  fn <- file.path(fig_path, paste0('Figure1_module_L',l,'_Heatmap.tiff'))
  tiff(fn, width = 500, height = 350)
  print(
    moduleHeatmap(data_sce_sub, featureModule = l, useAssay = 'counts', altExpName = 'featureSubset', #data_sce_sub
                  topCells = NULL, topFeatures = NULL, showFeatureNames = F, showHeatmapLegend = T)
  )
  dev.off()
}



sessionInfo()
