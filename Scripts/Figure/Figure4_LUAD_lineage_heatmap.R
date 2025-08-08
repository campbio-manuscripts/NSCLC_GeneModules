library(celda)
library(singleCellTK)
library(lme4)
library(dplyr)
library(insight)
library(ggplot2)
library(tidyr)
library(circlize)
library(vegan)
library(RColorBrewer)
library(tidyr)
library(viridis)
library(Polychrome)
library(ComplexHeatmap)
library(data.table)
setwd(here::here("./Scripts/Figure"))

source('./Figure_utils.R')


### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))


#### load corrected module score
cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')
module_score <- t(cr_module)

#### module probability
x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs


### 2. Creating figure 4A, lineage modules of LUAD
full_meta <- colData(data_sce)
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               nNSCLC = sum(Histology == 'NSCLC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells,
                                                                               perNSCLC = nNSCLC / nCells)

### only focus on LUAD cluster
his <- c('LUAD')
k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57, 63)) ### perLUAD in 63 is 87.5, 12.5 NSCLC
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')


target_modules <- paste0('L',
                         c(12, 10, 11,
                           #28, 26, 
                           21, 20, 22,
                           19, 16
                         ))





### 3. prepare for the heatmap
cf <- 0.0001 ### based on MAST thres on cell probabilty, seems 0 is a reasonable cutoff
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]



#### define color and dot-plot function
col_scale <- c(min(module_expr), 2.5,  4.5, 7, quantile(module_expr, 0.99)) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))



#### using the eigen-vector of the target modules to sort the cluster order
pc1 <- eigenGene(expr = module_expr, module = paste0('L', 10:12)) # rownames(module_expr)
col_order <- names(pc1)[order(pc1, decreasing=F)]


#### define top annotation data
stage_anno <- feature_summary(feature='Stage', data_sce, module_expr, k, his, color_brewer = 'GnBu')
colAnno <- HeatmapAnnotation(
  Stage = anno_barplot(stage_anno[[1]],
                       gp = gpar(fill = stage_anno[[2]]),
                       height = unit(2, "cm")),
  annotation_legend_param = NULL)


#### other plotting param
rownames(module_expr) <- c(
  'L12:Adeno_3 (NKX2-1)',
  'L10:Adeno_1 (NAPSA)',
  'L11:Adeno_2 (RNASE1)',
  'L21:Club_2 (MUC5B)',
  'L20:Club_1 (SLPI)',
  'L22:Club_3 (SCGB1A1)',
  'L19:Goblet (MUC5AC)',
  'L16:Basal_1 (KRT5)'
)
row_order <- rownames(module_expr)
col_name <- paste0('K', colnames(module_expr))
heatmap_legend_param <- list(title="Corrected\nModule\nScore)", at = c(0, round(max(module_expr), 2)), 
                             labels = c(0, round(quantile(module_expr, 0.99), 2))) #c("low", "median", "high")


hp_scale <-  Heatmap(module_expr, 
                     heatmap_legend_param=heatmap_legend_param,
                     column_title = "clustered dotplot", 
                     col=col_fun, row_order = row_order,  #col_fun
                     column_split = NULL, row_split = NULL,#sam_can$Stage
                     cluster_rows = F, cluster_columns = F,
                     row_names_gp = gpar(fontsize = 13),
                     row_km = NULL,
                     border = "black",
                     top_annotation = colAnno, #
                     show_column_dend = F, show_row_dend = F,
                     height = nrow(module_expr)*unit(4.7, "mm"),
                     column_order = col_order, column_labels = col_name)

stage_anno_col <- stage_anno[[2]]
lgd_list = list(
  Legend(labels = c("GGO", "I", 'II', 'III/IV'), title = "Stage", 
         legend_gp = gpar(fill = stage_anno_col))
)


fn <- file.path('../../Figure/Figure4', 'Figure4_LUAD_lineage_Heatmap.pdf')
pdf(file = fn, width = 12, height = 5)
draw(hp_scale, heatmap_legend_list = lgd_list, 
     heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(0.1, 4, 0.1, 10), "mm")) #bottom, left, Up, right
dev.off()
sessionInfo()










