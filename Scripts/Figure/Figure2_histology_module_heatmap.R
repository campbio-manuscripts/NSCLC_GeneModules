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


### 2. Creating figure 2A, marker modules of LUAD and LUSC
full_meta <- colData(data_sce)[, 1:14]
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               nNSCLC = sum(Histology == 'NSCLC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells,
                                                                               perNSCLC = nNSCLC / nCells)


### 2.1 define target modules
histo_m <- paste0('L', c(10:12, 
                         14:17))
target_modules <- histo_m
names(target_modules) <- rep(c('Histology'), c(length(histo_m)))

his <- c('LUAD', 'LUSC', 'NSCLC')
k <- unique(full_meta$celda_k) %>% sort()

cf <- 0.0001 ### ignore this parameter
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]




### define top annotation data
histology_anno <- feature_summary(feature='Histology', data_sce, module_expr, k, his, color_brewer = 'Dark2')
colAnno <- HeatmapAnnotation(Histology = anno_barplot(histology_anno[[1]], 
                                                      gp = gpar(fill = histology_anno[[2]]),
                                                      height = unit(2, "cm")),
                             annotation_legend_param = NULL)


row_order <- rownames(module_expr)
row_split <- names(target_modules)
col_name <- paste0('K', colnames(module_expr))
heatmap_legend_param <- list(title="Corrected\nModule\nScore", at = c(0, round(max(module_expr), 2)), 
                             labels = c(0, round(quantile(module_expr, 0.99), 2))) #c("low", "median", "high")




quantile(module_expr, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))
col_scale <- c(min(module_expr), 2.6, 5, 7, quantile(module_expr, 0.99)) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))


col_order <- NULL
hp<- Heatmap(module_expr, 
             heatmap_legend_param=heatmap_legend_param,
             column_title = "clustered dotplot", 
             col=col_fun, row_order = row_order,  #col_fun
             column_split = NULL, row_split = row_split,#sam_can$Stage
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 13),
             row_km = NULL,
             border = "black",
             top_annotation = colAnno, #
             show_column_dend = F, show_row_dend = F,
             height = nrow(module_expr)*unit(5, "mm"),
             column_order = NULL, column_labels = col_name)
draw(hp)








### 2.2 identify all modules associated with histology
nsclc_lm <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
nsclc_lm$Histology_fdr <- p.adjust(nsclc_lm$Histology_pval, 'fdr')

skip_m <- nsclc_lm$module[nsclc_lm$study_VarPer >= 0.3]
target_modules <- nsclc_lm$module[nsclc_lm$Histology_fdr < 0.05 & !nsclc_lm$module %in% skip_m]

target_modules <- target_modules[target_modules %in% paste0('L', 1:80)] #98


#### perform DE to determine target modules to include, and the row order
nsclc_m <- SingleCellExperiment(assay = list(module = cr_module[target_modules,  #lm_residual
                                                             colnames(data_sce)]),
                             rowData = NULL,
                             colData = colData(data_sce))
nsclc_m <- nsclc_m[, nsclc_m$Histology != 'NSCLC']

nsclc_m <- findMarkerDiffExp(nsclc_m,
                          useAssay = "module",
                          cluster = "Histology",
                          method = c('Limma'),
                          log2fcThreshold = 0.6) #0.6

hist_sum <- metadata(nsclc_m)$findMarker
luad_target <- hist_sum$Gene[hist_sum$Histology == 'LUAD' & hist_sum$clusterExprPerc > 0.3] #0.3
luad_order <- hist_sum %>% filter(Gene %in% luad_target) %>% arrange((clusterAveExpr)) %>% select(Gene) %>% unlist()
lusc_target <- hist_sum$Gene[hist_sum$Histology == 'LUSC' & hist_sum$clusterExprPerc > 0.3]
lusc_order <- hist_sum %>% filter(Gene %in% lusc_target) %>% arrange(desc(clusterAveExpr)) %>% select(Gene) %>% unlist()

target_modules <- c(luad_order, lusc_order)

his <- c('LUAD', 'LUSC', 'NSCLC')
k <- unique(full_meta$celda_k) %>% sort()

cf <- 0.001 ### ignore this parameter
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]




### define top annotation data
histology_anno <- feature_summary(feature='Histology', data_sce, module_expr, k, his, color_brewer = 'Dark2')
colAnno <- HeatmapAnnotation(Histology = anno_barplot(histology_anno[[1]], 
                                                      gp = gpar(fill = histology_anno[[2]]),
                                                      height = unit(2, "cm")),
                             annotation_legend_param = NULL)



### calculate the diversity index
#library(vegan)
#luad_diverse <- diversity(module_per[luad_order, ])
cf <- 0.05

luad_diverse <- apply(module_per[luad_order, ] > cf, 1, mean)
luad_order <- names(sort(luad_diverse))

#lusc_diverse <- diversity(module_per[lusc_order, ])
lusc_diverse <- apply(module_per[lusc_order, ] > cf, 1, mean)

lusc_order <- names(sort(lusc_diverse, decreasing = T))

row_order <- c(luad_order, lusc_order)
module_expr <- module_expr[row_order, ]

### relabel the rownames
### matching the module annotation as label
module_anno <- fread('../../Data/Module_annotation_Reorder_03272024.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))

rownames(module_expr) <- module_anno[match(rownames(module_expr), module_anno$Module), ] %>% pull(Label)
row_order <- module_anno[match(row_order, module_anno$Module), ] %>% pull(Label)

row_split <- rep(c('LUAD', 'LUSC'), c(length(luad_order), length(lusc_order)))
col_name <- paste0('K', colnames(module_expr))
heatmap_legend_param <- list(title="Corrected\nModule\nScore", at = c(0, round(max(module_expr), 2)), 
                             labels = c(0, round(quantile(module_expr, 0.99), 2))) #c("low", "median", "high")




quantile(module_expr, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))
col_scale <- c(min(module_expr), 2, 5, 7, quantile(module_expr, 0.99)) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))


col_order <- NULL
hp<- Heatmap(module_expr, 
             heatmap_legend_param=heatmap_legend_param,
             column_title = "", 
             col=col_fun, row_order = row_order,  #col_fun
             column_split = NULL, row_split = row_split,#sam_can$Stage
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             border = "black",
             top_annotation = colAnno, #
             show_column_dend = F, show_row_dend = F,
             #height = nrow(module_expr)*unit(2, "mm"),
             column_order = NULL, column_labels = col_name)


hist_anno_col <- histology_anno[[2]]
lgd_list = list(
  Legend(labels = c("LUAD",  'LUSC', "NSCLC"), title = "Histology", 
         legend_gp = gpar(fill = hist_anno_col))
)

fn <- file.path('../../Figure/Figure2', 'Figure2_Histology_Heatmap_10072024.pdf')
pdf(file = fn, width = 12, height = 11)
draw(hp, heatmap_legend_list = lgd_list, 
     heatmap_legend_side = "left", annotation_legend_side = "left") #bottom, left, Up, right
dev.off()
