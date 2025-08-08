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


### 2. Creating figure 3B, stage associated modules of LUAD
full_meta <- colData(data_sce)[, 1:14]
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               nNSCLC = sum(Histology == 'NSCLC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells,
                                                                               perNSCLC = nNSCLC / nCells)







ip <- '../../Mixed_model_output/'
luad_lm <- fread(file.path(ip, 'model_summary_Stage.csv'))
luad_lm$Stage_fdr <- p.adjust(luad_lm$Stage_pval, 'fdr')


stage_m <- luad_lm %>% filter(luad_lm$Stage_fdr <= 0.05, abs(luad_lm$Stage_coef) >= 0.55) %>%  #0.575, 0.6, 0.5
  arrange((Stage_coef)) %>% dplyr::select(module) %>% unlist()

#target_modules <- stage_m[!stage_m %in% c('L101', 'L102')]
target_modules <- stage_m
names(target_modules) <- ifelse(luad_lm[match(target_modules, luad_lm$module), 'Stage_coef'] > 0,
                                'Upregulated in late stage',
                                'Upregulated in early stage')



### only focus on LUAD cluster
his <- c('LUAD')
k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57, 63)) 
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')



### regenerate the module expression matrix
cf <- 0.0001 
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]


histology_anno <- feature_summary(feature='Histology', data_sce, module_expr, k, his, color_brewer = 'Dark2')
stage_anno <- feature_summary(feature='Stage', data_sce, module_expr, k, his, color_brewer = 'GnBu')
colAnno <- HeatmapAnnotation(Histology = anno_barplot(histology_anno[[1]], 
                                                      gp = gpar(fill = c("#1B9E77", "#7570B3")), #histology_anno[[2]]
                                                      height = unit(2, "cm")),
                             Stage = anno_barplot(stage_anno[[1]], 
                                                  gp = gpar(fill = stage_anno[[2]]),
                                                  height = unit(2, "cm")),
                             annotation_legend_param = NULL)
heatmap_legend_param <- list(title="Corrected\nModule\nScore", at = c(0, round(max(module_expr), 2)), 
                             labels = c(0, round(quantile(module_expr, 0.99), 2))) #c("low", "median", "high")

quantile(module_expr, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))
col_scale <- c(min(module_expr), 3, 4.6, 6, quantile(module_expr, 0.99)) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))

pc1 <- eigenGene(expr = module_expr, module = target_modules)
col_order <- names(pc1)[order(pc1, decreasing=F)]
col_name <- paste0('K', colnames(module_expr))

### matching the module annotation as label5
module_anno <- fread('../../Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))

rownames(module_expr) <- module_anno[match(rownames(module_expr), module_anno$Module), ] %>% pull(Label)
 #module_anno %>% filter(Module %in% rownames(module_expr)) %>% arrange(match(Module, rownames(module_expr))) %>% pull(Label)


stage_hp <- Heatmap(module_expr, 
             heatmap_legend_param=heatmap_legend_param,
             column_title = NULL, 
             col=col_fun, row_order = NULL,  #col_fun
             column_split = NULL, row_split = names(target_modules),#sam_can$Stage
             cluster_rows = T, cluster_columns = F,cluster_row_slices = F,
             row_names_gp = gpar(fontsize = 12),
             border = "black",
             top_annotation = colAnno, #
             show_column_dend = F, show_row_dend = F,
             height = nrow(module_expr)*unit(4.7, "mm"),
             column_order = col_order, column_labels = col_name)



### add customized legend color
hist_anno_col <- c("#1B9E77", "#7570B3")
stage_anno_col <- stage_anno[[2]]

lgd_list = list(
  Legend(labels = c("LUAD", "NSCLC"), title = "Histology", 
         legend_gp = gpar(fill = hist_anno_col)),
  Legend(labels = c("GGO", "I", 'II', 'III/IV'), title = "Stage", 
         legend_gp = gpar(fill = stage_anno_col))
  )


fn <- file.path('../../Figure/Figure4', 'Figure3_Stage_Heatmap_fc0.55.pdf') #Figure3_Stage_Heatmap.pdf
pdf(file = fn, width = 12, height = 12)#11
draw(stage_hp, heatmap_legend_list = lgd_list, 
     heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(0.1, 2, 0.1, 10), "mm")) #bottom, left, Up, right
dev.off()



sessionInfo()