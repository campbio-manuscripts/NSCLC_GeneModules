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


### 2. subset the LUSC cluster
full_meta <- colData(data_sce)
full_meta$stage_bin <- ifelse(full_meta$Stage %in% c('GGO', 'I'), 'Early', 'Late')
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells)



his <- c('LUSC')
k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.6]
target_modules <- paste0('L', 
                         c(15, 16, 14,
                           17, 18, 
                           # 34, 37, 
                           32, 35:36,  
                           31, 51:52, 
                           47, 
                           82, 59,27, 70))


### 3. Prepare for plotting
cf <- 0.0001 ### ignore this parameter
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]


stage_anno <- feature_summary(feature='Stage', data_sce, module_expr, k, his, color_brewer = 'GnBu')
colAnno <- HeatmapAnnotation(
                             Stage = anno_barplot(stage_anno[[1]], 
                                                  gp = gpar(fill = stage_anno[[2]]), #fontsize = 4
                                                  height = unit(0.6, "cm"), 
                                                  width = unit(0.6, 'cm'), 
                                                  labels_gp = gpar(fontsize = 7)),
                             annotation_legend_param = NULL, 
                             annotation_name_gp= gpar(fontsize = 7)
                             )

row_order <- rownames(module_expr)
row_split <- names(target_modules)
col_name <- paste0('K', colnames(module_expr))
heatmap_legend_param <- list(title="Corrected\nModule\nScore", at = c(0, round(max(module_expr), 2)), 
                             labels = c(0, round(quantile(module_expr, 0.99), 2)), 
                             labels_gp = gpar(fontsize = 5), 
                             grid_height = unit(3.5, "mm"), 
                             grid_width = unit(3, "mm")) #c("low", "median", "high")




quantile(module_expr, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))
col_scale <- c(min(module_expr), 3, 4.2, 6, quantile(module_expr, 0.99)) #, 2.6
# col_scale <- c(min(module_expr), 2.5, 4.6, 6, quantile(module_expr, 0.99)) #, 2.6
# col_scale <- c(min(module_expr), 0.8, 2.5, 6, quantile(module_expr, 0.99)) #, 2.6

col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))


### matching the module annotation as label5
module_anno <- fread('../../Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))

rownames(module_expr) <- module_anno[match(rownames(module_expr), module_anno$Module), ] %>% pull(Label)


hp<- Heatmap(module_expr, 
             heatmap_legend_param=heatmap_legend_param,
             column_title = "", 
             col=col_fun, row_order = NULL,  #col_fun
             column_split = NULL, row_split = row_split,#sam_can$Stage
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             
             row_km = NULL,
             border = "black",
             top_annotation = colAnno, #
             show_column_dend = F, show_row_dend = F,
             height = nrow(module_expr)*unit(3, "mm"),
             width = nrow(module_expr)*unit(4, "mm"),
             column_order = NULL, column_labels = col_name)


hist_anno_col <- c("#1B9E77", "#7570B3")
stage_anno_col <- stage_anno[[2]]

lgd_list = list(
  Legend(labels = c("I", 'II', 'III/IV'), title = "Stage", 
         legend_gp = gpar(fill = stage_anno_col, fontsize = 7), 
         title_gp = gpar(fontsize = 7, fontface = "bold"), 
         labels_gp = gpar(fontsize = 6))
)

fn <- file.path('../../Figure/Figure6', 'Figure7_LUSC_Squamous_DE_Heatmap.pdf')
ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm") ### this is global setting. Be careful. 
ht_opt$HEATMAP_LEGEND_PADDING = unit(0.5, "cm")
pdf(file = fn, width = 4, height = 3)
draw(hp, heatmap_legend_list = lgd_list, 
     heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(0, 0, 0, 0), "mm"), #bottom, left, Up, right
     legend_title_gp = gpar(fontsize = 4, fontface = "bold")) 
dev.off()


### 4. Identify squamous-low and squamous-high clusters. And perform DE analysis
squa_sce <- data_sce[, data_sce$celda_k %in% k & data_sce$Histology %in% his]

histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]


squa_sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                   colnames(squa_sce)]),
                                   rowData = NULL,
                                   colData = colData(squa_sce)
                                   )
squa_sce_m$Squamous_status <- ifelse(squa_sce_m$celda_k %in% c('53', '55', '60'), 'Squamous_low', 'Squamous_high')

squa_sce_m <- singleCellTK::getUMAP(squa_sce_m, useAssay = 'module',
                      reducedDimName = 'corrected_UMAP', useReducedDim = NULL, 
                      sample = NULL, minDist = 0.2, seed = 2022, logNorm = F)

table(squa_sce_m$sample, squa_sce_m$Squamous_status)
table(squa_sce_m$sample[squa_sce_m$celda_k == 53]) #55, 60
table(squa_sce_m$sample[squa_sce_m$celda_k == 65])

table(squa_sce_m$celda_k[squa_sce_m$sample == 'Wu_et-P34']) %>% sort()

### all squamous_low are in Wu_et study
### most squamous-low cells are in Wu_et-P18 (which has both high and low cells).  Wu_et-P30 and Wu_et-P34 has other squamous low cells

s <- 'Wu_et-P18'
wu_p18_m <- squa_sce_m[, squa_sce_m$sample == 'Wu_et-P18']
wu_p18_m <- findMarkerDiffExp(wu_p18_m,
                           useAssay = "module",
                           cluster = "Squamous_status",
                           method = c('wilcox'),
                           log2fcThreshold = 0.6)

wu_p18_mark <- metadata(wu_p18_m)$findMarker
wu_p18_mark %>% group_by(Squamous_status) %>% arrange(desc(Log2_FC))


### or try include all squamous cells
squa_sce_m <- findMarkerDiffExp(squa_sce_m,
                              useAssay = "module",
                              cluster = "Squamous_status",
                              method = c('wilcox'),
                              log2fcThreshold = 0.6)

squa_mark <- metadata(squa_sce_m)$findMarker
squa_mark %>% group_by(Squamous_status) %>% arrange(desc(Log2_FC))



squa_low_mar <- squa_mark %>% filter(Squamous_status == 'Squamous_low', Log2_FC > 0.8, FDR <= 0.05)
squa_high_mar <- squa_mark %>% filter(Squamous_status == 'Squamous_high', Log2_FC > 0.8, FDR <= 0.05)
  
write.csv(rbind(squa_high_mar, squa_low_mar), 
          "../../Table/Squamous_HighVSlow_DE.csv")
#squamous_m <- target_modules

target_modules <- c(squa_low_mar$Gene, squa_high_mar$Gene)
names(target_modules) <- rep(c('Upregulated_Squamous-low', 'Upregulated_Squamous-high'), 
                             c(nrow(squa_low_mar), nrow(squa_high_mar)))

### use the function above to make the plot











## 5. Visualize the DE modules

fc <- 0.8
fdr <- 0.05
object <- wu_p18_m
marker <- wu_p18_mark
s <- unique(object$new_sample)



sig_marker <-  marker[marker$FDR <= fdr & marker$Log2_FC >= fc, ]
target_modules <- sig_marker$Gene
names(target_modules) <- sig_marker[match(target_modules, sig_marker$Gene), 'Squamous_status']


### mean of z score in the cluster
module_info <- scale(module_score[colnames(object), ])
his <- c('LUSC')
data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf,
                   clusters = unique(object$Squamous_status), his = his, data_sce = object,
                   module = target_modules, clu_col = 'Squamous_status')
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]



col_scale <- c(-2, 0.1, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))

# quantile(module_expr, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))
# col_scale <- c(min(module_expr), 2.6, 4.6, 6, quantile(module_expr, 0.99)) #, 2.6
# col_fun = circlize::colorRamp2(col_scale,  c('white', 'white','orange','red', 'red3'))

#col_order <- c("Adeno-like (K15/56)",  "Squamous-like (K64)", "Club-like (K48/K62)", "Mucinous (K61)")
#row_order <- c("Adeno-like", "Squamous-like", "Club", "Mucinous")
#row_split <- factor(names(target_modules), levels = row_order)

row_split <- names(target_modules)
hm_expr <- module_expr
rownames <- module_anno[match(rownames(hm_expr), module_anno$Module), 'Label']
rownames(hm_expr) <- rownames$Label

cluster_col <- setNames(c("#1F78B4", "#FCBE6E"),
                        c("Squamous_low",  "Squamous_high")
)
topAnno <- HeatmapAnnotation('Cluster' = colnames(hm_expr),
                             col = list(Cluster = cluster_col))
col_order <- NULL

de_hm <- Heatmap(hm_expr, 
                 name = 'Relative Expression',
                 column_title = s,
                 row_title_rot  =90, 
                 col=col_fun, row_order = NULL, column_order = col_order,
                 column_split = NULL, row_split = row_split,
                 cluster_rows = F, cluster_columns = F, show_column_names = F,
                 row_gap = unit(2, "mm"),
                 row_names_gp = gpar(fontsize = 8),
                 row_km = NULL,
                 border = "black",
                 top_annotation = topAnno, #
                 # height = nrow(module_expr)*unit(4.7, "mm"),
                 # width = ncol(module_expr)*unit(9, "mm"),
                 heatmap_legend_param = list(
                   legend_direction = "horizontal", 
                   legend_width = unit(3, "cm")),
                 show_column_dend = F, show_row_dend = F)

draw(de_hm,
     heatmap_legend_side = "bottom", annotation_legend_side = "right")



##### compare the Squamous DE genes with what we saw in CPTAC data
mod_int <- intersect(sig_l, target_modules)
mod_int[!mod_int %in% paste0('L', 14:18)]

mod_int[!mod_int %in% paste0('L', 14:18)][mod_int[!mod_int %in% paste0('L', 14:18)] %in% squa_low_mar$Gene]


## 6. Create umap for Wu-P18 to show the heterogeneity
wu_p18_m <- squa_sce_m[, squa_sce_m$sample == 'Wu_et-P18']

### remove modules with high %study
histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]
wu_p18_m <- wu_p18_m[!rownames(wu_p18_m) %in% skip_m, ]

### generate umap
wu_p18_m <- getUMAP(wu_p18_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)
l <- 'Squamous_status'
cluster_col <- setNames(c("#1F78B4", "#FCBE6E"),
                        c("Squamous_low",  "Squamous_high")
)

p18_umap <- plotSCEDimReduceColData(wu_p18_m, 
                              colorBy = l, reducedDimName = 'module_UMAP', 
                              clusterLabelSize = 4, dotSize = 0.3,
                              colorScale = cluster_col, 
                              title = s) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )

pdf('../../Figure/Figure7/Wu_P18_umap.pdf', width = 3, height = 3)
print(p18_umap)
dev.off()



## 7. For DE between all squamous low vs high, check a few genes expression

#### several modules only high in squamous low cell clusters

lusc_pl <- lapply(paste0('L', c(82, 25, 90, 93, 31, 17)), function(l) {
  col_scale <- c(0, 0.3, 0.6, 0.8, 0.995)
  colors <- c('grey88', 'snow','orange','red', 'red3')
  
  
  
  plotSCEDimReduceFeatures(squa_sce_m, feature=l, reducedDimName = 'corrected_UMAP', useAssay = 'module',
                           #groupBy = 'sample',
                           legendTitle = 'Corrected\nScore',
                           dotSize = 0.1, title = l, titleSize = 10, legendSize = 6, legendTitleSize = 6,
                           combinePlot = 'none') + 
    scale_color_gradientn(values = col_scale, 
                          colours = colors,
                          na.value = "grey50") + 
    theme(
      axis.title = element_blank(),         # Change both x and y axis titles
      axis.text=element_blank(),
      axis.ticks = element_blank(),
      legend.key.width = unit(2, 'mm'), legend.key.height =unit(2.5, 'mm') ,legend.margin=margin(0,0,0,0), 
      legend.box.spacing = margin(2, 2, 2, 2)
    )
})

cowplot::plot_grid(plotlist = lusc_pl, nrow = 3)
