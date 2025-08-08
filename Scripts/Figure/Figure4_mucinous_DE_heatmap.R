library(celda)
library(singleCellTK)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
setwd(here::here("./Scripts/Figure"))


source('../../Scripts/Figure/Figure_utils.R')

### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

#### load module annotation
mod_anno <- fread('../../Data/Module_annotation_Reorder_1008.csv')
mod_anno$Label <- paste(mod_anno$Module, mod_anno$Name, sep = ':')


### 2. calculate module probability and corrected module score
x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs

#### load corrected module score
cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')
module_score <- t(cr_module)


### 3. selecting K that has high expression of Mucinous / club modules
full_meta <- colData(data_sce)
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               nNSCLC = sum(Histology == 'NSCLC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells,
                                                                               perNSCLC = nNSCLC / nCells)

#### only focus on LUAD cluster
his <- c('LUAD')
k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57, 63)) ### perLUAD in 63 is 87.5, 12.5 NSCLC
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')



### 4.  Select samples and cell clusters that has high expression of club and globlet modules
target_modules <- c('L22', 'L19')
cf <- 0.005 ### based on MAST thres on cell probabilty, seems 0 is a reasonable cutoff
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
module_per <- data[[3]]

club_k <- colnames(module_per)[module_per['L22', ] > 0.3]
mucinous_k <- colnames(module_per)[module_per['L19', ] > 0.3]
ks <- c(1, club_k, mucinous_k, 80) ### Cluster 1: oridinary LUAD cluster. Cluster 80: ordinary squamous cluster. 


###### further select samples and the clusters
sub_sce <- data_sce[, data_sce$celda_k %in% ks]
sub_sce$celda_k <- as.character(sub_sce$celda_k)
selected_sample <- names(table(sub_sce$sample))[table(sub_sce$sample) > 50] ### keep sample with > 10 cells in these clusters

sub_sce <- data_sce[, data_sce$sample %in% selected_sample]
select_k <- names(table(sub_sce$celda_k))[table(sub_sce$celda_k) > 20] ### keep cluster with > 10 cells
sub_sce <- sub_sce[, sub_sce$celda_k %in% select_k]
sub_sce$celda_k <- as.character(sub_sce$celda_k)



###### skip modules that has high %study
histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]



### 2. For Sinjab_P2
####### 2.1 only include P2
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Sinjab_et-P2-1-Epcam-pos')] 

# use original modules score to perform DE analysis
p2_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                kept_cid]),
                             rowData = NULL,
                             colData = colData(sub_sce)[kept_cid, ])

p2_m$cluster <- p2_m$celda_k
p2_m$cluster[p2_m$cluster %in% c('1', '16', '62', '48')] <- 'Club-like (K48/K62)'
p2_m$cluster[p2_m$cluster %in% c('59', '61', '63')] <- 'Mucinous (K61)'
p2_m$cluster[p2_m$cluster %in% c('64')] <- 'Intermediate (K64)'
p2_m$cluster[p2_m$cluster %in% c('3','4','15','22','56') & p2_m$sample == 'Sinjab_et-P2-1-Epcam-pos'] <- 'Adeno-like (K15/56)'

p2_m <- findMarkerDiffExp(p2_m,
                          useAssay = "module",
                          cluster = "cluster",
                          method = c('wilcox'),
                          log2fcThreshold = 0.6)

p2_marker <- metadata(p2_m)$findMarker
p2_marker <- p2_marker %>% mutate(
  cluster = case_when(
    cluster == 'Adeno-like (K15/56)' ~ 'Adeno-like',
    cluster == 'Club-like (K48/K62)' ~ 'Club',
    cluster == 'Mucinous (K61)' ~ 'Mucinous',
    cluster == 'Intermediate (K64)' ~ 'Intermediate'
    
  )
)

###### visualization
fc <- 1.5
fdr <- 0.05
object <- p2_m
s <- unique(object$new_sample)
marker <- p2_marker

module_info <- scale(module_score[colnames(object), ])
his <- c('LUAD', 'NSCLC')

sig_marker <-  marker[marker$FDR <= fdr & marker$Log2_FC >= fc, ]

###### select modules only DE in one cluster
dup_mods <-  marker %>% group_by(Gene) %>% summarise(numClu = length(unique(cluster))) %>% filter(numClu != 1)
target_modules <- sig_marker$Gene[!sig_marker$Gene %in% dup_mods$Gene]

###### adding some important modules into this heatmap
target_modules <- c(target_modules, 'L19')
names(target_modules) <- sig_marker[match(target_modules, sig_marker$Gene), 'cluster']

data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf, 
                   clusters = unique(object$cluster), his = his, data_sce = object,
                   module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]


col_scale <- c(-2, 0.1, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))

col_order <- c("Adeno-like (K15/56)",  "Intermediate (K64)", "Club-like (K48/K62)", "Mucinous (K61)")
row_order <- c("Adeno-like", "Intermediate", "Club", "Mucinous")
row_split <- factor(names(target_modules), levels = row_order)

hm_expr <- module_expr
rownames <- mod_anno[match(rownames(hm_expr), mod_anno$Module), 'Label']
rownames(hm_expr) <- rownames$Label


cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#A65528", "#FCBE6E"),
                        c("Intermediate (K64)", "Mucinous (K61)", "Club-like (K48/K62)", "Adeno-like (K15/56)")
)

topAnno <- HeatmapAnnotation('Cluster' = colnames(hm_expr),
                             col = list(Cluster = cluster_col[colnames(hm_expr)]))

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
                 #height = nrow(module_expr)*unit(3.5, "mm"),
                 height = unit(115, 'mm'), 
                 width = ncol(module_expr)*unit(12, "mm"),
                 heatmap_legend_param = list(
                   legend_direction = "horizontal", 
                   legend_width = unit(3, "cm")),
                 show_column_dend = F, show_row_dend = F)

pdf('../../Figure/Figure4/Sinjab_P2_DE_heatmap.pdf')
draw(de_hm,
     heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()








### 3. For p30 
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Bischoff_et-p030t')]  


###### adding ordinary luad clusters from the same study Bischoff_et
Bischoff_meta <- full_meta[full_meta$study == 'Bischoff_et', ]
table(Bischoff_meta$sample, Bischoff_meta$celda_k)

###### pick-up the cluster of oridinary LUAD cluster K29
kept_cid <- c(kept_cid, 
              colnames(data_sce)[data_sce$study == 'Bischoff_et' & data_sce$celda_k == 11])

p30_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                 kept_cid]),
                              rowData = NULL,
                              colData = colData(data_sce)[kept_cid, ]) ##sub_sce


p30_m$cluster <- as.character(p30_m$celda_k)
p30_m$cluster[p30_m$cluster %in% c('59', '61')] <- 'Mucinous (K59)'
p30_m$cluster[p30_m$cluster %in% c('63', '31')] <- 'Intermediate (K63)'
p30_m$cluster[p30_m$cluster %in% c('11')] <- 'Adeno (K11)'

p30_m <- findMarkerDiffExp(p30_m,
                           useAssay = "module",
                           cluster = "cluster",
                           method = c('wilcox'),
                           log2fcThreshold = 0.6)

p30_marker <- metadata(p30_m)$findMarker
p30_marker <- p30_marker %>% mutate(
  cluster = case_when(
    cluster == 'Mucinous (K59)' ~ 'Mucinous',
    cluster == 'Intermediate (K63)' ~ 'Intermediate',
    cluster == 'Adeno (K11)' ~ 'Adeno'
  )
)

####### visualization
fc <- 1.5
fdr <- 0.05
object <- p30_m
s <- "Bischoff_et-p030t"
marker <- p30_marker

module_info <- scale(module_score[colnames(object), ])
his <- c('LUAD', 'NSCLC')

sig_marker <-  marker[marker$FDR <= fdr & marker$Log2_FC >= fc, ]

###### select modules only DE in one cluster
dup_mods <-  marker %>% group_by(Gene) %>% summarise(numClu = length(unique(cluster))) %>% filter(numClu != 1)
target_modules <- sig_marker$Gene[!sig_marker$Gene %in% dup_mods$Gene]
names(target_modules) <- sig_marker[match(target_modules, sig_marker$Gene), 'cluster']


###### adding some important modules into this heatmap
if (!'L19' %in% target_modules) {
  target_modules <- c(target_modules, 'L19')
  names(target_modules)[target_modules == 'L19'] <- 'Mucinous'
}



data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf, 
                   clusters = unique(object$cluster), his = his, data_sce = object,
                   module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]

col_scale <- c(-2, 0.1, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))


col_order <- c("Adeno (K11)", "Intermediate (K63)",  "Mucinous (K59)" )
row_order <- c("Adeno", "Intermediate", "Mucinous")

row_split <- factor(names(target_modules), levels = row_order)

hm_expr <- module_expr
rownames <- mod_anno[match(rownames(hm_expr), mod_anno$Module), 'Label']
rownames(hm_expr) <- rownames$Label

cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#FCBE6E"),
                        c("Intermediate (K63)", "Mucinous (K59)", "Adeno (K11)")
)

topAnno <- HeatmapAnnotation('Cluster' = colnames(hm_expr),
                             col = list(Cluster = cluster_col[colnames(hm_expr)]))

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
                 #height = nrow(module_expr)*unit(3.5, "mm"),
                 height = unit(115, 'mm'), 
                 width = ncol(module_expr)*unit(12, "mm"),
                 heatmap_legend_param = list(
                   legend_direction = "horizontal", 
                   legend_width = unit(3, "cm")),
                 show_column_dend = F, show_row_dend = F)

pdf('../../Figure/Figure4/Bischofft_p030t_DE_heatmap.pdf')
draw(de_hm,
     heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()



### 4. For p33
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Bischoff_et-p033t')]  


###### adding ordinary luad clusters from the same study Bischoff_et
Bischoff_meta <- full_meta[full_meta$study == 'Bischoff_et', ]
table(Bischoff_meta$sample, Bischoff_meta$celda_k)

###### pick-up the cluster of oridinary LUAD cluster K29
kept_cid <- c(kept_cid, 
              colnames(data_sce)[data_sce$study == 'Bischoff_et' & data_sce$celda_k == 11])

p30_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                              kept_cid]),
                              rowData = NULL,
                              colData = colData(data_sce)[kept_cid, ]) ##sub_sce


p30_m$cluster <- as.character(p30_m$celda_k)
p30_m$cluster[p30_m$cluster %in% c('59', '61')] <- 'Mucinous (K59)'
p30_m$cluster[p30_m$cluster %in% c('63', '31')] <- 'Intermediate (K63)'
p30_m$cluster[p30_m$cluster %in% c('11')] <- 'Adeno (K11)'

p30_m <- findMarkerDiffExp(p30_m,
                           useAssay = "module",
                           cluster = "cluster",
                           method = c('wilcox'),
                           log2fcThreshold = 0.6)

p30_marker <- metadata(p30_m)$findMarker
p30_marker <- p30_marker %>% mutate(
  cluster = case_when(
    cluster == 'Mucinous (K59)' ~ 'Mucinous',
    cluster == 'Intermediate (K63)' ~ 'Intermediate',
    cluster == 'Adeno (K11)' ~ 'Adeno'
  )
)

####### visualization
fc <- 1.5
fdr <- 0.05
object <- p30_m
s <- "Bischoff_et-p033t"
marker <- p30_marker

module_info <- scale(module_score[colnames(object), ])
his <- c('LUAD', 'NSCLC')

sig_marker <-  marker[marker$FDR <= fdr & marker$Log2_FC >= fc, ]

###### select modules only DE in one cluster
dup_mods <-  marker %>% group_by(Gene) %>% summarise(numClu = length(unique(cluster))) %>% filter(numClu != 1)
target_modules <- sig_marker$Gene[!sig_marker$Gene %in% dup_mods$Gene]
names(target_modules) <- sig_marker[match(target_modules, sig_marker$Gene), 'cluster']


###### adding some important modules into this heatmap
if (!'L19' %in% target_modules) {
  target_modules <- c(target_modules, 'L19')
  names(target_modules)[target_modules == 'L19'] <- 'Mucinous'
}



data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf, 
                   clusters = unique(object$cluster), his = his, data_sce = object,
                   module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]

col_scale <- c(-2, 0.1, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))


col_order <- c("Adeno (K11)", "Intermediate (K63)",  "Mucinous (K59)" )
row_order <- c("Adeno", "Intermediate", "Mucinous")

row_split <- factor(names(target_modules), levels = row_order)

hm_expr <- module_expr
rownames <- mod_anno[match(rownames(hm_expr), mod_anno$Module), 'Label']
rownames(hm_expr) <- rownames$Label

cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#FCBE6E"),
                        c("Intermediate (K63)", "Mucinous (K59)", "Adeno (K11)")
)

topAnno <- HeatmapAnnotation('Cluster' = colnames(hm_expr),
                             col = list(Cluster = cluster_col[colnames(hm_expr)]))

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
                 #height = nrow(module_expr)*unit(3.5, "mm"),
                 height = unit(115, 'mm'), 
                 width = ncol(module_expr)*unit(12, "mm"),
                 heatmap_legend_param = list(
                   legend_direction = "horizontal", 
                   legend_width = unit(3, "cm")),
                 show_column_dend = F, show_row_dend = F)

pdf('../../Figure/Figure4/Bischofft_p033t_DE_heatmap.pdf')
draw(de_hm,
     heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()





### 5. For Zilionis_et-Zilionis_Klein_2019_p3
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Zilionis_et-p3t1', 'Zilionis_et-p3t2', 'Zilionis_et-p3t3')]  


###### adding ordinary luad clusters from the same study Zilionis_et
Zilionis_meta <- full_meta[full_meta$study == 'Zilionis_et', ]
table(Zilionis_meta$sample, Zilionis_meta$celda_k)

###### pick-up the cluster of oridinary LUAD cluster K29
kept_cid <- c(kept_cid, 
              colnames(data_sce)[data_sce$study == 'Zilionis_et' & data_sce$celda_k %in% c(6, 11, 15)]) ### mostly from p6t1/t2


p3_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                              kept_cid]),
                              rowData = NULL,
                              colData = colData(data_sce)[kept_cid, ]) ##sub_sce


p3_m$cluster <- as.character(p3_m$celda_k)
p3_m$cluster[p3_m$cluster %in% c('59', '61')] <- 'Mucinous (K59/K61)'
p3_m$cluster[p3_m$cluster %in% c('6', '11', '15')] <- 'Adeno (K11/K6/K15)'

p3_m <- findMarkerDiffExp(p3_m,
                           useAssay = "module",
                           cluster = "cluster",
                           method = c('wilcox'),
                           log2fcThreshold = 0.6)

p3_marker <- metadata(p3_m)$findMarker
p3_marker <- p3_marker %>% mutate(
  cluster = case_when(
    cluster == 'Mucinous (K59/K61)' ~ 'Mucinous',
    cluster == 'Adeno (K11/K6/K15)' ~ 'Adeno'
  )
)

####### visualization
fc <- 1.5
fdr <- 0.05
object <- p3_m
s <- "Zilionis_et-p3"
marker <- p3_marker

module_info <- scale(module_score[colnames(object), ])
his <- c('LUAD', 'NSCLC')

sig_marker <-  marker[marker$FDR <= fdr & marker$Log2_FC >= fc, ]

###### select modules only DE in one cluster
dup_mods <-  marker %>% group_by(Gene) %>% summarise(numClu = length(unique(cluster))) %>% filter(numClu != 1)
target_modules <- sig_marker$Gene[!sig_marker$Gene %in% dup_mods$Gene]
names(target_modules) <- sig_marker[match(target_modules, sig_marker$Gene), 'cluster']


data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf, 
                   clusters = unique(object$cluster), his = his, data_sce = object,
                   module = target_modules)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]

col_scale <- c(-2, 0.1, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))


col_order <- c("Adeno (K11/K6/K15)",  "Mucinous (K59/K61)" )
row_order <- c("Adeno", "Mucinous")

row_split <- factor(names(target_modules), levels = row_order)

hm_expr <- module_expr
rownames <- mod_anno[match(rownames(hm_expr), mod_anno$Module), 'Label']
rownames(hm_expr) <- rownames$Label

cluster_col <- setNames(c("#FF7E00",  "#FCBE6E"),
                        c("Mucinous (K59/K61)",  "Adeno (K11/K6/K15)")
)

topAnno <- HeatmapAnnotation('Cluster' = colnames(hm_expr),
                             col = list(Cluster = cluster_col[colnames(hm_expr)]))

de_hm <- Heatmap(hm_expr, 
                 name = 'Relative Expression',
                 column_title = s,
                 row_title_rot  =90, 
                 col=col_fun, row_order = NULL, column_order = col_order,
                 column_split = NULL, row_split = row_split,
                 cluster_rows = F, cluster_columns = F, show_column_names = F,
                 row_gap = unit(2, "mm"),
                 row_names_gp = gpar(fontsize = 6),
                 row_km = NULL,
                 border = "black",
                 top_annotation = topAnno, #
                 #height = nrow(module_expr)*unit(3.5, "mm"),
                 height = unit(145, 'mm'), 
                 width = ncol(module_expr)*unit(12, "mm"),
                 heatmap_legend_param = list(
                   legend_direction = "horizontal", 
                   legend_width = unit(3, "cm")),
                 show_column_dend = F, show_row_dend = F)

pdf('../../Figure/Figure4/Zilionis_p3t_DE_heatmap.pdf')
draw(de_hm,
     heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()





