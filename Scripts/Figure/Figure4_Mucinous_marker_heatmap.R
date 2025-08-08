library(stringr)
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


#### check a few Mucinous marker genes in Sinjab_P2 sample
p <- '/restricted/projectnb/camplab/home/rzhong/projects/luad_preprocess/Kadara_et/'
sinjab_sce <- readRDS(file.path(p, 'luad_Cells_SCE.rds'))


#### luckly, we got all these marker genes in this study! 
mad_marker <- c('MUC5AC', 'KRT7', "PRSS1", "CDX2", "MUC2", "PDX1", "TFF2", "MUC6", "REG4", "TFF1", "HNF4A") #"F2", "CPS1", "FGA", 
#mad_marker %in% rownames(sinjab_sce) 
scc_marker <- c('TP63', 'KRT5',  'KRT14', 'DSG3', 'KRT6A', 'KRT6B', 'PAX9', 'SOX2') #'KRT1',
lung_marker <- c('NKX2-1', 'SFTA3', 'LMO3', 'NAPSA', 'SFTPC', 'FOXA2', 'HNF1B')

target_marker <- c(scc_marker, 
                   lung_marker, 
                   mad_marker)
names(target_marker) <- c(rep('SCC', length(scc_marker)), 
                          rep('Lung', length(lung_marker)), 
                          rep('MAD', length(mad_marker)))


#### get the Sinjab_p2 index from our meta-analysis
kept_cid <- colnames(p2_m)
sinjab_p2_sce <- sinjab_sce[, kept_cid]

col_order <- c("Adeno-like (K15/56)",  "Intermediate (K64)", "Club-like (K48/K62)", "Mucinous (K61)")
sinjab_p2_sce$cluster <- factor(p2_m$cluster, levels = col_order)


#### plotting P2 Marker heatmap
# target_marker <- target_marker[!target_marker %in% c('F2', 'KRT1', 'FGA')]

object <- sinjab_p2_sce
assay <- 'logcounts'
column_split <- sinjab_p2_sce$cluster
s <- 'Sinjab_p2'
cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#A65528", "#FCBE6E"),
                        c("Intermediate (K64)", "Mucinous (K61)", "Club-like (K48/K62)", "Adeno-like (K15/56)")
)

hm_expr <- as.matrix( assay(object, assay)[target_marker, ])
hm_expr <- t(scale(t(hm_expr)))

col_scale <- c(-2, 0, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))

#row_order <- c("Adeno-like", "Intermediate", "Club", "Mucinous")



row_split <- factor(names(target_marker), levels = c('Lung', 'SCC', 'MAD'))
topAnno <- HeatmapAnnotation('Cluster' = column_split,
                             col = list(Cluster = cluster_col))


p2_marker_hm <- Heatmap(hm_expr, 
        name = 'Relative\nExpression',
        column_title = s,
        row_title_rot  = 90, 
        col=col_fun, 
        row_order = NULL, column_order = NULL,
        column_split = column_split, row_split = row_split,
        cluster_rows = F, cluster_columns = F, show_column_names = F,
        # row_gap = unit(2, "mm"),
        # row_names_gp = gpar(fontsize = 8),
        # row_km = NULL,
        border = "black",
        top_annotation = topAnno, #
        show_column_dend = F, show_row_dend = F)





#### look at p30t in Bischoff_et
fp <- '/restricted/projectnb/camplab/home/rzhong/projects/luad_preprocess/s41388-021-02054-3'
bischoff_sce <- readRDS(file.path(fp, 'rawCount_filter_SCE.rds'))

kept_cid <- colnames(p30_m)
bischoff_p30t_sce <- bischoff_sce[, kept_cid]


col_order <- c("Adeno (K11)",  "Intermediate (K63)",  "Mucinous (K59)")
bischoff_p30t_sce$cluster <- factor(p30_m$cluster, levels = col_order)



object <- bischoff_p30t_sce
assay <- 'logNormalized'
column_split <- bischoff_p30t_sce$cluster
s <- 'Bischoff_p30t' #Bischoff_p30t
cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#FCBE6E"),
                        c("Intermediate (K63)", "Mucinous (K59)", "Adeno (K11)")
)



hm_expr <- as.matrix( assay(object, assay)[target_marker[target_marker %in% rownames(object)], ])
hm_expr <- t(scale(t(hm_expr)))

col_scale <- c(-2, 0, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))
row_split <- factor(names(target_marker), levels = c('Lung', 'SCC', 'MAD'))
topAnno <- HeatmapAnnotation('Cluster' = column_split,
                             col = list(Cluster = cluster_col))



p30_marker_hm <- Heatmap(hm_expr, 
        name = 'Relative\nExpression',
        column_title = s,
        row_title_rot  = 90, 
        col=col_fun, 
        row_order = NULL, column_order = NULL,
        column_split = column_split, row_split = row_split,
        cluster_rows = F, cluster_columns = F, show_column_names = F,
        # row_gap = unit(2, "mm"),
        # row_names_gp = gpar(fontsize = 8),
        # row_km = NULL,
        border = "black",
        top_annotation = topAnno, #
        show_column_dend = F, show_row_dend = F)




#### look at p33t in Bischoff_et
fp <- '/restricted/projectnb/camplab/home/rzhong/projects/luad_preprocess/s41388-021-02054-3'
bischoff_sce <- readRDS(file.path(fp, 'rawCount_filter_SCE.rds'))

kept_cid <- colnames(p30_m)
bischoff_p30t_sce <- bischoff_sce[, kept_cid]


col_order <- c("Adeno (K11)",  "Intermediate (K63)",  "Mucinous (K59)")
bischoff_p30t_sce$cluster <- factor(p30_m$cluster, levels = col_order)



object <- bischoff_p30t_sce
assay <- 'logNormalized'
column_split <- bischoff_p30t_sce$cluster
s <- 'Bischoff_p33t' #Bischoff_p30t
cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#FCBE6E"),
                        c("Intermediate (K63)", "Mucinous (K59)", "Adeno (K11)")
)



hm_expr <- as.matrix( assay(object, assay)[target_marker[target_marker %in% rownames(object)], ])
hm_expr <- t(scale(t(hm_expr)))

col_scale <- c(-2, 0, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))
row_split <- factor(names(target_marker), levels = c('Lung', 'SCC', 'MAD'))
topAnno <- HeatmapAnnotation('Cluster' = column_split,
                             col = list(Cluster = cluster_col))



p33_marker_hm <- Heatmap(hm_expr, 
                         name = 'Relative\nExpression',
                         column_title = s,
                         row_title_rot  = 90, 
                         col=col_fun, 
                         row_order = NULL, column_order = NULL,
                         column_split = column_split, row_split = row_split,
                         cluster_rows = F, cluster_columns = F, show_column_names = F,
                         # row_gap = unit(2, "mm"),
                         # row_names_gp = gpar(fontsize = 8),
                         # row_km = NULL,
                         border = "black",
                         top_annotation = topAnno, #
                         show_column_dend = F, show_row_dend = F)



#### look at "Zilionis_et-p3"
fp <- '/restricted/projectnb/camplab/home/rzhong/projects/luad_preprocess/Salcher_Altas/Zilionis_Klein_2019'
Zilionis_sce <- readRDS(file.path(fp, 'allCells_SCE.rds'))

kept_cid <- colnames(p3_m)
Zilionis_p3t_sce <- Zilionis_sce[, kept_cid]
Zilionis_p3t_sce <- runNormalization(Zilionis_p3t_sce, useAssay = "counts", outAssayName = "logNormalized", 
                                     normalizationMethod = 'LogNormalize')

col_order <- c("Adeno (K11/K6/K15)", "Mucinous (K59/K61)")
Zilionis_p3t_sce$cluster <- factor(p3_m$cluster, levels = col_order)



object <- Zilionis_p3t_sce
assay <- 'logNormalized'
column_split <- Zilionis_p3t_sce$cluster
s <- 'Zilionis_p3'
cluster_col <- setNames(c("#FF7E00", "#FCBE6E"),
                        c("Mucinous (K59/K61)", "Adeno (K11/K6/K15)")
)



hm_expr <- as.matrix( assay(object, assay)[target_marker[target_marker %in% rownames(object)], ])
hm_expr <- t(scale(t(hm_expr)))

col_scale <- c(-2, 0, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))
row_split <- factor(names(target_marker[target_marker %in% rownames(object)]), levels = c('Lung', 'SCC', 'MAD'))
topAnno <- HeatmapAnnotation('Cluster' = column_split,
                             col = list(Cluster = cluster_col))


p3_marker_hm <- Heatmap(hm_expr, 
        name = 'Relative\nExpression',
        column_title = s,
        row_title_rot  = 90, 
        col=col_fun, 
        row_order = NULL, column_order = NULL,
        column_split = column_split, row_split = row_split,
        cluster_rows = F, cluster_columns = F, show_column_names = F,
        # row_gap = unit(2, "mm"),
        # row_names_gp = gpar(fontsize = 8),
        # row_km = NULL,
        border = "black",
        top_annotation = topAnno, #
        show_column_dend = F, show_row_dend = F)



pdf("../../Figure/Supplementary/Supplementary_Figure10.pdf", width = 7, height = 5)
print(p2_marker_hm)
print(p30_marker_hm)
print(p33_marker_hm)
print(p3_marker_hm)
dev.off()



#### 2.2 Try to look at Sinjab_P3, which is annotated as mucinous
p3_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                             colnames(data_sce)[data_sce$new_sample == 'Sinjab_et-P3-1-Epcam-pos']]),
                             rowData = NULL,
                             colData = colData(data_sce)[data_sce$new_sample == 'Sinjab_et-P3-1-Epcam-pos', ])

p3_m$cluster <- as.character(p3_m$celda_k)
p3_m$cluster[!p3_m$cluster %in% c(17, 23, 7)] <- 'others'
p3_m <- findMarkerDiffExp(p3_m,
                          useAssay = "module",
                          cluster = "cluster",
                          method = c('wilcox'),
                          log2fcThreshold = 0.6)

p3_m <- getUMAP(p3_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)

plotSCEDimReduceColData(p3_m, 
                        colorBy = 'cluster', reducedDimName = 'module_UMAP', 
                        clusterLabelSize = 8, dotSize = 0.8,
                        colorScale = NULL) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )


p3_marker <- metadata(p3_m)$findMarker

### unfortunately, no significant higher mucinous expression in any of the clusters. 
plotSCEViolinAssayData(p3_m, feature = 'L19', useAssay = 'module', groupBy = 'cluster')

plotSCEViolinAssayData(p2_m, feature = 'L19', useAssay = 'module', groupBy = 'cluster')
