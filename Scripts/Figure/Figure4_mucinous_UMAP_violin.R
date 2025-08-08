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




### 2. Select samples and cell clusters that has high expression of club and globlet modules
target_modules <- c('L22', 'L19')
full_meta <- colData(data_sce)[, 1:14]
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



#### 2.1 subset cell clusters with high Club and Goblet modules
cf <- 0.005 ### based on MAST thres on cell probabilty, seems 0 is a reasonable cutoff
data <- module_summary(module_score = module_score, module_prob = module_prob, cf=cf, 
                       clusters = k, his = his, data_sce = data_sce,
                       module = target_modules)
module_per <- data[[3]]

club_k <- colnames(module_per)[module_per['L22', ] > 0.3]
mucinous_k <- colnames(module_per)[module_per['L19', ] > 0.3]
ks <- c(1, club_k, mucinous_k, 80) ### Cluster 1: oridinary LUAD cluster. Cluster 80: ordinary squamous cluster. 


#### 2.2 subset samples
sub_sce <- data_sce[, data_sce$celda_k %in% ks]
sub_sce$celda_k <- as.character(sub_sce$celda_k)
selected_sample <- names(table(sub_sce$sample))[table(sub_sce$sample) > 50] ### keep sample with > 50 cells

sub_sce <- data_sce[, data_sce$sample %in% selected_sample]
select_k <- names(table(sub_sce$celda_k))[table(sub_sce$celda_k) > 20] ### keep cluster with > 20 cells
sub_sce <- sub_sce[, sub_sce$celda_k %in% select_k]
sub_sce$celda_k <- as.character(sub_sce$celda_k)





### 3. plotting
histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]



####  3.1 only include P2
kept_cid <- colnames(sub_sce)[sub_sce$sample %in% c('Sinjab_et-P2-1-Epcam-pos')] 
p2_sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  
                                                                   kept_cid]),
                                 rowData = NULL,
                                 colData = colData(sub_sce)[kept_cid, ])
p2_sce_m <- getUMAP(p2_sce_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)

p2_sce_m$cluster <- p2_sce_m$celda_k
p2_sce_m$cluster[p2_sce_m$cluster %in% c('1', '16', '62', '48')] <- 'Club-like (K48/K62)'
p2_sce_m$cluster[p2_sce_m$cluster %in% c('59', '61', '63')] <- 'Mucinous (K61)'
p2_sce_m$cluster[p2_sce_m$cluster %in% c('64')] <- 'Intermediate (K64)'
p2_sce_m$cluster[p2_sce_m$cluster %in% c('3','4','15','22','56') & p2_sce_m$sample == 'Sinjab_et-P2-1-Epcam-pos'] <- 'Adeno-like (K15/56)'


##### 3.1.1 cluster label UMAP
l <- 'cluster'
cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#A65528", "#FCBE6E"),
                        c("Intermediate (K64)", "Mucinous (K61)", "Club-like (K48/K62)", "Adeno-like (K15/56)")
)

p2 <- plotSCEDimReduceColData(p2_sce_m, 
                              colorBy = l, reducedDimName = 'module_UMAP', 
                              clusterLabelSize = 4, dotSize = 0.3,
                              colorScale = cluster_col) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )

pdf('../../Figure/Figure4/Sinjab_P2_umap.pdf', width = 4, height = 4)
print(p2)
dev.off()


# ##### 3.1.2 violin plot
# ks <- c(1, 15, 48, 56, 61, 62, 64,80)
# l <- paste0('L', c(10, 15))  #6
# cid <- colnames(data_sce)[data_sce$celda_k %in% ks & data_sce$sample == 'Sinjab_et-P2-1-Epcam-pos']
# 
# d <- cbind(colData(data_sce)[cid, c('sample', 'celda_k')] %>% as.data.frame(), 
#            module_score[cid, l])
# d$celda_k <- as.character(d$celda_k)
# d$celda_k[d$celda_k %in% c(1)] <- 'LUAD'
# d$celda_k[d$celda_k %in% c(80)] <- 'LUSC'
# d$celda_k[d$celda_k %in% c(15, 56)] <- 'Adeno-like'
# d$celda_k[d$celda_k %in% c(48, 62)] <- 'Club-like'
# d$celda_k[d$celda_k %in% c(61)] <- 'Mucinous'
# d$celda_k[d$celda_k %in% c(64)] <- 'Intermediate'
# 
# d <- reshape2::melt(d, id = 'celda_k', measure.vars = l)
# colnames(d) <- c('Cluster', 'Module', 'Module score')
# 
# d$Module <- as.character(d$Module)
# d$Module[d$Module == 'L10'] = 'Adeno_1 (L10)'
# d$Module[d$Module == 'L15'] = 'Squamous_2 (L15)'
# 
# 
# d$Cluster <- factor(d$Cluster, levels = c('LUAD', 'Adeno-like', 'Club-like', 'Mucinous', 'Intermediate', 'LUSC')) 
# cluster_col <- setNames(c("#FAB4AD", "#A6CEE2","#1F78B4", "#FF7E00", "#A65528", "#FCBE6E"),
#                         c("LUAD", "LUSC", "Intermediate", "Mucinous", "Club-like", "Adeno-like")
# )
# 
# p2v <- ggplot(aes(x = Cluster, y = `Module score`, fill = Cluster), data = d) + 
#   geom_violin(position = position_dodge(0.8), width = 1) + facet_wrap(~Module, nrow=1) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 14), # vjust = -0.3, 
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 22),
#         legend.position="none",
#         strip.text.x = element_text(size = 22)) + 
#   scale_fill_manual(values = cluster_col)
# 
# 
# 
# ### combine figure
# p2f <- plot_grid(p2, p2v, ncol = 2)








#### 3.2 only include p030 
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Bischoff_et-p030t')] 
##### adding adeno-like cells as reference?
p3_sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                   kept_cid]),
                                 rowData = NULL,
                                 colData = colData(sub_sce)[kept_cid, ])
p3_sce_m <- getUMAP(p3_sce_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)

p3_sce_m$cluster <- p3_sce_m$celda_k
p3_sce_m$cluster[p3_sce_m$cluster %in% c('59', '61')] <- 'Mucinous (K59)'
p3_sce_m$cluster[p3_sce_m$cluster %in% c('63', '31')] <- 'Intermediate (K63)'




### cluster label UMAP
l <- 'cluster'
cluster_col <- setNames(c("#1F78B4", "#FF7E00"),
                        c("Intermediate (K63)", "Mucinous (K59)")
)

p3 <- plotSCEDimReduceColData(p3_sce_m, 
                              colorBy = l, reducedDimName = 'module_UMAP', 
                              clusterLabelSize = 4, dotSize = 0.3,
                              colorScale = cluster_col) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )

pdf('../../Figure/Figure4/Bischofft_p030t_umap.pdf', width = 4, height = 4)
print(p3)
dev.off()

# ### gene expression UMAP
# p3_sce_r <- SingleCellExperiment(assay = list(counts = assay(data_sce)[,
#                                                                        colnames(p3_sce_m)]),
#                                  rowData = NULL,
#                                  colData = colData(data_sce)[colnames(p3_sce_m), ])
# reducedDim(p3_sce_r, 'module_UMAP') <-  reducedDim(p3_sce_m, 'module_UMAP')[colnames(p3_sce_r), ]
# p3_sce_r <- runNormalization(p3_sce_r, normalizationMethod = 'LogNormalize', outAssayName = 'logNormalized')
# 
# 
# ### violin plot
# s <- 'Bischoff_et-p030t'
# ref_id <- colnames(data_sce)[data_sce$celda_k %in% c(1, 80)]
# ks <- c(59, 63)
# s_id <- colnames(data_sce)[data_sce$celda_k %in% ks & data_sce$sample %in% s]
# 
# l <- paste0('L', c(10, 15)) 
# cid <- c(ref_id, s_id)
# 
# d <- cbind(colData(data_sce)[cid, c('sample', 'celda_k')] %>% as.data.frame(), 
#            module_score[cid, l])
# d$celda_k <- as.character(d$celda_k)
# d$celda_k[d$celda_k %in% c(1)] <- 'LUAD'
# d$celda_k[d$celda_k %in% c(80)] <- 'LUSC'
# d$celda_k[d$celda_k %in% c(59)] <- 'Mucinous'
# d$celda_k[d$celda_k %in% c(63)] <- 'Intermediate'
# 
# d <- reshape2::melt(d, id = 'celda_k', measure.vars = l)
# colnames(d) <- c('Cluster', 'Module', 'Module score')
# 
# d$Module <- as.character(d$Module)
# d$Module[d$Module == 'L10'] = 'Adeno_1 (L10)'
# d$Module[d$Module == 'L15'] = 'Squamous_2 (L15)'
# 
# 
# d$Cluster <- factor(d$Cluster, levels = c('LUAD', 'Mucinous', 'Intermediate', 'LUSC')) 
# cluster_col <- setNames(c("#FAB4AD", "#A6CEE2","#1F78B4", "#FF7E00"),
#                         c("LUAD", "LUSC", "Intermediate", "Mucinous")
# )
# 
# p3v <- ggplot(aes(x = Cluster, y = `Module score`, fill = Cluster), data = d) + 
#   geom_violin(position = position_dodge(0.8), width = 1) + facet_wrap(~Module, nrow=1) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 14), # vjust = -0.3, 
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 22),
#         legend.position="none",
#         strip.text.x = element_text(size = 22)) + 
#   scale_fill_manual(values = cluster_col)


#### 3.3 only include p033t
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c('Bischoff_et-p033t')]
##### adding adeno-like cells as reference?
p3_sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                 kept_cid]),
                                 rowData = NULL,
                                 colData = colData(sub_sce)[kept_cid, ])
p3_sce_m <- getUMAP(p3_sce_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)

p3_sce_m$cluster <- p3_sce_m$celda_k
p3_sce_m$cluster[p3_sce_m$cluster %in% c('59', '61')] <- 'Mucinous (K59)'
p3_sce_m$cluster[p3_sce_m$cluster %in% c('63', '31')] <- 'Intermediate (K63)'


### cluster label UMAP
l <- 'cluster'
cluster_col <- setNames(c("#1F78B4", "#FF7E00"),
                        c("Intermediate (K63)", "Mucinous (K59)")
)

p3 <- plotSCEDimReduceColData(p3_sce_m, 
                              colorBy = l, reducedDimName = 'module_UMAP', 
                              clusterLabelSize = 4, dotSize = 0.3,
                              colorScale = cluster_col) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )

pdf('../../Figure/Figure4/Bischofft_p033t_umap.pdf', width = 4, height = 4)
print(p3)
dev.off()



### 3.4 for Zillionis_P3 P3
kept_cid <- colnames(sub_sce)[sub_sce$new_sample %in% c("Zilionis_et-p3t2", "Zilionis_et-p3t1", "Zilionis_et-p3t3")] 
##### adding adeno-like cells as reference?
p3_sce_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                 kept_cid]),
                                 rowData = NULL,
                                 colData = colData(sub_sce)[kept_cid, ])
p3_sce_m <- getUMAP(p3_sce_m, useAssay = 'module',
                    reducedDimName = 'module_UMAP',
                    sample = NULL, minDist = 0.2, seed = 2022,
                    useReducedDim = NULL, logNorm = F)

p3_sce_m$cluster <- p3_sce_m$celda_k
p3_sce_m$cluster[p3_sce_m$cluster %in% c('59', '61')] <- 'Mucinous (K59)'



### cluster label UMAP
l <- 'cluster'
cluster_col <- setNames(c("#1F78B4", "#FF7E00"),
                        c("Intermediate (K63)", "Mucinous (K59)")
)

p3 <- plotSCEDimReduceColData(p3_sce_m, 
                              colorBy = l, reducedDimName = 'module_UMAP', 
                              clusterLabelSize = 8, dotSize = 0.8,
                              colorScale = cluster_col) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    legend.position="none"
  )

pdf('../../Figure/Figure4/Zilionis_et_p3_umap.pdf', width = 4, height = 4)
print(p3)
dev.off()
