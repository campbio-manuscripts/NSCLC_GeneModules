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
library(here)

setwd(here::here("./Scripts/Figure"))


source('./Figure_utils.R')

### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs


### try to use all stage associated modules here. Define level of progression
stage_lm <- fread('../../Mixed_model_output/model_summary_Stage.csv')
stage_lm$Stage_fdr <- p.adjust(stage_lm$Stage_pval, 'fdr')

early_modules <- stage_lm$module[stage_lm$Stage_fdr <= 0.05 & stage_lm$Stage_coef < 0]
late_modules <- stage_lm$module[stage_lm$Stage_fdr <= 0.05 & stage_lm$Stage_coef > 0]
target_modules <- c(early_modules, late_modules)


#### select the LUAD cells
his <- 'LUAD'
full_meta <- colData(data_sce)
full_meta$stage_bin <- ifelse(full_meta$Stage %in% c('GGO', 'I'), 'Early', 'Late')
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells)


k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57)) ### perLUAD in 63 is 87.5, 12.5 NSCLC
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')
select_s <- unique(full_meta$sample[full_meta$celda_k %in% k & full_meta$Histology %in% his])
select_c <- rownames(full_meta)[full_meta$sample %in% select_s]



### 2. divide by percent of stage-related modules
early_pert <- rowSums(cell_probs[select_c, early_modules]) / rowSums(cell_probs[select_c, target_modules])
early_pert[is.na(early_pert)] <- 0 ### cells that have 0 expression of the early modules
late_pert <- rowSums(cell_probs[select_c, late_modules]) / rowSums(cell_probs[select_c, target_modules])


### 3. grouping cells based on the early_pert
pert_s <- data.frame('cb' = select_c, 
                     full_meta[select_c, c('sample', 'study', 'Histology', 'celda_k', 'Stage')],
                     #module_info[select_c, target_modules],
                     early_pert, late_pert) #early_pc1, late_pc1, all_pc1, early_pc1, 
pert_s <- pert_s %>% filter(sample %in% select_s & Histology %in% his)


pert_s$cs <- 'Early'
pert_s$cs[pert_s$early_pert <= 0.6] <- 'Mid'
pert_s$cs[pert_s$early_pert <= 0.4] <- 'Late' #0.2

pert_early <- pert_s %>% group_by(sample) %>% summarise(Early_enriched = mean(cs == 'Early'), 
                                                        Intermediate = mean(cs == 'Mid'),
                                                        Late_enriched = mean(cs == 'Late'),
                                                        earlyN = sum(cs == 'Early'),
                                                        midN = sum(cs == 'Mid'),
                                                        lateN = sum(cs == 'Late'),
                                                        'study' = unique(study),
                                                        'Histology' = unique(Histology),
                                                        'Stage' = unique(Stage))


### 4. label samples based on this classification

pert_early <- pert_early %>% arrange(desc(Early_enriched)) 
early_percent <- reshape2::melt(pert_early[, c('Stage', 'sample', 'Early_enriched', 'Intermediate', 'Late_enriched')])
early_percent$variable <- factor(early_percent$variable, levels = c('Late_enriched', 'Intermediate', 'Early_enriched')) #
colnames(early_percent)[colnames(early_percent) == 'variable'] <- 'Progression'
early_percent <- lapply(pert_early$sample, function(s) {
  early_percent[early_percent$sample == s, ] %>% arrange(desc(Progression)) ### better order for the barplot below. desc
}) %>% do.call('rbind', .)


##### create density plot / histogram for the progression score
additional.cutoffs <- c(0.6, 0.4) # additional bins

pdf("../../Figure/Figure4/Figure3_Progression_score_density.pdf", 
    width = 2, height = 1.5)
ggplot(pert_s, aes(x=early_pert)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill = 'white')+
  geom_density(alpha=.3) + xlab('Progression Score') + 
  geom_vline(xintercept=additional.cutoffs, colour="red3") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text = element_text (size = 8), 
        axis.title = element_text (size = 9))

dev.off()



### 5. plot the barplot
early_percent$plotOrder <- rep(1:length(unique(early_percent$sample)), each=3)
progression_barplot <- ggplot(aes(x = plotOrder, y = value, fill = Progression), data = early_percent) +  #factor(sample, levels = unique(early_percent$sample))
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_blank(),         # Change both x and y axis titles
    strip.text.x = element_text(size = 20),
    legend.key.size = unit(0.6, 'cm'),
    legend.title = element_text(size=14), #change legend title font size
    legend.text = element_text(size=12),
    #panel.spacing = unit(0.1, "lines")
    #legend.position = 'none',
  ) +     scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + 

  scale_fill_manual(breaks = c("Early_enriched", "Intermediate", "Late_enriched"),
                    values=c("#3C99B1", "#F1BB7B", "#C93428"))

fn <- file.path('../../Figure/Figure3', 'Figure3_Progression_barplot.pdf')
pdf(file = fn, width = 10, height = 6)
plot(progression_barplot)
dev.off()


##### quantify % of late-enriched cells in stage I or GGO samples
early_percent %>% group_by(Stage) %>% summarise(
  n = n()/3, 
  nL30 = sum(value[Progression == "Late_enriched"] > 0.3)
)



#### 5.1.2 Separate the barplot for different stages
tp <- "../../Table"
stage_barplot <- lapply(c("GGO", "I", "II", "III/IV"), function(st) { #unique(early_percent$Stage)
  sub_early_percent <- early_percent %>% filter(Stage == st)
  sub_early_percent$plotOrder <- rep(1:length(unique(sub_early_percent$sample)), each=3)
  progression_barplot <- ggplot(aes(x = plotOrder, y = value, fill = Progression), data = sub_early_percent) +  #factor(sample, levels = unique(early_percent$sample))
    geom_bar(position="fill", stat="identity") + theme_bw() + 
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title = element_blank(),         # Change both x and y axis titles
      strip.text.x = element_text(size = 20),
      legend.key.size = unit(0.3, 'cm'),
      legend.title = element_text(size=7), #change legend title font size
      legend.text = element_text(size=6),
      #panel.spacing = unit(0.1, "lines")
      legend.position = 'none',
    ) +     scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    
    scale_fill_manual(breaks = c("Early_enriched", "Intermediate", "Late_enriched"),
                      values=c("#3C99B1", "#F1BB7B", "#C93428")) + 
    ggtitle(st)
  
  progression_barplot
  
  if (st == "III/IV") {st = "IIIorIV"}
  write.csv(sub_early_percent %>% dplyr::rename(Percentage = 'value') %>% dplyr::select(-plotOrder), 
            file.path(tp, paste0("Early_Percentage_", st, "_Figure3D.csv")))
})






fn <- file.path('../../Figure/Figure4', 'Figure3_Progression_barplot_Stage.pdf')
pdf(file = fn, width = 10, height = 6)
cowplot::plot_grid(plotlist = stage_barplot, nrow = 2)
dev.off()




### 6. UMAP  to visualize some heterogeneous sample
ip <- '../../Mixed_model_output/'
luad_lm <- fread(file.path(ip, 'model_summary_Stage.csv'))
skip_m <- luad_lm$module[luad_lm$study_VarPer > 0.3]

cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')

luad_res <- SingleCellExperiment(assay = list(corrected_module = cr_module[target_modules, select_c]), #!rownames(cr_module) %in% skip_m
                                 rowData = NULL,
                                 colData = colData(data_sce)[select_c,])
luad_res$cellDeType <- pert_s[colnames(luad_res), 'cs']
luad_res$early_pert <- early_pert
luad_res$late_pert <- late_pert



##### prepare data for visualization
# hs <- c('Xing_et-SSN22',
#         'Sinjab_et-P1-1', 
#         'Wu_et-P2')
#s <- hs[1]

s <- 'Xing_et-SSN22'
sce <- luad_res[, luad_res$sample %in% s]
sce <- singleCellTK::getUMAP(sce, useAssay = 'corrected_module',
                             reducedDimName = 'module_UMAP', useReducedDim = NULL, 
                             sample = NULL, minDist = 0.2, seed = 2022, logNorm = F, 
                             pca = T)

colData(sce) <- colData(sce) %>% as.data.frame() %>% mutate(
  CellStage = case_when(
    
    cellDeType == 'Early' ~ 'Early_enriched',
    cellDeType == 'Mid' ~ 'Intermediate',
    cellDeType == 'Late' ~ 'Late_enriched',
  ),
  
  `EarlyStage` = early_pert, 
  
  `LateStage` = late_pert
) %>% DataFrame(.)

d <- colData(sce)[, c('sample','CellStage', 'EarlyStage', 'LateStage'), drop=F] %>% as.data.frame()
d <- reshape2::melt(d, id = 'CellStage', measure.vars = c('EarlyStage', 'LateStage'))
colnames(d) <- c('CellStage', 'Module', '% counts')
d$Module <- factor(d$Module, levels = c('EarlyStage', 'LateStage'))
d$CellStage  <- factor(d$CellStage , levels = c('Early_enriched', 'Intermediate', 'Late_enriched')) 




##### plotting

cluster_col <- setNames(c("#3C99B1", "#F1BB7B", "#C93428"),
                        c('Early_enriched', 'Intermediate', 'Late_enriched')
)

# p2v <- ggplot(aes(x = CellStage, y = `% counts`, fill = CellStage), data = d) + 
#   geom_boxplot(position = position_dodge(0.8), width = 1) + facet_wrap(~Module, nrow=1) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 14), # vjust = -0.3, 
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 22),
#         legend.position="none",
#         strip.text.x = element_text(size = 22),) + 
#   scale_fill_manual(values = cluster_col)


#### plotting the cluster label
p2u_label <- plotSCEDimReduceColData(sce, 
                               colorBy = 'CellStage', reducedDimName = 'module_UMAP',
                               title = paste(s, unique(sce$Stage), sep='-Stage'), #
                               titleSize = 10,
                               colorScale = cluster_col,
                               axisSize = 0,
                               axisLabelSize = F,
                               dotSize = 1, 
                               baseSize = 10,
                               legendTitleSize = 10,
                               clusterLabelSize = 8,
                               labelClusters = F) +  theme_bw() +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position= c(0.2, 0.3),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'),
        legend.text = element_text(size=13)
  ) + labs(y= "UMAP_2", x = "UMAP_1") + 
  guides(colour = guide_legend(override.aes = list(size=4)))


fn <- file.path('../../Figure/Figure4', 'Figure3_Progression_UMAP.pdf')
pdf(file = fn, width = 7, height = 7)
plot(p2u_label)
dev.off()


#### plotting the progression score / % early module
# col_scale <- c(0, 0.3, 0.6, 0.8, 1)
# colors <- c('white', 'snow','orange','red', 'red3')
# col_scale <- c(0, 0.3, 0.7, 1)
# colors <- c('#0042AD','#006DCA','#F96B4A','red3')
col_scale <- c(0, 0.4, 0.8, 1)
colors <- c('#0042AD','orange','red','red3')

p2u_score <- plotSCEDimReduceColData(sce, 
                                     colorBy = 'LateStage', reducedDimName = 'module_UMAP',
                                     title = paste(s, unique(sce$Stage), sep='-Stage'), #
                                     titleSize = 10,
                                     axisSize = 0,
                                     axisLabelSize = F,
                                     dotSize = 1, 
                                     baseSize = 10,
                                     legendTitle = 'Progression\nScore',
                                     legendTitleSize = 8,
                                     legendSize = 6,
                                     labelClusters = F) +  
  scale_color_gradientn(values = col_scale, 
                        colours = colors,
                        na.value = "grey50") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position= c(0.2, 0.3),
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text = element_text(size=13)
  ) + labs(y= "UMAP_2", x = "UMAP_1") 





#### further check the 2 late-enriched clusters in sample Xing_et-SSN22
sce$late_cluster <- ifelse(sce$celda_k %in% c('22', '64'), 'late', 'other')


#### redefine cluster within this sample
sce <- runSeuratPCA(inSCE = sce, useAssay = "corrected_module", reducedDimName = "module_pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = rownames(sce))
sce <- runSeuratFindClusters(inSCE = sce, useAssay = "corrected_module", useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
sce$seurat_cluster <- sce$Seurat_louvain_Resolution0.8
# sce <- runSeuratUMAP(sce, useReduction = c("pca"))
# plotSeuratReduction(sce, useReduction = "umap", showLegend = TRUE)
# 
# 
# plotSCEDimReduceColData(sce,
#                         colorBy = 'seurat_cluster', reducedDimName = 'module_UMAP',
#                         title = paste(s, unique(sce$Stage), sep='-Stage'), #
#                         titleSize = 10,
#                         colorScale = NULL,
#                         axisSize = 0,
#                         axisLabelSize = F,
#                         dotSize = 1,
#                         legendSize = 6,
#                         legendTitleSize = 5,
#                         clusterLabelSize = 8,
#                         labelClusters = T) +  theme_bw()

late_clu <- c('2', '6')
late_de <- list()

for (clu in late_clu) {
  target <- clu
  control <- as.character(unique(sce$Seurat_louvain_Resolution0.8));control <- control[!control %in% target]
  sce_sub <- sce[, sce$seurat_cluster %in% c(target, control)]
  
  sce_sub$seurat_cluster <- as.character(sce_sub$Seurat_louvain_Resolution0.8)
  sce_sub <- runDEAnalysis(inSCE = sce_sub, method = "wilcox", useAssay = "corrected_module",
                           class = "seurat_cluster", classGroup1 = target, classGroup2 = control,
                           groupName1 = target, groupName2 = "Low-Enriched", 
                           analysisName = target, overwrite = T)
  plotDEGHeatmap(sce_sub, useResult = target, log2fcThreshold = 0.3, rowLabel = TRUE)
  
  late_de[[clu]] <- metadata(sce_sub)$diffExp[[1]]$result
}

clu2_de <- late_de[[1]] %>% filter(FDR < 0.01, Log2_FC > 0.6)
clu6_de <- late_de[[2]] %>% filter(FDR < 0.01, Log2_FC > 0.6)



### plotting some distinct module in these two clusters

setdiff(clu2_de$Gene, clu6_de$Gene)
#l <- 'L32'

### high in cluster 2
l <- 'L47'
t <- 'L47 Chemokine_1: CXCl1,CXCL3'

### high in cluster 6

setdiff(clu6_de$Gene, clu2_de$Gene)
#l <- 'L61'
l <- 'L59'
t <- 'L59 ECM_4: COL1A1,MMP9'

col_scale <- c(0, 0.3, 0.6, 0.8, 1)
colors <- c('white', 'snow','orange','red', 'red3')

plotSCEDimReduceFeatures(sce, feature=l, reducedDimName = 'module_UMAP', useAssay = 'corrected_module',
                         legendTitle = 'Corrected\nScore',
                         dotSize = 0.5, title = t, titleSize = 16, legendSize = 10) + 
  scale_color_gradientn(values = col_scale, 
                        colours = colors,
                        na.value = "grey50") + 
  theme(
    #axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(5, 'mm'), 
    legend.margin=margin(0,0,0,0)
  )




#### use heatmap to visualize the DE modules

##### 1. get de modules
all_de <- union(clu2_de$Gene, clu6_de$Gene)



#hm_expr <- as.matrix( assay(sce, 'corrected_module')[all_de, ])
#hm_expr <- t(scale(t(hm_expr)))

col_scale <- c(-2, 0, 2) #, 2.6
col_fun = circlize::colorRamp2(col_scale,  c('royalblue', 'white', 'red2'))


##### 2. define column split
sce$Seurat_late_clu <- sce$Seurat_louvain_Resolution0.8 %>% as.character
sce$Seurat_late_clu[!sce$Seurat_late_clu %in% c('2', '6')] <- 'Early_enriched'
sce$Seurat_late_clu[sce$Seurat_late_clu %in% c('2')] <- 'Late_enriched_1'
sce$Seurat_late_clu[sce$Seurat_late_clu %in% c('6')] <- 'Late_enriched_2'

# sce$Progress_Clu_fine <- sce$CellStage
# sce$Progress_Clu_fine[sce$seurat_cluster == '2'] <- 'Late_enriched_1' 
# sce$Progress_Clu_fine[sce$seurat_cluster == '6'] <- 'Late_enriched_2' 


##### 3. summarize the module expression in three different clusters
cf <- 0.005
object <- sce
his <- 'LUAD'
clus <- 'Seurat_late_clu'
target_modules <- all_de
names(target_modules) <- ifelse(target_modules %in% clu2_de$Gene, 'Up_LateEnriched_1', 'Up_LateEnriched_2')
names(target_modules)[target_modules %in% clu2_de$Gene & target_modules %in% clu6_de$Gene] <- 'Up_Both_LateEnriched'


module_info <- scale(t(cr_module)[colnames(object), target_modules])

data <- hm_summary(module_score = module_info, module_prob = module_prob, cf=cf, 
                   clusters = unique(object[[clus]]), his = his, data_sce = object,
                   module = target_modules, clu_col = clus)
sample_module <- data[[1]]; module_expr <- data[[2]]; module_per <- data[[3]]
hm_expr <- module_expr



column_split <- colnames(hm_expr)
cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#A65528"),
                        c("Early_enriched", "Late_enriched_1", "Late_enriched_2")
)

topAnno <- HeatmapAnnotation('Prgression' = column_split,
                             #'Cluster' = sce$Progress_Clu_fine,
                             col = list(Prgression = cluster_col) #, Cluster = cluster_col2)
                                        )
row_split <- names(target_modules)


### adding the module annotation to the rowname
module_anno <- fread('../../Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))

rownames(hm_expr) <- module_anno[match(rownames(hm_expr), module_anno$Module), ] %>% pull(Label)

late_hm <- Heatmap(hm_expr, 
        name = "Relative Expression",
        #column_title = s,
        row_title_rot  = 90, 
        col=col_fun, 
        row_order = NULL, column_order = NULL,
        column_split = column_split, row_split = row_split,
        cluster_rows = T, cluster_columns = F, show_column_names = F,
        # row_gap = unit(2, "mm"),
        # row_names_gp = gpar(fontsize = 8),
        # row_km = NULL,
        border = "black",
        top_annotation = topAnno, #
        show_column_dend = F, show_row_dend = F, 
        column_title_rot = 0)

fn <- '../../Figure/Figure4/Figure3_Progression_heatmap.pdf'
pdf(file = fn, width = 11, height = 7)
print(late_hm)
dev.off()





#### plotting the matched annotation umap
sce$Progress_Clu_fine <- sce$Seurat_late_clu
#sce$Progress_Clu_fine <- sce$CellStage

#sce$Progress_Clu_fine[sce$Seurat_late_clu == 'Early_enriched' & sce$CellStage == 'Intermediate'] <- 'Intermediate'
#sce$Progress_Clu_fine[sce$seurat_cluster == '4'] <- 'Intermediate' ### only 18 late enriched cells in cluster 4, which looks like intermediate cluster

# sce$Progress_Clu_fine[sce$seurat_cluster == '2'] <- 'Late_enriched_1' 
# sce$Progress_Clu_fine[sce$seurat_cluster == '6'] <- 'Late_enriched_2' 


# cluster_col <- setNames(c("#3C99B1", "#F1BB7B", 'red',"#C93428", 'tan4'),
#                         c('Early_enriched', 'Intermediate', 'Late_enriched' ,'Late_enriched_1', 'Late_enriched_2')
# )

cluster_col <- setNames(c("#1F78B4", "#FF7E00", "#A65528"),
                        c("Early_enriched", "Late_enriched_1", "Late_enriched_2")
)


p2ProU <- plotSCEDimReduceColData(sce, 
                          colorBy = 'Progress_Clu_fine', reducedDimName = 'module_UMAP',
                          title = paste(s, unique(sce$Stage), sep='-Stage'), #
                          titleSize = 10,
                          colorScale = cluster_col,
                          axisSize = 0,
                          axisLabelSize = F,
                          dotSize = 1, 
                          baseSize = 10,
                          legendTitleSize = 10,
                          clusterLabelSize = 8,
                          labelClusters = F) +  theme_bw() +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position= c(0.2, 0.3),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'),
        legend.text = element_text(size=13)
  ) + labs(y= "UMAP_2", x = "UMAP_1") + 
  guides(colour = guide_legend(override.aes = list(size=4)))

p2ProU


