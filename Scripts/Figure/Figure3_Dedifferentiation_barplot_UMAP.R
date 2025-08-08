
### late stage moduleL

#### L45:COL1A1; PXDN; MMP9; COL13A1; ADAMTS14; COL23A1
#### L46:ITGB1; LAMB3;  ITGB4 
#### L64:HMGA2


### Early stage module:
# L1	Surfactant_C	SFTPC
# L2	Surfactant_A	SFTPA1, SFTPA2
# L4	Surfactant_BD	SFTPB, SFTPD
# L5	Adeno_1	NAPSA, RNASE1
# L6	Adeno_2	NKX2-1
#L19	Secretory_9	AGR2
#L21	Secretory_11	MUC1, JUND
####

library(celda)
library(singleCellTK)
library(data.table)
library(dplyr)
library(ggplot2)

setwd(here::here("./Scripts/Figure"))


### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs


### Define early-stage and late-stage associated modules
### early modules: SFTPC, SFTPA1/2, SFTPB/D, NAPSA, NKX2-1, SCGB3A2, SCGB3A1
### late modules: COL1A1/MMP9, ITGB1/LAMB3/ITGB4, HMGA2

early_modules <- paste0('L', c(1:4, 10, 12, 26, 28))
late_modules <- paste0('L', c(59, 60, 64)) # 26,
target_modules <- c(early_modules, late_modules)


### try to use all stage associated modules here
stage_lm <- fread('../../Mixed_model_output/model_summary_Stage.csv')
stage_lm$Stage_fdr <- p.adjust(stage_lm$Stage_pval, 'fdr')

early_modules <- stage_lm$module[stage_lm$Stage_fdr <= 0.05 & stage_lm$Stage_coef < 0]
#early_modules <- early_modules[!early_modules %in% c('L51', 'L82', 'L50')]

late_modules <- stage_lm$module[stage_lm$Stage_fdr <= 0.05 & stage_lm$Stage_coef > 0]
#late_modules <- late_modules[!late_modules %in% c('L22')]

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

# pert_s$cs <- 'Early'
# pert_s$cs[pert_s$early_pert <= 0.8] <- 'Mid'
# pert_s$cs[pert_s$early_pert <= 0.2] <- 'Late' #0.2
pert_s$cs <- 'Early'
pert_s$cs[pert_s$early_pert <= 0.6] <- 'Mid'
pert_s$cs[pert_s$early_pert <= 0.4] <- 'Late' #0.2

pert_early <- pert_s %>% group_by(sample) %>% summarise(Differentiated = mean(cs == 'Early'), 
                                                        Intermediate = mean(cs == 'Mid'),
                                                        `De-differentiated` = mean(cs == 'Late'),
                                                        earlyN = sum(cs == 'Early'),
                                                        midN = sum(cs == 'Mid'),
                                                        lateN = sum(cs == 'Late'),
                                                        'study' = unique(study),
                                                        'Histology' = unique(Histology),
                                                        'Stage' = unique(Stage))



# ##### create density plot / histogram for the progression score
# additional.cutoffs <- c(0.6, 0.4) # additional bins
# 
# pdf("../../Figure/Figure4/Figure3_Progression_score_density.pdf")
# ggplot(pert_s, aes(x=early_pert)) + 
#   #geom_histogram(aes(y=..density..), colour="black", fill = 'white')+
#   geom_density(alpha=.3) + xlab('Progression Score') + 
#   geom_vline(xintercept=additional.cutoffs, colour="red3") + 
#   theme_classic()
# dev.off()




### 4. label samples based on this classification

pert_early <- pert_early %>% arrange((Differentiated)) 
early_percent <- reshape2::melt(pert_early[, c('Stage', 'sample', 'Differentiated', 'Intermediate', 'De-differentiated')])
early_percent$variable <- factor(early_percent$variable, levels = c('De-differentiated', 'Intermediate', 'Differentiated')) #
colnames(early_percent)[colnames(early_percent) == 'variable'] <- 'Level of differentiation'
early_percent <- lapply(pert_early$sample, function(s) {
  early_percent[early_percent$sample == s, ] %>% arrange(desc(`Level of differentiation`)) ### better order for the barplot below
}) %>% do.call('rbind', .)


### 5. plot the barplot
early_percent$plotOrder <- rep(1:length(unique(early_percent$sample)), each=3)
p2 <- ggplot(aes(x = plotOrder, y = value, fill = `Level of differentiation`), data = early_percent) +  #factor(sample, levels = unique(early_percent$sample))
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_blank(),         # Change both x and y axis titles
    strip.text.x = element_text(size = 20),
    legend.key.size = unit(1.2, 'cm'),
    legend.title = element_text(size=14), #change legend title font size
    legend.text = element_text(size=12)
    #panel.spacing = unit(0.1, "lines")
    #legend.position = 'none',
  ) + 
  scale_fill_manual(breaks = c("Differentiated", "Intermediate", "De-differentiated"),
                    values=c("#3C99B1", "#F1BB7B", "#C93428"))

plot(p2)




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
hs <- c('Xing_et-SSN22',
        'Sinjab_et-P1-1', 
        'Wu_et-P2')
s <- hs[1]
sce <- luad_res[, luad_res$sample %in% s]
sce <- singleCellTK::getUMAP(sce, useAssay = 'corrected_module',
               reducedDimName = 'module_UMAP', useReducedDim = NULL, 
               sample = NULL, minDist = 0.2, seed = 2022, logNorm = F, 
               pca = T)

colData(sce) <- colData(sce) %>% as.data.frame() %>% mutate(
  CellStage = case_when(
    
      cellDeType == 'Early' ~ 'Differentiated',
      cellDeType == 'Mid' ~ 'Intermediate',
      cellDeType == 'Late' ~ 'Dedifferentiated',
  ),
  
  `EarlyStage` = early_pert, 
  
  `LateStage` = late_pert
) %>% DataFrame(.)

d <- colData(sce)[, c('sample','CellStage', 'EarlyStage', 'LateStage'), drop=F] %>% as.data.frame()
d <- reshape2::melt(d, id = 'CellStage', measure.vars = c('EarlyStage', 'LateStage'))
colnames(d) <- c('CellStage', 'Module', '% counts')
d$Module <- factor(d$Module, levels = c('EarlyStage', 'LateStage'))
d$CellStage  <- factor(d$CellStage , levels = c('Differentiated', 'Intermediate', 'Dedifferentiated')) 




##### plotting

cluster_col <- setNames(c("#3C99B1", "#F1BB7B", "#C93428"),
                        c('Differentiated', 'Intermediate', 'Dedifferentiated')
)

p2v <- ggplot(aes(x = CellStage, y = `% counts`, fill = CellStage), data = d) + 
  geom_boxplot(position = position_dodge(0.8), width = 1) + facet_wrap(~Module, nrow=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 14), # vjust = -0.3, 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22),
        legend.position="none",
        strip.text.x = element_text(size = 22),) + 
  scale_fill_manual(values = cluster_col)

p2u <- plotSCEDimReduceColData(sce, 
                               colorBy = 'CellStage', reducedDimName = 'module_UMAP',
                               title = paste(s, unique(sce$Stage), sep='-Stage'), #
                               titleSize = 10,
                               colorScale = cluster_col,
                               axisSize = 0,
                               axisLabelSize = F,
                               dotSize = 1, 
                               legendSize = 6, 
                               legendTitleSize = 5,
                               clusterLabelSize = 8,
                               labelClusters = F) +  theme_bw() +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        #legend.position="none"
        )


plot(
  cowplot::plot_grid(p2u, p2v, ncol=2,
                     labels = NULL)
)

