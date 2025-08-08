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

#### module probability
x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs

cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')

### 2. load linear model summary result
histo_m_summary <- read.csv('../..//Mixed_model_output/model_summary_Histology.csv')

module_anno <- fread('../..//Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))



### 3. Generate % of module express in each cluster
full_meta <- colData(data_sce)
cell_prob_cf <- 0.01
sample_prob_cf <- 0.1
cluster_prob_cf <- 0.05
  
  
mod_sample_prob <- lapply(colnames(module_prob), function(l) {
  ls <- lapply(unique(full_meta$celda_k), function(s) { #sample
    #sid <- rownames(full_meta)[full_meta$sample == s]
    sid <- rownames(full_meta)[full_meta$celda_k == s]
    mean(module_prob[sid, l] > cell_prob_cf)
  }) %>% unlist()
  ifelse(ls > cluster_prob_cf, 1, 0) #sample_prob_cf
})

mod_sample_p <- do.call('rbind', mod_sample_prob)
rownames(mod_sample_p) <- rownames(cr_module); colnames(mod_sample_p) <- unique(full_meta$celda_k) #

mod_sample_p <- rowMeans(mod_sample_p == 1)
#module_diversity <- apply(mod_sample_p, 1, function(x) {diversity(x,"shannon")})


### 3. compute the Gini Index for modules
y <- factorizeMatrix(data_sce, type = "counts")
cell_m_counts <- y[[1]][['cell']]; cell_m_counts <- t(cell_m_counts)
colSum_sizeFactor <- colSums(t(cell_m_counts)) / 10000 
normalized_m_counts <- sweep(cell_m_counts, 1, colSum_sizeFactor, "/")
module_activity <- log2(normalized_m_counts + 1)

library(DescTools)
Gini(module_activity[, 2])

mod_Gini <- lapply(colnames(module_activity), function(l) {
  Gini(module_activity[, l])
})
mod_Gini <- unlist(mod_Gini); names(mod_Gini) <- colnames(module_activity)

### 4. summarize all data together
module_d <- data.frame('module' = histo_m_summary$module, 
                       'Category' = module_anno[match(histo_m_summary$module, module_anno$Module), 'Category'],
                       'Per_Study' = histo_m_summary$study_VarPer, 
                       'Per_Histology+Stage' = histo_m_summary$Fixed_VarPer, 
                       'Per_Sample' = histo_m_summary$sample_VarPer, 
                       'Per_Cluster' = histo_m_summary$celda_VarPer, 
                       'Per_Unexplained' = histo_m_summary$Residual_VarPer, 
                       'Per_expressed' = mod_sample_p[histo_m_summary$module],
                       'Gini_coefficient' = mod_Gini[histo_m_summary$module])
                       #'Diversity' = module_diversity[ histo_m_summary$module])

module_dl <- reshape2::melt(module_d)
module_dl <- module_dl %>% filter(Category != 'Unannotated')
module_dl$Category <- factor(module_dl$Category, levels = c( 
                                                            'Lineage', 'Xenobiotic metabolism', 
                                                            'Proliferation', 'Cellular process', 
                                                            'Antigen Presentation', 'Signalling', 
                                                            'Extracellular matrix organization', 
                                                            'Stress response', 
                                                            'Housekeeping') %>% rev()) #, 'Unannotated'



xcf <- 0.75
study_box <- ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_Study')) + 
  geom_boxplot(fill = '#278EA2') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle('%Study')

fix_box <-  ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_Histology.Stage')) + 
  geom_boxplot(fill = '#71482A') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle('%Histology+Stage')

# diver_box <-  ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Diversity')) + 
#   geom_boxplot(fill = '#FFE9D0')  + theme_bw() + 
#   theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


exprPer_box <- ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_expressed')) + 
  geom_boxplot(fill = '#007F61') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(),axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('%cluster expressed')


sample_box <- ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_Sample')) + 
  geom_boxplot(fill = '#324B51') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(),axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('%Sample')

cluster_box <-  ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_Cluster')) + 
  geom_boxplot(fill = '#95B0B7') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(),axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle('%Cluster')

residual_box <-  ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Per_Unexplained')) + 
  geom_boxplot(fill = 'grey') + xlim(0, xcf) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle('%Unexplained')

gini_box <- ggplot(aes(x = value, y = Category), data = module_dl %>% filter(variable == 'Gini_coefficient')) + 
  geom_boxplot(fill = '#95B0B7') + xlim(0, 1) + theme_bw() + 
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle('Gini coefficient')

#cowplot::plot_grid(study_box, fix_box, diver_box, ncol=3, rel_widths = c(1.7,1,1))


#cowplot::plot_grid(study_box, fix_box, sample_box, cluster_box, residual_box, gini_box, ncol=6, rel_widths = c(1.7,1,1,1,1))

cowplot::plot_grid(study_box, fix_box, sample_box, cluster_box, residual_box, ncol=5, rel_widths = c(1.7,1,1,1,1))

pdf("../../Figure/Figure2/Figure2_Variability_boxplot.pdf", width = 5, height = 3)
print(
  cowplot::plot_grid(study_box, fix_box, ncol=2, rel_widths = c(2.5,1))
)
dev.off()


pdf("../../Figure2/Figure2_Variability_boxplot_full.pdf", width = 9, height = 3)
print(
  cowplot::plot_grid(study_box, fix_box, sample_box, cluster_box, residual_box, ncol=5, rel_widths = c(2.5,1,1,1,1))
)
dev.off()

