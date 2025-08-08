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
library(ggpubr)
library(stringr)

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

#### original module score
y <- factorizeMatrix(data_sce, type = "counts")
cell_m_counts <- y[[1]][['cell']]; cell_m_counts <- t(cell_m_counts)
colSum_sizeFactor <- colSums(t(cell_m_counts)) / 10000 
normalized_m_counts <- sweep(cell_m_counts, 1, colSum_sizeFactor, "/")
module_activity <- log2(normalized_m_counts + 1)

#### load the corrected module score
cr_module <-  readRDS('../../Mixed_model_output/corrected_module_score.rds')

#### look at module heatmap of L44
moduleHeatmap(data_sce, featureModule = 44, useAssay = 'counts', altExpName = 'featureSubset', 
              topCells = 500, topFeatures = 20)


### 2. identify studies/samples that has high module L44 (CDK4/MDM2) expression
colData(data_sce) <- cbind(colData(data_sce), 
                           data.frame('L44_prob' = module_prob[, c('L44')],
                                      'L44_expr' = module_activity[, c('L44')],
                                      'L44_corrected_expr' = cr_module['L44', ]))


plotSCEViolinColData(data_sce, coldata = 'L44_expr', groupBy = 'celda_k')
plotSCEViolinColData(data_sce, coldata = 'L44_corrected_expr', groupBy = 'celda_k')





#### 2.1. Prepare for the plot
# create violin for 5 groups. CDK4+/- in two studies. 5th group for all other studies
ks <- unique(data_sce$celda_k)
l <- paste0('L', c(44))


cid <- colnames(data_sce)[data_sce$celda_k %in% ks]
d <- colData(data_sce)[cid, c('study','new_sample', 'celda_k', 'L44_expr')] %>% as.data.frame()


d$celda_k <- as.character(d$celda_k)
d$CDK4_status<- 'Neutral'
d$CDK4_status[d$new_sample %in% c('Wu_et-P2', 'Bischoff_et-p019t')] <- 'Amplified'
d$study[!d$study %in% c("Wu_et", "Bischoff_et")] <- 'Other_studies'

d$CDK4 <- paste(d$study, d$CDK4_status, sep = '-')
rd <- reshape2::melt(d, id = c('CDK4'), measure.vars = 'L44_expr')
colnames(rd) <- c('CDK4_status', 'Module', 'log2(Normalized score+1)')


rd$CDK4_status <- factor(rd$CDK4_status, levels = c('Other_studies-Neutral', 'Wu_et-Neutral', 'Bischoff_et-Neutral',  
                                                  'Wu_et-Amplified', 'Bischoff_et-Amplified'))
rd$`CDK4 CNV` <- str_split(rd$CDK4_status, '-', simplify = T)[,2] %>% factor(., levels = c('Neutral', 'Amplified'))
cluster_col <- setNames(c("#1F78B4", "red3"), #"#1F78B4", "#1F78B4", "red3",  
                        levels(rd$`CDK4 CNV`)
)


cdk4_ov <- ggplot(aes(x = CDK4_status, y = `log2(Normalized score+1)`, fill = `CDK4 CNV`), data = rd) +  
  geom_violin(position = position_dodge(0.8), width = 1) + 
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 45, hjust= 1, size = 8), 
        axis.text.x = element_blank(), 
        #axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        #legend.position="none",
        strip.text.x = element_text(size = 12)) + 
  scale_fill_manual(values = cluster_col)

cdk4_ov

### adding FC annotation
cdk4ov_map <- data.frame('CDK4_status' = rd$CDK4_status, 'x' = as.numeric(rd$CDK4_status)) %>% distinct(CDK4_status, .keep_all = T)
cdk4ov_sum <- rd %>% group_by(CDK4_status) %>% summarise(meanScore = mean(`log2(Normalized score+1)`), 
                                                        maxScore = max(`log2(Normalized score+1)`))
cdk4ov_sum$x <- cdk4ov_map[match(cdk4ov_sum$CDK4_status, cdk4ov_map$CDK4_status), 'x']
cdk4ov_sum <- as.data.frame(cdk4ov_sum); rownames(cdk4ov_sum) <- cdk4ov_sum$x

comparison <- list(c(4,5), c(2,4), c(3,5), c(1,2), c(1,3))


### plot
bar_data <- lapply(comparison, function(com) {
  
  fc <- (cdk4ov_sum[com[2], 'meanScore'] - cdk4ov_sum[com[1], 'meanScore'])
  
  tribble(
    ~group1, ~group2, ~label, ~y.position,
    com[1], com[2], paste0('log2FC = ', round(fc, 3)), max(cdk4ov_sum[com, 'maxScore']) + 0.3
  )
  
}) %>% do.call('rbind', .)

cdk4_ov <- cdk4_ov + geom_signif(y_position = c(9, 11, 12, 9, 8) + 0.5, xmin = bar_data$group1, 
                     xmax = bar_data$group2, annotation = bar_data$label,
                     tip_length = 0)




#### also plot the corrected module score
d <- colData(data_sce)[cid, c('study','new_sample', 'celda_k', 'L44_corrected_expr')] %>% as.data.frame()

d$celda_k <- as.character(d$celda_k)
d$CDK4_status<- 'Neutral'
d$CDK4_status[d$new_sample %in% c('Wu_et-P2', 'Bischoff_et-p019t')] <- 'Amplified'
d$study[!d$study %in% c("Wu_et", "Bischoff_et")] <- 'Other_studies'

d$CDK4 <- paste(d$study, d$CDK4_status, sep = '-')
rd <- reshape2::melt(d, id = c('CDK4'), measure.vars = 'L44_corrected_expr')
colnames(rd) <- c('CDK4_status', 'Module', 'log2(Corrected score+1)')


rd$CDK4_status <- factor(rd$CDK4_status, levels = c('Other_studies-Neutral', 'Wu_et-Neutral', 'Bischoff_et-Neutral',  
                                                    'Wu_et-Amplified', 'Bischoff_et-Amplified'))
rd$`CDK4 CNV` <- str_split(rd$CDK4_status, '-', simplify = T)[,2] %>% factor(., levels = c('Neutral', 'Amplified'))
cluster_col <- setNames(c("#1F78B4", "red3"), #"#1F78B4", "#1F78B4", "red3",  
                        levels(rd$`CDK4 CNV`)
)


cdk4_v <- ggplot(aes(x = CDK4_status, y = `log2(Corrected score+1)`, fill = `CDK4 CNV`), data = rd) +  #`Module score`
  geom_violin(position = position_dodge(0.8), width = 1) + #facet_wrap(~Study, nrow=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8), # vjust = -0.3, 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        #legend.position="none",
        strip.text.x = element_text(size = 12)) + 
  scale_fill_manual(values = cluster_col)

cdk4_v



### adding FC annotation
cdk4v_map <- data.frame('CDK4_status' = rd$CDK4_status, 'x' = as.numeric(rd$CDK4_status)) %>% distinct(CDK4_status, .keep_all = T)
cdk4v_sum <- rd %>% group_by(CDK4_status) %>% summarise(meanScore = mean(`log2(Corrected score+1)`), 
                                           maxScore = max(`log2(Corrected score+1)`))
cdk4v_sum$x <- cdk4v_map[match(cdk4v_sum$CDK4_status, cdk4v_map$CDK4_status), 'x']
cdk4v_sum <- as.data.frame(cdk4v_sum); rownames(cdk4v_sum) <- cdk4v_sum$x

comparison <- list(c(4,5), c(2,4), c(3,5), c(1,2), c(1,3))


### plot
bar_data <- lapply(comparison, function(com) {
  
  fc <- (cdk4v_sum[com[2], 'meanScore'] - cdk4v_sum[com[1], 'meanScore'])
  
  tribble(
    ~group1, ~group2, ~label, ~y.position,
    com[1], com[2], paste0('log2FC = ', round(fc, 3)), max(cdk4v_sum[com, 'maxScore']) + 0.3
    )
  
}) %>% do.call('rbind', .)
cdk4_v <- cdk4_v + geom_signif(y_position = c(9, 12, 13, 9, 8) + 0.5, xmin = bar_data$group1, 
            xmax = bar_data$group2, annotation = bar_data$label,
            tip_length = 0)


fig_path <- '../../Figure/Figure2'
fn <- file.path(fig_path, 'Figure2_CDK4_violin.pdf')
pdf(fn, width = 6, height = 8)
print(cowplot::plot_grid(cdk4_ov, cdk4_v, nrow = 2, rel_heights = c(0.8, 1)))
dev.off()




### 3. Select modules which has high % variability explained by study id (batch effect).
histo_m_summary <- read.csv(file.path('../../Mixed_model_output/model_summary_AllCells.csv'))
m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]


l <- 'L158'
ks <- unique(data_sce$celda_k)

#### identify studies/samples that has high module L158 expression
colData(data_sce) <- cbind(colData(data_sce), 
                           data.frame(
                                      'L158_expr' = module_activity[, l]))

#### also include the corrected module score
cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')

colData(data_sce) <- cbind(colData(data_sce), 
                           data.frame(
                             'L158_corrected_expr' = cr_module[l, colnames(data_sce)]))

##### 3.1 showing the module expression at the same study above
cid <- colnames(data_sce)[data_sce$celda_k %in% ks]
d <- colData(data_sce)[cid, c('study','new_sample', 'celda_k', 'L158_expr')] %>% as.data.frame()

d$celda_k <- as.character(d$celda_k)
d$CDK4_status<- 'Neutral'
d$CDK4_status[d$new_sample %in% c('Wu_et-P2', 'Bischoff_et-p019t')] <- 'Amp'
d$study[!d$study %in% c("Wu_et", "Bischoff_et")] <- 'Other_studies'

d$CDK4 <- paste(d$study, d$CDK4_status, sep = '-')
d <- reshape2::melt(d, id = c('CDK4'), measure.vars = 'L158_expr')
colnames(d) <- c('CDK4_status', 'Module', 'Module score')

d$CDK4_status <- factor(d$CDK4_status, levels = c('Other_studies-Neutral', 'Wu_et-Neutral', 'Bischoff_et-Neutral',  
                                                  'Wu_et-Amp', 'Bischoff_et-Amp'))
cluster_col <- setNames(c("#1F78B4", "#1F78B4", "#1F78B4", "#FF7E00",  "#FF7E00"),
                        levels(d$CDK4_status)
)


l158_v <- ggplot(aes(x = CDK4_status, y = `Module score`, fill = CDK4_status), data = d) + 
  geom_violin(position = position_dodge(0.8), width = 1) + #facet_wrap(~Study, nrow=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 10), # vjust = -0.3, 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        #legend.position="none",
        strip.text.x = element_text(size = 15)) + 
  scale_fill_manual(values = cluster_col)

l158_v


##### 3.2 plot it across study
cid <- colnames(data_sce)[data_sce$celda_k %in% ks]
d <- colData(data_sce)[cid, c('study','new_sample', 'celda_k', 'L158_corrected_expr')] %>% as.data.frame()

d <- reshape2::melt(d, id = c('study'), measure.vars = 'L158_corrected_expr')
colnames(d) <- c('Study', 'Module', 'Corrected score')

ggplot(aes(x = Study, y = `Corrected score`, fill = Study), data = d) + 
  geom_violin(position = position_dodge(0.8), width = 1) + #facet_wrap(~Study, nrow=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 10), # vjust = -0.3, 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position="none",
        strip.text.x = element_text(size = 15)) #+ 
  #scale_fill_manual(values = cluster_col)
