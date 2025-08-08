### figure 5, LUSC specific umap with corrected module scores

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


source('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Scripts/Figure/Figure_utils.R')

### 1. Load dataset
## load sce object
op <- '/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

#### load corrected module score
cr_module <- readRDS('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/corrected_module_score.rds')
module_score <- t(cr_module)


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
lusc_id <- colnames(data_sce)[data_sce$celda_k %in% k & data_sce$Histology %in% his]

### 3. create LUSC sce object
histo_m_summary <- read.csv('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/model_summary_Histology.csv')
skip_m <- histo_m_summary$module[histo_m_summary$study_VarPer >= 0.3]


lusc_m <- SingleCellExperiment(assay = list(module = cr_module[!rownames(cr_module) %in% skip_m,  #lm_residual
                                                                        lusc_id]),
                                        rowData = NULL,
                                        colData = colData(data_sce)[lusc_id, ])


lusc_m <- singleCellTK::getUMAP(lusc_m, useAssay = 'module',
                                reducedDimName = 'corrected_UMAP', useReducedDim = NULL, 
                                sample = NULL, minDist = 0.2, seed = 2022, logNorm = F)



### 4. Plot UMAP with sample label
l <- 'new_sample'
st_umap <- plotSCEDimReduceColData(lusc_m, 
                        colorBy = l, reducedDimName = 'corrected_UMAP',
                        titleSize = 10,
                        colorScale = NULL,
                        axisSize = 0,
                        axisLabelSize = F,
                        dotSize = 1, 
                        baseSize = 3.5,
                        legendTitleSize = 7,
                        clusterLabelSize = 0.6,
                        labelClusters = T) +  theme_bw() +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
        axis.text.y = element_blank(),

        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position= 'none', #c(0.2, 0.3)
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'),
        legend.text = element_text(size=13)
  ) + labs(y= "UMAP_2", x = "UMAP_1") #+ ggrepel::geom_text_repel(max.overlaps = Inf) 


### 5. Plot the Keritinization module L31
l <- 'L31'
t <- 'L31 (Keratinization)'

col_scale <- c(0, 0.2, 0.4, 0.6, 0.995)
colors <- c('grey88', 'snow','orange','red', 'red3')

# fig_path <- '/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure1'
# fn <- file.path(fig_path, paste0('Figure1_module_',l,'_UMAP.tiff'))
# tiff(fn, width = 400, height = 350)

umap_ker <- plotSCEDimReduceFeatures(lusc_m, feature=l, reducedDimName = 'corrected_UMAP', useAssay = 'module',
                         #groupBy = 'sample',
                         legendTitle = 'Corrected\nScore',
                         dotSize = 0.1, title = t, titleSize = 5, legendSize = 5, legendTitleSize = 5,
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
# dev.off()

umap_ker$layers[[2]] <- st_umap$layers[[3]]
umap_ker$layers[[2]]$constructor$`max.overlap` = 100 ### hard-code the max.overlap here to show all label

fn <- file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure6', 'Figure7_LUSC_Keratinization_UMAP.pdf')
pdf(file = fn, width = 4, height = 4)
print(umap_ker)
dev.off()



### 6. check expression of Keritinization module L31
# #### module probability
# x <- factorizeMatrix(data_sce[, lusc_id], type = "proportion")
# cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
# lusc_module_prob <- cell_probs
# 
# l <- 'L31'
# plot(density(lusc_module_prob[, l]))


# ### subset to K68, 75, 79
# subk <- c(68, 75, 79)
# ker_sce_m <- lusc_m[, lusc_m$celda_k %in% subk]
#   
# ker_sce_m <- singleCellTK::getUMAP(ker_sce_m, useAssay = 'module',
#                         reducedDimName = 'corrected_UMAP', useReducedDim = NULL, 
#                         sample = NULL, minDist = 0.01, seed = 2022, logNorm = F)
# 
# 
# 
# ### plot cluster
# l <- 'new_sample'
# l <- 'celda_k'
# plotSCEDimReduceColData(ker_sce_m, 
#                         colorBy = l, reducedDimName = 'corrected_UMAP',
#                         titleSize = 10,
#                         colorScale = NULL,
#                         axisSize = 0,
#                         axisLabelSize = F,
#                         dotSize = 1, 
#                         baseSize = 10,
#                         legendTitleSize = 10,
#                         clusterLabelSize = 3,
#                         labelClusters = T) +  theme_bw() +
#   theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust= 1)
#         axis.text.y = element_blank(),
#         
#         axis.ticks.x=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position= 'none', #c(0.2, 0.3)
#         legend.key.size = unit(1, 'cm'), #change legend key size
#         legend.key.height = unit(0.6, 'cm'),
#         legend.text = element_text(size=13)
#   ) + labs(y= "UMAP_2", x = "UMAP_1")
# 
# 
# 
# ### plot the modules
# l <- 'L31'
# 
# col_scale <- c(0, 0.2, 0.4, 0.6, 0.995)
# colors <- c('white', 'snow','orange','red', 'red3')
# 
# plotSCEDimReduceFeatures(ker_sce_m, feature=l, reducedDimName = 'corrected_UMAP', useAssay = 'module',
#                          legendTitle = 'Corrected\nScore',
#                          dotSize = 0.1, title = t, titleSize = 28, legendSize = 10, legendTitleSize = 15) + 
#   scale_color_gradientn(values = col_scale, 
#                         colours = colors,
#                         na.value = "grey50") + 
#   theme(
#     axis.title = element_blank(),         # Change both x and y axis titles
#     axis.text=element_blank(),
#     axis.ticks = element_blank(),
#     legend.key.size = unit(5, 'mm'), legend.margin=margin(0,0,0,0)
#   )
# 
# 
# 
# ### count % of cells in K68, 75, 79. Use that as the intra-tumor heterogeneity
# lusc_meta <- colData(lusc_m) %>% as.data.frame()
# keri_per <- lusc_meta %>% group_by(sample) %>% summarise(per_Keri = mean(celda_k %in% subk)) %>% as.data.frame
# keri_per <- keri_per %>% arrange((per_Keri))
# ggplot(aes(x = 1:nrow(keri_per), y = per_Keri, fill = per_Keri), d = keri_per) + geom_bar(stat="identity", fill="steelblue") +
#   theme_bw() + 
#   theme(
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     #axis.title = element_blank(),         # Change both x and y axis titles
#     strip.text.x = element_text(size = 20),
#     legend.key.size = unit(1.2, 'cm'),
#     legend.title = element_text(size=14), #change legend title font size
#     legend.text = element_text(size=12)
#     #panel.spacing = unit(0.1, "lines")
#     #legend.position = 'none',
#   ) + ylab("% Keritinization cells") + xlab('Sample')


################
# use a different way to define whether keritinization module is on
################
kert_d <- data.frame('L31' = assay(lusc_m)['L31',])

pdf("/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure7/Figure7_LUSC_L31Keratinazation_density.pdf", 
    width = 2, height = 1.5)
ggplot(kert_d, aes(x=L31)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill = 'white')+
  geom_density(alpha=.3) + xlab('L31 Score') + 
  geom_vline(xintercept=0.5, colour="blue") + 
  geom_vline(xintercept=4.5, colour="red3") + 
  
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text = element_text (size = 8), 
        axis.title = element_text (size = 9))
dev.off()


lusc_m$keriti_status <- ifelse(assay(lusc_m)['L31',] >= 0.5, 'Intermediate', 'Low')
lusc_m$keriti_status[assay(lusc_m)['L31',] > 4.5] <- 'High'

lusc_meta <- colData(lusc_m) %>% as.data.frame()
keri_order <- lusc_meta %>% group_by(sample) %>% summarise(Keriti_low = mean(keriti_status == 'Low')) %>% arrange(desc(Keriti_low))
keri_order$order <- 1:nrow(keri_order)

keri_meta <- as.data.frame(table(lusc_meta$sample, lusc_meta$keriti_status))
colnames(keri_meta) <- c('Sample', "Keritinization", 'Freq')
keri_meta$order <- keri_order[match(keri_meta$Sample, keri_order$sample), 'order'] %>% unlist()


fn <- file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure7', 'Figure7_LUSC_Keratinization_barplot_LowMedHigh.pdf')
pdf(file = fn, width = 2, height = 1.5)
ggplot(aes(x = order, y = Freq, fill = Keritinization), data = keri_meta) +  
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y = element_text (size = 4),
    axis.title = element_blank(),         # Change both x and y axis titles
    strip.text.x = element_text(size = 4),
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size=3), #change legend title font size
    legend.text = element_text(size=3),
    legend.margin=margin(c(0,0,0,0)),
    plot.margin=unit(x=c(0.3,0.3,0.3,0.3),units="mm"),
    # legend.spacing.x = unit(0.5, 'cm')
    # panel.spacing = unit(0.1, "lines")
    #legend.position = 'none',
  ) +     scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + 
  
  scale_fill_manual(breaks = c("Low", "Intermediate", "High"),
                    values=c("#3C99B1", "#F1BB7B", "#C93428"))
dev.off()



#### just checking percent of High or intermediate cells
keri_per <- lusc_meta %>% group_by(sample) %>% summarise(per_Keri = mean(keriti_status %in% c('High', 'Intermediate'))) %>% as.data.frame
keri_per <- keri_per %>% arrange(per_Keri)
summary(keri_per$per_Keri)

# fn <- file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure7', 'Figure7_LUSC_Keratinization_barplot.pdf')
# pdf(file = fn, width = 7, height = 4)
# ggplot(aes(x = 1:nrow(keri_per), y = per_Keri, fill = per_Keri), d = keri_per) + geom_bar(stat="identity", fill="red4") +
#   theme_bw() + 
#   theme(
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     #axis.title = element_blank(),         # Change both x and y axis titles
#     strip.text.x = element_text(size = 20),
#     legend.key.size = unit(1.2, 'cm'),
#     legend.title = element_text(size=14), #change legend title font size
#     legend.text = element_text(size=12)
#     #panel.spacing = unit(0.1, "lines")
#     #legend.position = 'none',
#   ) + ylab("% Keritinization cells") + xlab('Sample') + 
#   scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) + 
#   scale_x_continuous(expand = c(0, 0))
# dev.off()






### look at within-sample heterogeneity
s <- 'Wu_et-P17'

s <- 'Goveia_et-41_tumor_primary'
p6_m <- lusc_m[, lusc_m$sample == s]
p6_m <- p6_m[target_modules, ]

p6_m <-  singleCellTK::getUMAP(p6_m, useAssay = 'module',
                               reducedDimName = 'corrected_UMAP', useReducedDim = NULL, 
                               sample = NULL, minDist = 0.01, seed = 2022, logNorm = F)

l <- 'L31'
col_scale <- c(0, 0.2, 0.4, 0.6, 0.995)
colors <- c('steelblue', 'snow','orange','red', 'red3')
plotSCEDimReduceFeatures(p6_m, feature=l, reducedDimName = 'corrected_UMAP', useAssay = 'module',
                         
                         legendTitle = 'Corrected\nScore',
                         dotSize = 0.3, title = t, titleSize = 28, legendSize = 10, legendTitleSize = 15) + 
  scale_color_gradientn(values = col_scale, 
                        colours = colors,
                        na.value = "grey50") + 
  theme(
    axis.title.y = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(5, 'mm'), legend.margin=margin(0,0,0,0)
  ) + xlab(s)




l <- 'keriti_status'
plotSCEDimReduceColData(p6_m, 
                        colorBy = l, reducedDimName = 'corrected_UMAP',
                        titleSize = 10,
                        colorScale = NULL,
                        axisSize = 3,
                        axisLabelSize = 2,
                        dotSize = 0.3, 
                        baseSize = 10,
                        legendTitleSize = 10,
                        clusterLabelSize = 3,
                        labelClusters = T) 




l <- 'L31'
plotSCEViolinAssayData(lusc_m, feature = l, useAssay = 'module', groupBy = 'celda_k') 


