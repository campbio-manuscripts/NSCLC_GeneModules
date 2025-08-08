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

setwd(here::here("./Scripts/Figure"))

source('../../Scripts/Figure/Figure_utils.R')

### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

#### load corrected module score
cr_module <- readRDS('../../Mixed_model_output/corrected_module_score.rds')
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
histo_m_summary <- read.csv('../../Mixed_model_output/model_summary_Histology.csv')
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
                                   clusterLabelSize = 2,
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
  ) + labs(y= "UMAP_2", x = "UMAP_1")


### 5. Plot the sqaumous module L14
l <- 'L14'
t <- 'L14 (Sqiamous_1)'


col_scale <- c(0, 0.2, 0.4, 0.6, 0.995)
colors <- c('grey88', 'snow','orange','red', 'red3')

umap_ker <- plotSCEDimReduceFeatures(lusc_m, feature=l, reducedDimName = 'corrected_UMAP', useAssay = 'module',
                                     #groupBy = 'sample',
                                     legendTitle = 'Corrected\nScore',
                                     dotSize = 0.1, title = t, titleSize = 10, legendSize = 6, legendTitleSize = 6,
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

fn <- file.path('../../Figure/Figure6', 'Figure7_LUSC_Squamous1K14_UMAP.pdf')
pdf(file = fn, width = 3, height = 3)
print(umap_ker)
dev.off()





################
# use a different way to define whether keritinization module is on
################
kert_d <- data.frame('L14' = assay(lusc_m)['L14',])

pdf("../../Figure/Figure6/Figure7_LUSC_L14Squamous_density.pdf", 
    width = 2, height = 1.5)
ggplot(kert_d, aes(x=L14)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill = 'white')+
  geom_density(alpha=.3) + xlab('L14 Score') + 
  geom_vline(xintercept=1, colour="blue") + 
  geom_vline(xintercept=5, colour="red3") +
  
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text = element_text (size = 8), 
        axis.title = element_text (size = 9))
dev.off()


lusc_m$squamous_status <- 'Low'
lusc_m$squamous_status[assay(lusc_m)['L14',] >= 1] <- 'Intermediate'
lusc_m$squamous_status[assay(lusc_m)['L14',] >= 5] <- 'High'

# lusc_m$squamous_status <- ifelse(assay(lusc_m)['L14',] >= 1, 'High', 'Low') #0.5

# lusc_m$keriti_status <- 'Low'
# lusc_m$keriti_status[assay(lusc_m)['L17',] >= 0.5] <- 'Intermediate'
# lusc_m$keriti_status[assay(lusc_m)['L17',] >= 8] <- 'High'



lusc_meta <- colData(lusc_m) %>% as.data.frame()
squa_order <- lusc_meta %>% group_by(sample) %>% summarise(Squamous_low = mean(squamous_status == 'Low')) %>% arrange(desc(Squamous_low)) #Low
squa_order$order <- 1:nrow(squa_order)

squa_meta <- as.data.frame(table(lusc_meta$sample, lusc_meta$squamous_status))
colnames(squa_meta) <- c('Sample', "Squamous", 'Freq')
squa_meta$order <- squa_order[match(squa_meta$Sample, squa_order$sample), 'order'] %>% unlist()


fn <- file.path('../../Figure/Figure6', 'Figure7_LUSC_Squamous_barplot_L14.pdf')
pdf(file = fn, width = 2, height = 1.5)
ggplot(aes(x = order, y = Freq, fill = Squamous ), data = squa_meta) +  
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


#### quantify % squamous-low, intermediate or high cells in each sample
lusc_squamous_summary <- lusc_meta %>% group_by(sample) %>% summarize(
  low_pert = mean(squamous_status == 'Low'), 
  mid_pert = mean(squamous_status == 'Intermediate'), 
  high_pert = mean(squamous_status == 'High'), 
  
)
 
sum(lusc_squamous_summary$low_pert > 0 &lusc_squamous_summary$high_pert > 0 )

