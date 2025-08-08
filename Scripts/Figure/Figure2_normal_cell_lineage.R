### this script is used to look at healthy tissue
library(celda)
library(singleCellTK)
library(ggplot2)
library(dplyr)
library(eulerr)
library(data.table)
library(stringr)

setwd(here::here("./Scripts/Figure"))

### 1. load data
fp <- '../../Data'
nor_sce <- readRDS(file.path(fp, 'Basil_Nature_2022_celda_cg_decontXcounts_min3_filter_mito10_detected500_sum750_sub_K40_L125.rds'))
nor_sce <- runNormalization(nor_sce, outAssayName = 'logNormalized', normalizationMethod = 'LogNormalize')

### 2. plot UMAP with celda cluster
nor_sce$celda_k <- colData(altExp(nor_sce))$celda_cell_cluster
reducedDim(nor_sce, 'celda_umap') <- reducedDim(altExp(nor_sce), 'celda_UMAP')

plotSCEDimReduceColData(nor_sce, 
                        colorBy = 'celda_k', reducedDimName = 'celda_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = T, 
                        colorScale = NULL, title = "celda_k", titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank()
    #legend.position="none"
  ) 

### 3. double check the cell type identity
g <- 'ENSG00000168484_SFTPC'
t <- 'SFTPC'

t <- 'SFTPB'
g <- 'ENSG00000168878_SFTPB'

t <- 'AGER'
g <- 'ENSG00000204305_AGER'

t <- 'SCGB3A2'
g <- 'ENSG00000164265_SCGB3A2'

t <- 'SCGB1A1'
g <- 'ENSG00000149021_SCGB1A1'


t <- 'SCGB3A1'
g <- 'ENSG00000161055_SCGB3A1'

t <- 'MUC5AC'
g <- 'ENSG00000215182_MUC5AC'

t <- 'MUC5B' 
g <- 'ENSG00000117983_MUC5B'

t <- 'PTPRC'
g <- 'ENSG00000081237_PTPRC'

t <- 'SFTPA1'
g <- 'ENSG00000122852_SFTPA1'
rowData(nor_sce) %>% as.data.frame() %>% filter(Symbol == t)

col_scale <- c(0, 0.3, 0.6, 0.8, 0.995)
colors <- c('white', 'snow','orange','red', 'red3')


plotSCEDimReduceFeatures(nor_sce, feature=g, reducedDimName = 'celda_umap', useAssay = 'logNormalized',
                         legendTitle = 'Expression',
                         dotSize = 0.1, title = t, titleSize = 28, legendSize = 10, legendTitleSize = 15) + 
  scale_color_gradientn(values = col_scale, 
                        colours = colors,
                        na.value = "grey50") + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(5, 'mm'), legend.margin=margin(0,0,0,0)
  )


plotSCEViolinAssayData(nor_sce, feature = g, useAssay = "logNormalized", groupBy = "celda_k")


## 4. Get the list of module gene
nor_ta <- featureModuleTable(nor_sce)

gt <- c('SFTPC', 'SFTPD','SFTPB', 'SFTPA1', 'SFTPA2', 'SCGB3A1', 'SCGB3A2', "AGER", "AGR3", "HOPX")
gs <- rownames(nor_sce)[match(gt, rowData(nor_sce)$Symbol)]
names(gs) <- gt
featureModuleLookup(nor_sce, features = gs)


#### here is the list of modules that our original lineage module is in
#### SFTPC: 53
#### SFTPD: 56
#### SFTPB: 57
#### SFTPA1/2: 55
#### SCGB3A1: 80
#### SCGB3A2: 17



# ## 5. First, look at AT2 cell and SCGB3A2+/SCGB1A1- population
# nor_at2_scgb3a2 <- nor_sce[, nor_sce$celda_k %in% as.character(c(11:20, 22, 25, 30, 31))]
# nor_at2_scgb3a2$celda_k <- as.character(nor_at2_scgb3a2$celda_k)
# nor_at2_scgb3a2$cell_type <- 'AT2'; nor_at2_scgb3a2$cell_type[nor_at2_scgb3a2$celda_k == 25] <- 'SCGB3A1+/SCGB1A1-'
# 
# plotSCEDimReduceColData(nor_at2_scgb3a2, 
#                         colorBy = 'cell_type', reducedDimName = 'celda_umap', 
#                         clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = T, 
#                         colorScale = NULL, title = NULL, titleSize = 8) + 
#   theme(
#     axis.title = element_blank(),         # Change both x and y axis titles
#     axis.text=element_blank(),
#     axis.ticks = element_blank()
#     #legend.position="none"
#   ) 
# 
# 
# ### 5.1. Compute the probability matrix
# x <- factorizeMatrix(nor_sce, type = "proportion")
# cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
# module_prob <- cell_probs
# 
# ### 5.2. Compute the on/off of some modules
# cf <- 0.001 ##0.001
# nor_on_mat <- 1* t(module_prob[colnames(nor_at2_scgb3a2), ] > cf) #target_module
# lm_on <- t(nor_on_mat) %>% as.data.frame()
# 
# 
# lm_on <- lm_on %>% mutate(
#   SFTPC = 1*(rowSums(lm_on[, paste0('L', 53), drop=F]) != 0),
#   `SFTPA1/2` = 1*(rowSums(lm_on[, paste0('L', 55), drop=F]) != 0),
#   `SFTPB` = 1*(rowSums(lm_on[, paste0('L', 57), drop=F]) != 0),
#   `SFTPD` = 1*(rowSums(lm_on[, paste0('L', 56), drop=F]) != 0),
#   SCGB3A2 = 1*(rowSums(lm_on[, paste0('L', 17), drop=F]) != 0),
#   SCGB1A1 = 1*(rowSums(lm_on[, paste0('L', 80), drop=F]) != 0)
# )
# lm_on <- lm_on[, -grep('^L', names(lm_on))]
# 
# 
# 
# ### 5.3 plot eulerr plot
# lm_on$cell_type <- NULL
# m_x <- lapply(lm_on, function(l) {
#   rownames(lm_on)[l == 1]
# })
# 
# library(eulerr)
# m <- grep('^L',  colnames(lm_on), invert = T, value=T)
# 
# md <- m_x[m]
# set.seed(2020)
# fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,
# 
# p <- plot(fit,
#      quantities = list(type = c("counts"), cex = 1), #TRUE #, "counts", "percent"
#      #quantities = T,
#      fill = RColorBrewer::brewer.pal(length(m), "Set3"),
#      lty = 1:length(m),
#      #labels = list(cex = 3),
#      #adjust_labels = T,
#      #labels = NULL,
#      legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"),
#                    labels = rownames(fit$ellipses),
#                    cex = 1))
# p


# # Create a grid grob
# grid.grob <- grobTree(textGrob("Hello, world!", x = 0.5, y = 0.5))
# 
# p
# # Add text labels
# geom_text_repel(data = data.frame(x = 0.5, y = 0.5, label = "Hello, world!"), aes(x = x, y = y, label = label))
# 
# grid.draw(
#   
# )
# # Draw the grid grob
# grid.draw(grid.grob)



# check what intersections are missing
# miss_op <- data.frame('OP' = names(fit$fitted.values)[fit$fitted.values == 0],
#                       'NumOP' = fit$original.values[fit$fitted.values == 0])
# miss_op$OP_Per <- miss_op$NumOP / length(Reduce('union', m_x[m])) * 100
# miss_op %>% filter(OP_Per != 0) %>% arrange(OP_Per)




### simplify the eulerr plot
# overall_fit <- fit$original.values[fit$original.values >=50]
# of <- euler(overall_fit, quantities = TRUE,  shape = "ellipse") #,
# 
# p <- plot(of,
#      quantities = list(type = c("counts", 'percent'), cex = 0.8), #TRUE #, "counts", "percent"
#      #quantities = T,
#      fill = list(RColorBrewer::brewer.pal(nrow(of$ellipses), "Set2"), alpha = 1),
#      lty = 1:nrow(of$ellipses),
#      labels = list(cex = 2),
#      #labels = NULL,
#      legend = list(cols = RColorBrewer::brewer.pal(nrow(of$ellipses), "Set2"),
#                    labels = rownames(of$ellipses),
#                    cex = 1.5))
# p

# check what intersections are missing
# miss_op <- data.frame('OP' = names(of$fitted.values)[of$fitted.values == 0],
#                       'NumOP' = of$original.values[of$fitted.values == 0])
# miss_op$OP_Per <- miss_op$NumOP / length(Reduce('union', m_x[m])) * 100
# of_missOP <- miss_op %>% filter(OP_Per != 0) %>% arrange(OP_Per)
# sum(of_missOP$NumOP); sum(of_missOP$OP_Per)


### plot the missing number and the overlap
# missof <- euler(setNames(of_missOP$NumOP, rownames(of_missOP)), quantities = TRUE,  shape = "ellipse")
# plot(missof, quantities = list(type = c("counts"), cex = 1.3),
#      fill = list(RColorBrewer::brewer.pal(nrow(missof$ellipses), "Set2"), alpha = 1),
#      lty = 1:nrow(missof$ellipses),
#      labels = list(cex = 2),
#      legend = list(cols = RColorBrewer::brewer.pal(nrow(missof$ellipses), "Set2"),
#                    labels = rownames(missof$ellipses),
#                    cex = 1.5))

# library(UpSetR)
# 
# p <- upset(fromExpression(setNames(of_missOP$NumOP, rownames(of_missOP))), 
#       #nintersects = 40, 
#       nsets = 10, 
#       order.by = "degree", 
#       decreasing = T, 
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0, 
#       text.scale = 1.1, 
#       point.size = 2.8, 
#       line.size = 1
# )



### 5.4 check % expression on each cell type
# lm_on_check <- lm_on
# lm_on_check$cell_type <- nor_at2_scgb3a2$cell_type
# lm_on_check$'SCGB3A2'
# 
# l <- 'SFTPC'
# lm_on_l <- reshape2::melt(table(lm_on_check[[l]], lm_on_check$cell_type))
# ggplot(aes(x = Var2, y = value, fill = as.character(Var1)), data = lm_on_l) + geom_bar(stat = 'identity', position = 'fill') + ggtitle(l)




### 6. A lot of overlapping. Plot the module expression on the umap
# y <- factorizeMatrix(nor_sce, type = "counts")
# cell_m_counts <- y[[1]][['cell']]; cell_m_counts <- t(cell_m_counts)
# colSum_sizeFactor <- colSums(t(cell_m_counts)) / 10000 
# normalized_m_counts <- sweep(cell_m_counts, 1, colSum_sizeFactor, "/")
# module_activity <- log2(normalized_m_counts + 1)
# 
# #### SFTPC: 53
# #### SFTPD: 56
# #### SFTPB: 57
# #### SFTPA1/2: 55
# #### SCGB3A1: 80
# #### SCGB3A2: 17
# 
# 
# nor_at2_scgb3a2_m <- SingleCellExperiment(assay = list(module = t(module_activity[colnames(nor_at2_scgb3a2), ])),
#                                           rowData = NULL, 
#                                           colData = colData(nor_at2_scgb3a2), 
#                                           reducedDims = reducedDims(nor_at2_scgb3a2))
# 
# 
# col_scale <- c(0, 0.3, 0.6, 0.8, 0.995)
# colors <- c('white', 'snow','orange','red', 'red3')
# 
# t <- 'L17'
# 
# plotSCEDimReduceFeatures(nor_at2_scgb3a2_m, feature=t, reducedDimName = 'celda_umap', useAssay = 'module',
#                          legendTitle = 'Expression',
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







### 7. Include normal AT1, AT2 cells together. Check the expression of our nsclc_meta AT1 and AT2 modules / markers
nor_at2at1 <- nor_sce[, nor_sce$celda_k %in% as.character(c(11:20, 22, 25, 30, 31, 23, 24))]
nor_at2at1$celda_k <- as.character(nor_at2at1$celda_k)
nor_at2at1$cell_type <- 'AT2'
nor_at2at1$cell_type[nor_at2at1$celda_k == 25] <- 'SCGB3A1+/SCGB1A1-'
nor_at2at1$cell_type[nor_at2at1$celda_k %in% c(23, 24)] <- 'AT1'



#### 7.1 marker level
lineage_marker <- c('NKX2-1', 'RNASE1', 'NAPSA', 
                    'SFTPC', 'SFTPB', 'SFTPD', 'SFTPA1', 'SFTPA2', "SCGB3A2", "SCGB3A1",
                    'HOPX', 'AGER', 'AGR3')
lineage_gmar <- rowData(nor_sce) %>% as.data.frame() %>% filter(Symbol %in% lineage_marker)

gene_assay <- assay(nor_at2at1, 'counts')[rownames(lineage_gmar), ] %>% as.data.frame()
rownames(gene_assay) <- lineage_gmar$Symbol

gcf <- 5
gene_on <- 1*(gene_assay > gcf)

lm_on <- data.frame(
  AT2 = 1*(colSums(gene_on[c('SFTPC', 'SFTPB', 'SFTPD', 'SFTPA1', 'SFTPA2', "SCGB3A2", "SCGB3A1"), ]) > 0),
  AT1 = 1*(colSums(gene_on[c('HOPX', 'AGER', 'AGR3'), ]) > 0),  # 
  Alveolar = 1*(colSums(gene_on[c('NKX2-1', 'RNASE1', 'NAPSA'), ]) > 0)
)


m_x <- lapply(lm_on, function(l) {
  rownames(lm_on)[l == 1]
})

set.seed(2020)
mar_fit <-euler(m_x, quantities = TRUE,  shape = "circle", loss  = 'square', loss_aggregator = 'max') #, shape = "ellipse",

plot(mar_fit,
     quantities = list(type = c("counts"), cex = 1), #TRUE #, "counts", "percent"
     #quantities = F, 
     fill = RColorBrewer::brewer.pal(length(m_x), "Set2"),
     lty = 1:length(m_x),
     labels = list(cex = 1.5),
     #labels = NULL,
     legend = list(cols = RColorBrewer::brewer.pal(length(m_x), "Set2"), 
                   labels = rownames(mar_fit$ellipses)))



#### 7.2 use the module define in this study first. 
x <- factorizeMatrix(nor_at2at1, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs

lineage_mod <- featureModuleLookup(nor_at2at1, features = rownames(lineage_gmar)) ### look up the module list
names(lineage_mod) <- lineage_gmar$Symbol

# cf <- 0.01
cf <- 0.005

lm_on <- 1* t(module_prob[, unique(lineage_mod)] > cf) #target_module
lm_on <- t(lm_on) %>% as.data.frame()


lm_on <- lm_on %>% mutate(
  AT2 = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SFTPB", "SFTPC", "SFTPA2", "SFTPA1", "SFTPD", "SCGB3A2", "SCGB3A1")] %>% unique)]) != 0), 
  AT1 = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("AGER", "AGR3")]%>% unique)]) != 0), 
  Alveolar = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("RNASE1", "NKX2-1", "NAPSA")]%>% unique), drop=F]) != 0)
)



m_x <- lapply(lm_on, function(l) {
  rownames(lm_on)[l == 1]
})
m <- c('AT2', 'AT1', 'Alveolar')
md <- m_x[m]
fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,


plot(fit,
     quantities = list(type = c("counts"), cex = 1), #TRUE #, "counts", "percent"
     fill = RColorBrewer::brewer.pal(length(m), "Set2"),
     lty = 1:length(m),
     labels = list(cex = 1.5),
     #labels = NULL,
     legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"), 
                   labels = rownames(fit$ellipses)))



##### or separate for each module
lm_on <- lm_on %>% mutate(
  `SFTPA1/2` = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SFTPA2", "SFTPA1")] %>% unique), drop = F]) != 0), 
  SFTPC = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SFTPC")]%>% unique), drop = F]) != 0), 
  SFTPB = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SFTPB")]%>% unique), drop = F]) != 0), 
  SFTPD = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SFTPD")]%>% unique), drop = F]) != 0), 
  AGER = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("AGER")]%>% unique), drop = F]) != 0), 
  
  SCGB3A2 = 1*(rowSums(lm_on[, paste0('L', lineage_mod[c("SCGB3A2"), drop = F]%>% unique), drop=F]) != 0)
)


m_x <- lapply(lm_on, function(l) {
  rownames(lm_on)[l == 1]
})
m <- c("SFTPA1/2", "SFTPC", "SFTPB", "SFTPD", "AGER", "SCGB3A2")
md <- m_x[m]
fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,

pdf("../../Figure/Figure4/Figure4_Lung_lineageVenn_probbased_cf0.005.pdf", width = 4, height = 4)
plot(fit,
     quantities = list(type = c("counts"), cex = 0.5), #TRUE #, "counts", "percent"
     fill = RColorBrewer::brewer.pal(length(m), "Set2"),
     lty = 1:length(m),
     labels = list(cex = 0.5),
     #labels = NULL,
     legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"), 
                   labels = rownames(fit$ellipses), 
                   cex = 0.5))
dev.off()

##### investigate whether this cutoff = 0.01 is too high. Check the density plot of probability
cf <- 0.01
cf <- 0.005

prob_meta <- module_prob[colnames(nor_at2at1), 
                         paste0('L', lineage_mod[c("SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "SFTPD", "AGER", "SCGB3A2")] %>% unique)] %>% as.data.frame()
colnames(prob_meta) <- c("SFTPA1/2", "SFTPC", "SFTPB", "SFTPD", "AGER", "SCGB3A2")
prob_meta <- reshape2::melt(prob_meta)
colnames(prob_meta) <- c('Module', 'Probability')
ggplot(aes(x = Probability, color = Module), data = prob_meta) + geom_density() + 
  geom_vline(xintercept = cf, color = 'red3') + xlim(0, 0.25) + ylim(0, 300)

lapply(colnames(prob_meta), function (l) {
  ggplot(aes(x = Probability, color = Module), data = prob_meta) + geom_density() + 
    geom_vline(xintercept = cf, color = 'red3') + xlim(0, 0.25) + ylim(0, 300)
})




### 8. check specificity of these marker gene in AT1/AT2 cells
g <- "AGR3"
g <- "SFTPB"
plotSCEViolinAssayData(nor_at2at1, gs[g], useAssay = 'logNormalized', groupBy = "cell_type", title = g)
plotSCEViolinAssayData(nor_at2at1, gs[g], useAssay = 'counts', groupBy = "cell_type", title = g)



### 9. Perform DE analysis across these cell types. 
plotSCEDimReduceColData(nor_at2at1, 
                        colorBy = 'celda_k', reducedDimName = 'celda_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = T, 
                        colorScale = NULL, title = "celda_k", titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank()
    #legend.position="none"
  ) 

nor_at2at1$celltype <- 'AT2'
nor_at2at1$celltype[nor_at2at1$celda_k %in% c(25)] <- 'Transient'
nor_at2at1$celltype[nor_at2at1$celda_k %in% c(23, 24)] <- 'AT1'


nor_at2at1 <- runFindMarker(nor_at2at1, 
                            useAssay = 'logNormalized', 
                            cluster = 'celltype', 
                            log2fcThreshold = 0.5)

nor_at2at1_cm <- metadata(nor_at2at1)$findMarker

### 10. Overlapping between DE genes of each cell type and the module
nor_ta <- featureModuleTable(nor_at2at1)

gt <- c('SFTPC', 'SFTPD','SFTPB', 'SFTPA1', 'SFTPA2', 'SCGB3A1', 'SCGB3A2', "AGER") #"AGR3", "HOPX"
gs <- rownames(nor_at2at1)[match(gt, rowData(nor_at2at1)$Symbol)]
names(gs) <- gt
g_mods <- featureModuleLookup(nor_at2at1, features = gs)


#### here is the list of modules that our original lineage module is in
#### SFTPC: 53
#### SFTPD: 56
#### SFTPB: 57
#### SFTPA1/2: 55
#### SCGB3A1: 80
#### SCGB3A2: 17
#### AGER: 12

module_ls <- lapply(g_mods, function(l) {
  gl <- nor_ta[, l]; gl <- gl[gl != '']
  gl
})
names(module_ls) <- names(g_mods)
lapply(module_ls, length)


#### 10.1 Most gene modules are pretty small. Maybe need to check for SFTPD, and the SCGB3A1
at2_mark <- nor_at2at1_cm %>% filter(FDR < 0.05, Log2_FC > 0.5, celltype == 'AT2') %>% select(Gene) %>% unlist()
at1_mark <- nor_at2at1_cm %>% filter(FDR < 0.05, Log2_FC > 0.5, celltype == 'AT1') %>% select(Gene) %>% unlist()
transient_mark <- nor_at2at1_cm %>% filter(FDR < 0.05, Log2_FC > 0.5, celltype == 'Transient') %>% select(Gene) %>% unlist()


intersect(at2_mark, module_ls[['ENSG00000133661_SFTPD']]) ### 8 / 38 module genes in SFTPD modules have overlapped with AT2 marker
intersect(at2_mark, module_ls[['ENSG00000161055_SCGB3A1']]) ##0
intersect(at2_mark, module_ls[['ENSG00000164265_SCGB3A2']]) ##0

intersect(transient_mark, module_ls[['ENSG00000161055_SCGB3A1']]) ## 1/ 10 module genes overlapped with transient marker (SCGB3A1)
intersect(transient_mark, module_ls[['ENSG00000164265_SCGB3A2']]) ## 1/ 10 module genes overlapped with transient marker (SCGB3A1)



### 11. Use the marker defined in normal tissue, check the modules defined in tumor dataset
tumor_module_l <- fread("../../Data/Final_module_list.txt") %>% as.data.frame()
gm <- c(1,2,3,4,9,26,28)
names(gm) <- c('SFTPC', 'SFTPA1/2', 'SFTPB', 'SFTPD', 'AGER', 'SCGB3A2', 'SCGB3A1')

tumor_module_ls <- lapply(gm, function(l) {
  l <- paste0('L', l)
  gl <- tumor_module_l[, l]; gl <- gl[gl != '']
  gl
})
names(tumor_module_ls) <- names(gm)
lapply(tumor_module_ls, length)



##### most module only have a few genes. But need to check the AGER (25 genes)
intersect(str_split(at1_mark, '_', simplify = T)[, 2], tumor_module_ls[['AGER']]) ## 4 / 25 overlap. these 4 genes are within the top 7 genes in the module list. 

