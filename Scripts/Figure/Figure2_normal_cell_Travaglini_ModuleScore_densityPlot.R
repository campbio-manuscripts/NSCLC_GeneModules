### this script is used to look at healthy tissue
library(celda)
library(singleCellTK)
library(ggplot2)
library(dplyr)
library(eulerr)
library(stringr)
setwd(here::here("./Scripts/Figure"))

### 1. load the sce object: /restricted/projectnb/camplab/projects/Single_Cell_Public/Travaglini_Nature_2020_Lung_Atlas/Analysis/Celda/2023.01.15
nor_sce <- readRDS("../../Data/Travaglini_Nature_2020_celda_cg_counts_filter_detected500_sum750_sub_K30_L100.rds")
nor_sce <- runNormalization(nor_sce, outAssayName = 'logNormalized', normalizationMethod = 'LogNormalize')

### 2. plot UMAP with celda cluster
nor_sce$celda_k <- colData(altExp(nor_sce))$celda_cell_cluster
reducedDim(nor_sce, 'celda_umap') <- reducedDim(altExp(nor_sce), 'celda_UMAP')


g <- 'SFTPC'
g <- 'AGER'
g <- 'MUC5AC'
g <- 'MUC5B'
g <- 'SCGB3A2'
g <- 'KRT5'
g <- 'FOXJ1'
g <- 'SCGB1A1'
g <- 'SCGB3A1'
g <- 'HOPX'

col_scale <- c(0, 0.3, 0.6, 0.8, 0.995)
colors <- c('white', 'snow','orange','red', 'red3')


plotSCEDimReduceFeatures(nor_sce, feature=g, reducedDimName = 'celda_umap', useAssay = 'logNormalized',
                         legendTitle = 'Expression',
                         dotSize = 0.1, title = g, titleSize = 28, legendSize = 10, legendTitleSize = 15) + 
  scale_color_gradientn(values = col_scale, 
                        colours = colors,
                        na.value = "grey50") + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(5, 'mm'), legend.margin=margin(0,0,0,0)
  )


plotSCEDimReduceColData(nor_sce, 
                        colorBy = 'sample', reducedDimName = 'celda_umap', 
                        clusterLabelSize = 4.5, dotSize = 0.1, labelClusters = T, 
                        colorScale = NULL, title = "celda_k", titleSize = 8) + 
  theme(
    axis.title = element_blank(),         # Change both x and y axis titles
    axis.text=element_blank(),
    axis.ticks = element_blank()
    #legend.position="none"
  )



### 3. generate the score
y <- factorizeMatrix(nor_sce, type = "counts")
cell_m_counts <- y[[1]][['cell']]; cell_m_counts <- t(cell_m_counts)
colSum_sizeFactor <- colSums(t(cell_m_counts)) / 10000 
normalized_m_counts <- sweep(cell_m_counts, 1, colSum_sizeFactor, "/")
module_activity <- log2(normalized_m_counts + 1)



## 4. Get the list of module gene
nor_ta <- featureModuleTable(nor_sce)

gt <- c('SFTPC', 'SFTPD','SFTPB', 'SFTPA1', 'SFTPA2', 'SCGB3A1', 'SCGB3A2', "AGER", "AGR3") ###"HOPX" module is too big, overally expressed
gs_mod <- featureModuleLookup(nor_sce, features = gt)
target_modues <- paste0('L', gs_mod) %>% unique()


### 5 . use the josh script to find turning points on the density plot
mod <- t(module_activity)


pdf("../../Figure/Figure4/Travaglini_Lung_Normal_density_test_bw1.pdf")
density.res <- list()
num.mins <- 5
bw = 1 #0.2
minDiff <- 0.03

for(i in target_modues) {
  x <- mod[i,]
  
  adjust <- 1
  mins <- rep(NA, num.mins + 1)
  while(length(mins) >= num.mins) {
    den <- density(x, adjust=adjust, n=2^10)
    dd <- diff(den$y)
    
    mins <- which(dd[-1] >= 0 & dd[-length(dd)] < 0) + 1
    maxs <- which(dd[-1] <= 0 & dd[-length(dd)] > 0) + 1
    
    adjust <- adjust + bw # 1
  }  
  
  nd <- den$y / max(den$y)    
  diffs1 <- nd[maxs[-1]] - nd[mins]
  diffs2 <- nd[maxs[-length(maxs)]] - nd[mins]
  
  mins.final <- mins[which(diffs2 > minDiff)]
  maxs.final <- maxs[which(diffs2 > minDiff)]
  
  breaks <- cut(x, breaks=c(min(x), den$x[mins.final], max(x)), include.lowest=TRUE, right=TRUE)
  
  hist(x, breaks=100, freq=FALSE, main = paste(i, adjust-1, sep = ": "))
  lines(den$x, den$y, col="black", lwd=2)
  points(den$x[maxs.final], den$y[maxs.final], pch=16, col="red")
  points(den$x[mins.final], den$y[mins.final], pch=16, col="blue")
  # density.res[[i]] <- list(min=mins, max=maxs, breaks=breaks)
  density.res[[i]] <- list(max = den$y[maxs.final], min = den$x[mins.final])
  print(i)
}
dev.off()




### 6. combine all modules together. Check the cutoff point
x <- as.numeric(unlist(mod))
adjust <- 1
mins <- rep(NA, num.mins + 1)
while(length(mins) >= num.mins) {
  den <- density(x, adjust=adjust, n=2^10)
  dd <- diff(den$y)
  
  mins <- which(dd[-1] >= 0 & dd[-length(dd)] < 0) + 1
  maxs <- which(dd[-1] <= 0 & dd[-length(dd)] > 0) + 1
  
  adjust <- adjust + bw # 1
}  

nd <- den$y / max(den$y)    
diffs1 <- nd[maxs[-1]] - nd[mins]
diffs2 <- nd[maxs[-length(maxs)]] - nd[mins]

mins.final <- mins[which(diffs2 > minDiff)]
maxs.final <- maxs[which(diffs2 > minDiff)]

breaks <- cut(x, breaks=c(min(x), den$x[mins.final], max(x)), include.lowest=TRUE, right=TRUE)

hist(x, breaks=100, freq=FALSE, main = paste(i, adjust-1, sep = ": "))
lines(den$x, den$y, col="black", lwd=2)
points(den$x[maxs.final], den$y[maxs.final], pch=16, col="red")
points(den$x[mins.final], den$y[mins.final], pch=16, col="blue")
list(max = den$x[maxs.final], min = den$x[mins.final])


### 7. define on and off
lm_on <- module_activity[, ]

mod_on_status <- lapply(names(density.res), function(i) {
  print(i)
  cutoff <- den$x[mins.final][1] ### unified cutoff
  print(cutoff)
  1*(rowSums(lm_on[, paste0('L', 53), drop=F]) > cutoff)
}) %>% do.call('cbind', .)
colnames(mod_on_status) <- names(density.res)
mod_on_status <- as.data.frame(mod_on_status)


m_x <- lapply(colnames(mod_on_status), function(l) {
  sd <- mod_on_status[mod_on_status[[l]] == 1, ]
  rownames(sd)
})
names(m_x) <- names(gs_mod)[match(str_sub(colnames(mod_on_status), start = 2L), gs_mod)]

library(eulerr)
m <- names(m_x)

md <- m_x[m]
set.seed(2020)
fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,

p <- plot(fit,
          quantities = list(type = c("counts"), cex = 1), #TRUE #, "counts", "percent"
          #quantities = T,
          fill = RColorBrewer::brewer.pal(length(m), "Set3"),
          lty = 1:length(m),
          #labels = list(cex = 3),
          #adjust_labels = T,
          #labels = NULL,
          legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"),
                        labels = rownames(fit$ellipses),
                        cex = 1))
p
