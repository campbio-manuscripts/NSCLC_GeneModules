### this script is used to look at healthy tissue
library(celda)
library(singleCellTK)
library(ggplot2)
library(dplyr)
library(eulerr)
library(data.table)
library(stringr)

### 1. load data
# fp <- '/restricted/projectnb/camplab/projects/Single_Cell_Public/Basil_Nature_2022_Lung_Terminal/Analysis/Celda/2023.01.22'
# nor_sce <- readRDS(file.path(fp, 'celda_cg_decontXcounts_min3_filter_mito10_detected500_sum750_sub_K40_L125.rds'))
fp <- "../../Data"
nor_sce <- readRDS(file.path(fp, 'Basil_Nature_2022_celda_cg_decontXcounts_min3_filter_mito10_detected500_sum750_sub_K40_L125.rds'))
nor_sce <- runNormalization(nor_sce, outAssayName = 'logNormalized', normalizationMethod = 'LogNormalize')
nor_sce$celda_k <- colData(altExp(nor_sce))$celda_cell_cluster


## 2. Get the list of module gene
nor_ta <- featureModuleTable(nor_sce)

gt <- c('SFTPC', 'SFTPD','SFTPB', 'SFTPA1', 'SFTPA2', 'SCGB3A1', 'SCGB3A2', "AGER", "AGR3", "HOPX")
gs <- rownames(nor_sce)[match(gt, rowData(nor_sce)$Symbol)]
names(gs) <- gt
featureModuleLookup(nor_sce, features = gs)


lineage_marker <- c('NKX2-1', 'RNASE1', 'NAPSA', 
                    'SFTPC', 'SFTPB', 'SFTPD', 'SFTPA1', 'SFTPA2', "SCGB3A2", "SCGB3A1",
                    'HOPX', 'AGER', 'AGR3')
lineage_gmar <- rowData(nor_sce) %>% as.data.frame() %>% filter(Symbol %in% lineage_marker)




## 3. Include normal AT1, AT2 cells together. Check the expression of our nsclc_meta AT1 and AT2 modules / markers
nor_at2at1 <- nor_sce[, nor_sce$celda_k %in% as.character(c(11:20, 22, 25, 30, 31, 23, 24))]
nor_at2at1$celda_k <- as.character(nor_at2at1$celda_k)
nor_at2at1$cell_type <- 'AT2'
nor_at2at1$cell_type[nor_at2at1$celda_k == 25] <- 'SCGB3A1+/SCGB1A1-'
nor_at2at1$cell_type[nor_at2at1$celda_k %in% c(23, 24)] <- 'AT1'

#### probability based
nx <- factorizeMatrix(nor_at2at1, type = "proportion")
cell_probs <- nx[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs

lineage_mod <- featureModuleLookup(nor_at2at1, features = rownames(lineage_gmar)) ### look up the module list
names(lineage_mod) <- lineage_gmar$Symbol

# cf <- 0.01
cf <- 0.005
lm_on <- 1* t(module_prob[, unique(lineage_mod)] > cf) #target_module
lm_on <- t(lm_on) %>% as.data.frame()



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
m <- c("SFTPC", "SFTPA1/2",  "SFTPB", "SFTPD", "AGER", "SCGB3A2")
md <- m_x[m]
fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,
fill_col <- setNames(RColorBrewer::brewer.pal(length(m), "Set3"), m)


pdf("../../Figure/Figure4/Figure4_Lung_lineageVenn_probbased_cf0.005.pdf", width = 4, height = 4)
plot(fit,
     quantities = list(type = c("percent"), cex = 0.5), #TRUE #, "counts", "percent"
     fill = fill_col,
     lty = 1:length(m),
     labels = list(cex = 0.5),
     #labels = NULL,
     legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"), 
                   labels = rownames(fit$ellipses), 
                   cex = 0.5))
dev.off()


miss_op <- data.frame('OP' = names(fit$original.values), 
                      'Number' = fit$original.values)
miss_op$OP_Per <- miss_op$Number / sum(miss_op$Number) * 100
missOP <- miss_op %>% filter(OP_Per != 0) %>% arrange(OP_Per)
# sum(missOP$NumOP); sum(missOP$OP_Per)

write.csv(miss_op, 
          "../..//Table/Supplementary_Normal_Lung_VennDiagram.csv")


