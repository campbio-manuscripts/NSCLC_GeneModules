#### analyze ALI dataset
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

source('../../Scripts/Figure/Figure_utils.R')


### 1. load dataset
# ip <- '/restricted/projectnb/camplab/projects/Single_Cell_Public/Wohnhaas_IJMS_2022_Airway_Smoke/Analysis/2023.01.12'
# ali_sce <- readRDS(file.path(ip, 'celda_cg_counts_min5_filter_mito10_detected500_sum750_sub_K40_L125.rds'))

ip <- "../../Data/"
ali_sce <- readRDS(file.path(ip, 'Wohnhaas_IJMS_2022_Airway_Smoke_K40_L125.rds'))
reducedDim(ali_sce, 'celda_UMAP') <- reducedDim(altExp(ali_sce), 'celda_UMAP')

# l <- 'Seurat_louvain_Resolution0.8'
# l <- 'disease'
# l <- 'duration'
# l <- 'exposure'
# plotSCEDimReduceColData(ali_sce, colorBy = l, 
#                         reducedDimName = 'celda_UMAP')



### 2. subset the data
#### only keep healthy donor
hl <- ali_sce$disease == 'healthy'
ali_healthy <- ali_sce[, hl]
ali_healthy <- runNormalization(ali_healthy, outAssayName = 'logNormalized', normalizationMethod = 'LogNormalize')

# g <- 'KRT5'
# g <- 'KRT6A'
# g <- 'TP63'
# g <- 'NFE2L2'
# 
# gid <- rowData(ali_healthy)[match(g, rowData(ali_healthy)$feature_name), 'feature_ID']
# plotSCEDimReduceFeatures(ali_healthy, feature = gid, useAssay = 'logNormalized', title = g,
#                          reducedDimName = 'celda_UMAP')



ali_healthy$Seurat_louvain_Resolution0.8 <- as.character(ali_healthy$Seurat_louvain_Resolution0.8)


# l <- 'Seurat_louvain_Resolution0.8'
# plotSCEDimReduceColData(ali_healthy, colorBy = l, 
#                         reducedDimName = 'celda_UMAP')
# 
# g <- 'KRT5'
# gid <- rowData(ali_healthy)[match(g, rowData(ali_healthy)$feature_name), 'feature_ID']
# 
# plotSCEViolinAssayData(ali_healthy, feature = gid, useAssay = 'logNormalized', 
#                        groupBy = 'Seurat_louvain_Resolution0.8')


#### only keep basal cells
cid <- c(0, 19, 14, 22, 29, 8, 11, 31) #23
ali_healthy_basal <- ali_healthy[, ali_healthy$Seurat_louvain_Resolution0.8 %in% cid]
ali_healthy_basal$Seurat_louvain_Resolution0.8 <- as.character(ali_healthy_basal$Seurat_louvain_Resolution0.8)


plotSCEDimReduceFeatures(ali_healthy_basal, feature = gid, useAssay = 'logNormalized', title = g,
                         reducedDimName = 'celda_UMAP')

l <- 'exposure'
l <- 'duration'
plotSCEDimReduceColData(ali_healthy_basal, colorBy = l, 
                        reducedDimName = 'celda_UMAP')


### 3. Quantify the expression of our NSCLC-meta modules
nsclc_module <- fread('../../Data/Final_module_list.txt') %>% as.data.frame()

ali_healthy_basal_m <- lapply(colnames(nsclc_module), function(l) {
  print(l)
  g <- unlist(nsclc_module[, l]); g <- g[g != '']
  gid <- rowData(ali_healthy)[match(g, rowData(ali_healthy)$feature_name), 'feature_ID']
  gid <- na.omit(gid)
  colSums(assay(ali_healthy_basal)[gid, , drop=F])
  
})
ali_healthy_basal_m <- do.call('rbind', ali_healthy_basal_m)
rownames(ali_healthy_basal_m) <- colnames(nsclc_module)


colSum_sizeFactor <- colSums((assay(ali_healthy_basal))) / 10000 
ali_he_bal_norm <- sweep(ali_healthy_basal_m, 2, colSum_sizeFactor, "/")
ali_he_bal_norm_log <- log2(ali_he_bal_norm + 1)




### 4. Identify differentially expressed modules between air and smoking, taking exposure time into consideration. 
basal_meta <- colData(ali_healthy_basal)
basal_meta$cluster <- basal_meta$Seurat_louvain_Resolution0.8

anova_lm <- list()
for (l in rownames(ali_healthy_basal_m)) {
  print(l)
  d <- data.frame('module' = ali_he_bal_norm_log[l, rownames(basal_meta)],
                  basal_meta[, c('subject', 'cluster', 'duration', 'exposure')])
  
  m <- lme4::lmer(module ~ exposure + duration + (1|subject), data=d) #cluster
  
  ### Get summary statistics
  m_sum <- summary(m)$coefficients
  af <- car::Anova(m, type = 'III')
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  anova_lm[[l]] <- data.frame('module' = l,  'subject_VarPer' = as.numeric(vara$var.intercept['subject']/all_var),  #'celda_VarPer' = as.numeric(vara$var.intercept['cluster:subject']/all_var),
                              'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                              'Residual_VarPer' = vara[['var.residual']]/all_var,
                              'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                              'Exposure_pval' = af['exposure', 3],
                              'Smoke_coef' =  m_sum['exposuresmoke', 1],
                              'Duration_pval' = af['duration', 3],
                              'Chronic_coef' = m_sum['durationchronic', 1])
  
}
anova_lm_s <- do.call('rbind', anova_lm)
anova_lm_s$smoke_fdr <- p.adjust(anova_lm_s$Exposure_pval, 'fdr')
anova_lm_s$chronic_fdr <- p.adjust(anova_lm_s$Duration_pval, 'fdr')
write.csv(anova_lm_s, 
          '../../Scripts/Linear_Mixed_Model/ali_healthy_exposure_module_summary.csv')

anova_lm_s %>% filter(smoke_fdr < 0.05, Smoke_coef > 0.5) %>% View()




### 5. Plot heatmap of DE modules
source('../../Scripts/Figure/Figure_utils.R')

mod_anno <- fread('../../Data/Module_annotation_Reorder_01202025.csv')
mod_anno$Label <- paste(mod_anno$Module, mod_anno$Name, sep = ':')


### adding some interesting modules in ALI study
target_modules <- anova_lm_s %>% filter(smoke_fdr < 0.05, abs(Smoke_coef) > 0.2) %>% pull(module) #0.5
# target_modules <- c(target_modules, 
#                     c('L52', 'L47', 'L49'))

# hm_expr <- t(scale(t(ali_he_bal_norm_log[target_modules, ])))
# rownames(hm_expr) <- mod_anno[match(rownames(hm_expr), mod_anno$Module), 'Label'] %>% unlist()
# 
# 
# 
# 
# colSplit <- basal_meta[colnames(ali_he_bal_norm_log), 'phenotype'] #exposure
# colSplit <- factor(colSplit, levels = c('smoke_chronic', 'air_chronic', 'smoke_acute', 'air_acute'))
#   
# rowSplit <- ifelse(anova_lm_s[target_modules, 'Smoke_coef'] > 0, 'Up_Smoke', 'Up_Air')
# 
# colAnno <- HeatmapAnnotation(
#   df = basal_meta[colnames(ali_he_bal_norm_log), c('duration', 'subject')]
# )
# 
# 
# ### plot
# Heatmap(hm_expr, column_split = colSplit, cluster_columns = F, row_split = rowSplit, 
#         cluster_rows = T, top_annotation = colAnno, 
#         show_column_names = F, show_column_dend = F, show_row_dend = F)














### 6.summarize the z score within each phenotype

##### need to generate z-score within each subject. 
subjects <-  as.character(unique(basal_meta[['subject']]))
z_subject <- lapply(subjects, function(s) {
  sz <- t(scale(t(ali_he_bal_norm_log[target_modules, basal_meta$subject == s])))
  sz
})
subject_z_module <- do.call('cbind', z_subject)

##### summarize by sample
# fea <- 'sample'
# features <- as.character(unique(basal_meta[[fea]]))

###### also try to average all three subjects together
fea <- 'phenotype'
features <- as.character(unique(basal_meta[[fea]]))

  
  
fea_summary <- lapply(features, function(f) {
  cid <- rownames(basal_meta)[basal_meta[[fea]] == f] ##& data_sce$Histology %in% his
  
  
  d <- data.frame('feature' = f, 
                  'avg.expr' = rowMeans(as.matrix(subject_z_module[, cid])))
  #d$module <- rownames(d)
  d$module <- rownames(d)
  d
})
fea_summary <- do.call('rbind', fea_summary)

module_expr <- fea_summary %>%
  pivot_wider(names_from = feature, values_from = avg.expr) %>% as.data.frame()
rownames(module_expr) <- module_expr$module
module_expr <- module_expr[,-1] %>% as.matrix()


# colSplit <- basal_meta[match(colnames(new_hm_expr), basal_meta$sample), 'subject']
# rowSplit <- ifelse(anova_lm_s[target_modules, 'Smoke_coef'] > 0, 'Up_Smoke', 'Up_Air')
# 
# colAnno <- HeatmapAnnotation(
#   df = basal_meta[match(colnames(new_hm_expr), basal_meta[[fea]]), c('exposure','duration', 'subject')]
# )


### average of three subjects
new_hm_expr <- module_expr

lev <- c('air_acute', 'air_chronic', 'smoke_acute', 'smoke_chronic')
new_hm_expr <- new_hm_expr[, lev]

colSplit <- NULL
rowSplit <- ifelse(anova_lm_s[target_modules, 'Smoke_coef'] > 0, 'Up_Smoke', 'Up_Air')

# colAnno <- HeatmapAnnotation(
#   df = basal_meta[match(colnames(new_hm_expr), basal_meta[['phenotype']]), c('exposure','duration')], 
#   annotation_name_gp= gpar(fontsize = 5)
# )


colAnno <- HeatmapAnnotation(
  Duration = basal_meta[match(colnames(new_hm_expr), basal_meta[['phenotype']]), c('duration')],
  Exposure = basal_meta[match(colnames(new_hm_expr), basal_meta[['phenotype']]), c('exposure')],
  annotation_name_gp= gpar(fontsize = 6), ###font size of name on the top annotation bar
  col = list(Duration = c("acute" = "red3", "chronic" = "steelblue"),
             Exposure = c("air" = "orange", "smoke" = "limegreen")
  ), 
  annotation_legend_param = list(labels_gp = gpar(fontsize = 5), 
                                 grid_height = unit(3, "mm"),
                                 grid_width = unit(3, "mm")), 
  simple_anno_size  = unit(0.3, "cm")
  
)


rownames(new_hm_expr) <- mod_anno[match(rownames(new_hm_expr), mod_anno$Module), 'Label'] %>% unlist()

hm <- Heatmap(new_hm_expr, 
              name = "Module\nScore", 
              row_split = rowSplit, column_split = colSplit, 
              cluster_rows = F, cluster_columns = F, top_annotation = colAnno, 
              row_title_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              height = nrow(new_hm_expr)*unit(2.7, "mm"),
              width = ncol(new_hm_expr)*unit(3.7, "mm"),
              column_names_rot = 45,
              heatmap_legend_param = list(
                # title_position = "lefttop",
                # legend_direction = "horizontal", 
                title_position = "topleft",
                legend_direction = "vertical", 
                labels_gp = gpar(fontsize = 5),
                grid_width  = unit(0.3, "cm")))


pdf("../../Figure/Figure6/Figure7_Ali_module_heatmap_FC0.2.pdf", width = 3, height = 7)
draw(hm, 
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(0.1, 1, 0.1, 3), "mm"), 
     legend_title_gp = gpar(fontsize = 8, fontface = "bold") ### legend title font size
)

dev.off()




##### Also summarize by sample
fea <- 'sample'
features <- as.character(unique(basal_meta[[fea]]))

fea_summary <- lapply(features, function(f) {
  cid <- rownames(basal_meta)[basal_meta[[fea]] == f] ##& data_sce$Histology %in% his
  
  
  d <- data.frame('feature' = f, 
                  'avg.expr' = rowMeans(as.matrix(subject_z_module[, cid])))
  #d$module <- rownames(d)
  d$module <- rownames(d)
  d
})
fea_summary <- do.call('rbind', fea_summary)

module_expr <- fea_summary %>%
  pivot_wider(names_from = feature, values_from = avg.expr) %>% as.data.frame()
rownames(module_expr) <- module_expr$module
module_expr <- module_expr[,-1] %>% as.matrix()


new_hm_expr <- module_expr


colSplit <- basal_meta[match(colnames(new_hm_expr), basal_meta$sample), 'subject']
rowSplit <- ifelse(anova_lm_s[target_modules, 'Smoke_coef'] > 0, 'Up_Smoke', 'Up_Air')

colAnno <- HeatmapAnnotation(
  df = basal_meta[match(colnames(new_hm_expr), basal_meta[[fea]]), c('exposure','duration', 'subject')]
)


rownames(new_hm_expr) <- mod_anno[match(rownames(new_hm_expr), mod_anno$Module), 'Label'] %>% unlist()

hm <- Heatmap(new_hm_expr, 
              name = "Module\nScore", 
              row_split = rowSplit, column_split = colSplit, 
              cluster_rows = F, cluster_columns = F, top_annotation = colAnno, 
              row_title_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              height = nrow(new_hm_expr)*unit(2.7, "mm"),
              width = ncol(new_hm_expr)*unit(3.7, "mm"),
              column_names_rot = 45,
              heatmap_legend_param = list(
                # title_position = "lefttop",
                # legend_direction = "horizontal", 
                title_position = "topleft",
                legend_direction = "vertical", 
                labels_gp = gpar(fontsize = 5),
                grid_width  = unit(0.3, "cm")))

pdf("../../Figure/Figure6/Figure7_Ali_module_heatmap_subject_level_FC0.2.pdf", width = 5, height = 7)
draw(hm, 
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(0.1, 1, 0.1, 3), "mm"), 
     legend_title_gp = gpar(fontsize = 5, fontface = "bold") ### legend title font size
)

dev.off()


### 7.Visualize some marker genes across these samples. 
markers <- c('SOX2', 'TP63', 'KRT6A', 'KRT13', 'CEACAM5', 'CEACAM6', 'SCGB1A1', 'MUC5AC')

marker_expr <- assay(ali_healthy_basal, 'logNormalized')[rowData(ali_healthy_basal)[match(markers,rowData(ali_healthy_basal)$feature_name), 'feature_ID'], ]
rownames(marker_expr) <- markers


##### need to generate z-score within each subject. 
subjects <-  as.character(unique(basal_meta[['subject']]))
z_subject <- lapply(subjects, function(s) {
  sz <- t(scale(t(marker_expr[, basal_meta$subject == s])))
  sz
})
subject_z_marker <- do.call('cbind', z_subject)



##### summarize by sample
fea <- 'sample'
features <- as.character(unique(basal_meta[[fea]]))

###### also try to average all three subjects together
# fea <- 'phenotype'
# features <- as.character(unique(basal_meta[[fea]]))



fea_summary <- lapply(features, function(f) {
  cid <- rownames(basal_meta)[basal_meta[[fea]] == f] ##& data_sce$Histology %in% his
  
  
  d <- data.frame('feature' = f, 
                  'avg.expr' = rowMeans(as.matrix(subject_z_marker[, cid])))
  #d$module <- rownames(d)
  d$module <- rownames(d)
  d
})
fea_summary <- do.call('rbind', fea_summary)

marker_expr <- fea_summary %>%
  pivot_wider(names_from = feature, values_from = avg.expr) %>% as.data.frame()
rownames(marker_expr) <- marker_expr$module
marker_expr <- marker_expr[,-1] %>% as.matrix()

new_hm_expr <- marker_expr

lev <- c("acute_A_air", "chronic_A_air", "acute_A_smoke",  "chronic_A_smoke",
         "acute_B_air", "chronic_B_air", "acute_B_smoke",  "chronic_B_smoke",
         "acute_C_air", "chronic_C_air", "acute_C_smoke",  "chronic_C_smoke")
new_hm_expr <- new_hm_expr[, lev]

colSplit <- basal_meta[match(colnames(new_hm_expr), basal_meta$sample), 'subject']
rowSplit <- NULL

colAnno <- HeatmapAnnotation(
  df = basal_meta[match(colnames(new_hm_expr), basal_meta[[fea]]), c('exposure','duration', 'subject')]
)

Heatmap(new_hm_expr, row_split = rowSplit, column_split = colSplit, 
        cluster_rows = F, cluster_columns = F, top_annotation = colAnno)



#### weird, look the marker expression using violin plot at single cell level
for (g in markers) {
  gid <- rowData(ali_healthy_basal)[match(g, rowData(ali_healthy_basal)$feature_name), 'feature_ID']
  print(
    plotSCEViolinAssayData(ali_healthy_basal, feature = gid, useAssay = 'logNormalized', title = g,
                           groupBy = "phenotype")
  )
}



# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8305830/

### study design
#### disease: healthy vs COPD
#### smoking: smoking vs air
#### smoking duration: acute versus chronic


#### Acute response
# Iron homeostasis and oxidative stress associated ferritin subunits were induced in all cell types [49]
# increased expression of MUC1 and MUC5AC, although the latter was only significantly induced in COPD donors, and a down-regulation of SC4 associated SFTPB
#  Expression of the early barrier alarmins KRT16, KRT17 and KRT6 by keratinocytes is associated with barrier breach and hyperproliferation, suggesting that the expression of these markers in the BC3 cluster is associated with impaired barrier integrity of the small airway epithelium


#### chronic
# Moreover, KRT6A+ cells were significantly increased within the basal cell population after chronic smoke exposure compared to acute smoking but not the respective air controls
# As observed for acute smoke exposure, shifts of cell subpopulation frequencies were overall similar across HC and COPD donors, suggesting a similar response across both groups 
#  Secretory cells were more strongly affected by chronic than acute smoke exposure when compared to air controls
# Chronic smoke exposure induced the expression of oxidative phosphorylation (e.g., COX5A, NDUFB2, ATP5J, UQCRQ) and translation associated gene modules, including mitochondrial ribosome (MRPS15, MRPL20, MRPL51) and endoplasmic reticulum processing and transport (SEC61B, SEC11C, SRP14, SSR4) associated gene
#  specific or increased induction of immunity related molecular processes, including antigen processing and presentation (HLA-DRA, HLA, DRB1, HLA, HLA-B) as well as antimicrobial and immune cell recruitment associated genes (LTF, S100A8, CXCL6, CXCL1, SAA1, SAA2) 
# several keratins associated with epithelial cell and squamous differentiation (KRT6A, KRT13, KRT14) and SPRR1B [3,35,36,52] were selectively induced after chronic smoke exposure, indicating the onset of basal cell hyperplasia and squamous





### report
### # no ciliated cells. 
### basal on the left, secretory cells on the right. 
### 3 healthy, 3 copd. Maybe remove COPD. 
### cluster on the bottom right: chronic and smooke samples, versus on the top right


### small cluster with S100A8/SPRR cells might mapped to our keratinization cells in LUSC. and they are from smokers. Metapiasia. 
### smoking won't separate well from air. SPRR1B more in basal cells. 
### VCAN EMT high in basal of one patient E. 
### gsva to project score. Or quick: counts / # total counts. 
### More interested: KRT13/ KRT6 associaetd with smoking. Cycokine signaling, or other associated with smoking. 
### big picture: smoking causing mutation, mutation cause cancer. A lot of passenger mutation. 
### smoking change the transcriptome changes. Also expresseed in cancer cells. Somehow get maintained in tumor. 
### smoking shift cell state. Even the patient quit smoking. 


### small airway basal: terminal airway. Not sure. But at least they are all KRT5. 
### rownames of the heatmap. Variability figures. 


