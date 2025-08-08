library(ggrepel)
library(ggplot2)
library(dplyr)
library(data.table)
setwd(here::here("./Scripts/Figure"))

ip <- '../..//Mixed_model_output/'


nsclc_lm <- fread(file.path(ip,'model_summary_Histology.csv')) 
rownames(nsclc_lm) <- nsclc_lm$V1; nsclc_lm$V1 <- NULL
nsclc_lm$study_VarPer <- nsclc_lm$study_VarPer * 100
nsclc_lm$Fixed_VarPer <- nsclc_lm$Fixed_VarPer * 100

### matching the module annotation as label5
module_anno <- fread('../..//Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))


#### 1.1 Study VS Fixed
nsclc_lm$label <- nsclc_lm$module
nsclc_lm$label <- module_anno[match(nsclc_lm$label, module_anno$Module), ] %>% pull(Label)
#nsclc_lm$label[!nsclc_lm$module %in% paste0('L', 1:98)] <- nsclc_lm[!nsclc_lm$module %in% paste0('L', 1:98), ] %>% pull(module)
  

fn <- "../../Figure/Figure2/Figure2_Variability_scatter.pdf"
pdf(fn, width = 2.5, height = 2.5)
plot(
  ggplot(aes(x = study_VarPer, y = Fixed_VarPer, label = label), data = nsclc_lm) +
    geom_point(color = dplyr::case_when(
      nsclc_lm$study_VarPer > 20 ~ "black", 
      nsclc_lm$Fixed_VarPer >= 20 ~ "red4",
      TRUE ~ "grey"), alpha = 0.8, size = 0.8) + 
    ggrepel::geom_label_repel(data = dplyr::filter(nsclc_lm, study_VarPer >= 30 | Fixed_VarPer >= 30),
                              box.padding   = 0.5, 
                              point.padding = 0.5,
                              size = 1, xlim = c(0.2, NA),
                              max.overlaps = 16,
                              segment.color = 'grey50') +
    theme_classic() + 
    theme(axis.text = element_text(size = 4),
          axis.title = element_text(size = 6)
    )  +  
    labs(x= "% Study variability", y = "% Stage+Histology variability")
)
dev.off()



#### 2 Volcano plot

##### 2.1 Histology associated
nsclc_lm <- fread(file.path(ip,'model_summary_Histology.csv')) 
rownames(nsclc_lm) <- nsclc_lm$V1; nsclc_lm$V1 <- NULL
nsclc_lm$Histology_fdr <- p.adjust(nsclc_lm$Histology_pval, 'fdr')
nsclc_lm$label <- nsclc_lm$module
nsclc_lm$label[!nsclc_lm$module %in% paste0('L', 1:98)] <- ''

nsclc_lm$label <- module_anno[match(nsclc_lm$label, module_anno$Module), ] %>% pull(Label)


ggplot(aes(x = `Histology_coef`, y = -log(Histology_pval, 10), label = label), data = nsclc_lm) +
  geom_hline(yintercept = -log(0.05, 10)) +
  geom_point(color = dplyr::case_when(
    nsclc_lm$`Histology_coef` >= log2(1.5) ~ "red4", 
    nsclc_lm$`Histology_coef` <= -log2(1.5) ~ "royalblue",
    nsclc_lm$Histology_fdr < 0.05 ~ 'black',
    TRUE ~ "grey"), alpha = 0.8) + xlim(-5,5) + 
  ggrepel::geom_label_repel(data = dplyr::filter(nsclc_lm, `Histology_coef` > log2(1.5) ), #& module %in% paste0('L', 1:67)
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 4, xlim = c(1, NA),
                            segment.color = 'grey50') +
  ggrepel::geom_label_repel(data = dplyr::filter(nsclc_lm, `Histology_coef` < -1*log2(1.5) ), #& module %in% paste0('L', 1:67)
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 4, xlim = c(NA, -1),
                            segment.color = 'grey50') +
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 20)
  ) + ggtitle(label = 'LUSC_VS_LUAD') + 
  labs(x= "log2(FC)", y = "-log10(fdr)")




##### 2.2 Stage associated
luad_lm <- fread(file.path(ip, 'model_summary_Stage.csv')) ### single-cell level checking, luad_singleCell_randomEffect_StudyOnly.csv
rownames(luad_lm) <- luad_lm$V1; luad_lm$V1 <- NULL
luad_lm$Stage_fdr <- p.adjust(luad_lm$Stage_pval, 'fdr')
luad_lm$label <- luad_lm$module
luad_lm$label[!luad_lm$module %in% paste0('L', 1:98)] <- ''

luad_lm$label <- module_anno[match(luad_lm$label, module_anno$Module), ] %>% pull(Label)



fn <- file.path('../../Figure/Figure3', 'Figure3_Stage_Volcano.pdf')
pdf(file = fn, width = 11, height = 7)
ggplot(aes(x = `Stage_coef`, y = -log(Stage_fdr, 10), label = label), data = luad_lm) +
  geom_hline(yintercept = -log(0.05, 10)) +
  geom_point(color = dplyr::case_when(
    luad_lm$`Stage_coef` >= log2(1.5) ~ "red4", 
    luad_lm$`Stage_coef` <= -log2(1.5) ~ "royalblue",
    luad_lm$Stage_fdr < 0.05 ~ 'black',
    TRUE ~ "grey"), alpha = 0.8) + xlim(-4,4) + 
  ggrepel::geom_label_repel(data = dplyr::filter(luad_lm, `Stage_coef` > log2(1.5) & Stage_fdr < 0.05), #& module %in% paste0('L', 1:67) | nsclc_lm$module %in% c('L35', 'L36')
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 4, xlim = c(1, NA), 
                            segment.color = 'grey50', seed = 1) +
  ggrepel::geom_label_repel(data = dplyr::filter(luad_lm, `Stage_coef` < -1*log2(1.5) & Stage_fdr < 0.05), #& module %in% paste0('L', 1:67) 
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 4, xlim = c(NA, -1), 
                            segment.color = 'grey50', seed = 1) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20)
  ) + ggtitle(label = 'Late_VS_Early') + 
  labs(x= "log2(FC)", y = "-log10(fdr)")
dev.off()
