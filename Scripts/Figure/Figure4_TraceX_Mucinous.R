### analyze association between module score and stage or histology subtype
library(GSVA)
library(ASSIGN)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(SummarizedExperiment)
library(dplyr)
#library(CaDrA)
library(circlize)
library(tidyr)
library(stringr)
library(SingleCellExperiment)
library(celda)
library(fst)

setwd(here::here("./Scripts/Figure"))


### 1. load data
tp <- "../../Data"

# tp_meta <- read_fst(file.path(tp, '2022-10-18all_metadata.fst'))
tp_meta <- read_fst(file.path(tp, 'Tracex_metadata.fst'))

tp_module <- fread("../../Tracex_output/Tracex_Module_Score.csv") %>% as.data.frame()
rownames(tp_module) <- tp_module$V1; tp_module$V1 <- NULL

### 2. subset to the sample that has annotation in tp_meta as LUAD and LUSC
nsclc_meta <- tp_meta %>% filter(Histology %in% c("LUAD", 'LUSC', 'Collision LUAD and LUSC', 'Adenosquamous carcinoma'))
tp_id <- colnames(tp_module)[colnames(tp_module) %in% nsclc_meta$region] ###541 samples in total
nsclc_meta <- nsclc_meta %>% filter(region %in% tp_id)



### 3. process the meta data
nsclc_meta <- nsclc_meta %>% mutate(
  Stage = case_when(
    stage %in% c("IA", "IB") ~ "I",
    stage %in% c("IIA", "IIB") ~ "II",
    stage %in% c("IIIA", "IIIB") ~ "III" 
  )
)



### 4. Perform association with luad histological subtype
luad_meta <- nsclc_meta %>% filter(Histology == "LUAD", 
                                   !luad_subtype %in% c("not possible (too crushed)", "cribriform")) ### mucious later, "mucinous"
luad_meta <- luad_meta %>% filter(!is.na(luad_subtype))
luad_module <- tp_module[, luad_meta$region]



mucinous_lm <- lapply(rownames(luad_module), function(l) {
  d <- luad_meta %>% filter(!is.na(luad_subtype))
  d <- luad_meta %>% filter(luad_subtype %in% c('acinar', 'lepidic', 'micropapillary', 'mucinous', 'papillary')) 
  d$Score <- luad_module[l, d$region] %>% unlist()
  
  d$mucinous <- factor(ifelse(d$luad_subtype == 'mucinous', 'Mucinous', "Others"), levels = c("Others", "Mucinous"))
  module_fit <- lme4::lmer(Score ~ mucinous + (1|patient), d)
  
  
  model_summary <- car::Anova(module_fit, type = 'III')
  coef <- summary(module_fit)$coefficients
  data.frame('module' = l, 'p' = model_summary['mucinous', 'Pr(>Chisq)'], coef = coef[2, 'Estimate'])
})
mucinous_de <- do.call('rbind', mucinous_lm)
mucinous_de$fdr <- p.adjust(mucinous_de$p, 'fdr')


#### 5. generate avg expression in each subtype
eg_Pathology_AvgExpr <- lapply(unique(luad_meta$luad_subtype), function(s) {
  
  if (s == 'LUSC') {
    ### adding some LUSC samples as reference
    s_idx <- nsclc_meta$region[nsclc_meta$Histology == 'LUSC']
  } else {
    s_idx <- luad_meta$region[luad_meta$luad_subtype == s]
  }
  data.frame(rowMeans(tp_module[,s_idx, drop=F]))
  
  
})

eg_Pathology_AvgExpr <- do.call('cbind', eg_Pathology_AvgExpr)
colnames(eg_Pathology_AvgExpr) <- unique(luad_meta$luad_subtype) 
eg_Pathology_AvgExpr_scale <- t(scale(t(eg_Pathology_AvgExpr)))


#### matching rownames of the matrix to module annotation
module_anno <- fread('../../Data/Module_annotation_Reorder_1008.csv')
module_anno$Name[module_anno$Name == ''] <- 'Unannotated'
module_anno <- module_anno %>% mutate(Label = paste(Module, Name, sep = ':'))



#### 6. visualization
sig_l <- mucinous_de$module[mucinous_de$fdr < 0.05]

hm_expr <- eg_Pathology_AvgExpr_scale[sig_l, ]
rownames(hm_expr) <- module_anno[match(rownames(hm_expr), module_anno$Module), ] %>% pull(Label)

##### if we remove solid in finding mucinous define modules
hm_expr <- hm_expr[, colnames(hm_expr) != "solid"]
hm_expr <- hm_expr[grep("Ribosomal", rownames(hm_expr), invert = T), ]


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_split <- colnames(eg_Pathology_AvgExpr_scale)

col_order <- c('lepidic', 'papillary', 'acinar', 'micropapillary', 'mucinous')
hm <- Heatmap(hm_expr, col = col_fun, name = "Scale\nModule\nScore", top_annotation = NULL,  #  sig_l
              left_annotation = NULL, cluster_rows=T, show_column_names = T, cluster_columns = F, 
              column_order =col_order, column_split = NULL, 
              column_names_rot  = 45, show_row_dend = F,
              column_gap = unit(2, "mm"),
              column_names_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 8),
              height = nrow(hm_expr)*unit(6, "mm"),
              # height = unit(145, 'mm'), 
              width = ncol(hm_expr)*unit(9, "mm"),
              heatmap_legend_param = list(legend_gp = gpar(fontsize = 8))
)
hm_d <- draw(hm)


pdf('../../Figure/Figure4/TraceX_Mucinous_DE_heatmap.pdf')
hm_d
dev.off()




### 7. calculate correlation between Goblet module and all other modules
muc_meta <- nsclc_meta %>% filter(luad_subtype %in% c("mucinous")) 

muc_module <- luad_module[, muc_meta$region]

muc_cor <- lapply(rownames(muc_module), function(l) {
  if (l != 'L19') {
    cor.t <- cor.test(muc_module[l, ] %>% unlist(), muc_module['L19', ] %>% unlist())
    data.frame('module' = l, 'cor' = cor.t$estimate, 'P' = cor.t$p.value)
  }
})
muc_cor <- do.call('rbind', muc_cor)
muc_cor$Label <- module_anno[match(muc_cor$module, module_anno$Module), ] %>% pull(Label)
muc_cor$FDR <- p.adjust(muc_cor$P, 'fdr')


int_mod <- c(14, 16, 17, 18, 57, 38, 41, 55, 70)
muc_mod <- c(15, 19, 24, 25, 31, 53, 54, 61, 98, 20, 59, 79, 81, 30)
l <- paste0('L', 
            c(int_mod, muc_mod)
) %>% unique()

label_l <- paste0('L', 
                  c(79, 14:17))

pdf('../../Figure/Figure4/TraceX_L19_Module_Correlation.pdf', width = 3, height = 3)

ggplot(aes(x = cor, y = -log(FDR, 10), label = Label), data = muc_cor) +
  geom_hline(yintercept = -log(0.05, 10)) +
  geom_point(color = dplyr::case_when(
    muc_cor$module %in% l ~ "red4",
    muc_cor$FDR < 0.05 ~ 'black',
    TRUE ~ "grey"), alpha = 0.8) + xlim(-1,1) +
  ggrepel::geom_label_repel(data = dplyr::filter(muc_cor, `cor` > 0 & muc_cor$module %in% label_l), 
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 2, xlim = c(0.1, NA),
                            segment.color = 'grey50', max.overlaps = 20) + 
  ggrepel::geom_label_repel(data = dplyr::filter(muc_cor, `cor` < 0 & muc_cor$module %in% label_l), 
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            size = 2, xlim = c(-0.3, NA),
                            segment.color = 'grey50') + 
theme_classic() + 
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 5)
  ) + ggtitle(label = 'Association with L19: Goblet') + 
  labs(x= "Correlation", y = "-log10(FDR)")


dev.off()
