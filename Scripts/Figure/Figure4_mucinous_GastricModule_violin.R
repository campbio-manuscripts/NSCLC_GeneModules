
#### script to perform analysis using Tracex dataset
library(fst)
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


### functions
eigenGene <- function(expr, module, gs, scale = T) {
  g <- gs[[module]]
  g <- g[g %in% colnames(expr)]
  if (length(g) < 2) {return()}
  
  mod_expr <- expr[, g]
  gsd <- colSds(mod_expr)
  g <- g[gsd != 0]
  if (length(g) < 2) {return()}
  
  mod_expr <- expr[, g]
  pca <- prcomp(mod_expr, scale = T, center = T)
  eg <- pca$x[, 1]
  names(eg) <- rownames(expr)
  per_var <- pca$sdev[1]^2 / sum(pca$sdev^2)
  return(list(eg, per_var))
}



### quantify Pascal's module on the tracex dataset
#hurley_csc_module <- fread("/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Analysis/Validation_Datasets/Hurley_CellStemCell_2021/module_features.csv")
kathiriya_ncb_module <- fread("/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Analysis/Validation_Datasets/Kathiriya_NatureCellBiology_2022/output/module_features.csv")



### 1. load the module
gastric_mod <- kathiriya_ncb_module$L79; gastric_mod <- gastric_mod[gastric_mod != ""]
basal_mod <- kathiriya_ncb_module$L91; basal_mod <- basal_mod[basal_mod != ""]



### 2. load data
tp <- "/restricted/projectnb/pulmseq/dylan/TRACERx_processed/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA"

tp_meta <- read_fst(file.path(tp, '2022-10-18all_metadata.fst'))
rsem_tpm <- read_fst(file.path(tp, '2022-10-17_rsem_tpm_mat.fst'))

gene_id <- rsem_tpm$gene_id
rsem_tpm$gene_id <- NULL


### 3. subset to the sample that has annotation in tp_meta as LUAD and LUSC
nsclc_meta <- tp_meta %>% filter(Histology %in% c("LUAD", 'LUSC', 'Collision LUAD and LUSC', 'Adenosquamous carcinoma'))
tp_id <- colnames(rsem_tpm)[colnames(rsem_tpm) %in% nsclc_meta$region] ###541 samples in total
sub_tpm <- rsem_tpm[, tp_id]
sub_log2tpm <- log2(sub_tpm + 1)
rownames(sub_log2tpm) <- gene_id


### 4. calculate module score
gs_list <- list('Basal' = basal_mod, 
                'Gastric' = gastric_mod)

eigen_gene <- lapply(names(gs_list), function(l){print(l); eigenGene(t(sub_log2tpm), l, gs_list, scale=T)}) ### sample as row
module_score <- lapply(eigen_gene, function(x) {x[[1]]})
module_score <- do.call('rbind', module_score)
rownames(module_score) <- names(gs_list)


### 5. Eigengene (PC1) is chosen arbitrarily. Need to align the eigengene direction based on the average gene expression
for (l in rownames(module_score)) {
  g <- gs_list[[l]]
  avg_g <- colMeans(sub_log2tpm[g[g %in% rownames(sub_log2tpm)], colnames(module_score)])
  if (cor(module_score[l,], avg_g) < 0) {module_score[l,] <- module_score[l,] * -1 }
  
}




### 6. Plot these two modules across luad histological subtype
luad_meta <- nsclc_meta %>% filter(Histology == "LUAD", 
                                   !luad_subtype %in% c("not possible (too crushed)", "cribriform")) ### mucious later, "mucinous"
luad_meta <- luad_meta %>% filter(!is.na(luad_subtype), region %in% colnames(module_score))


luad_meta$luad_subtype <- factor(luad_meta$luad_subtype, levels =  c('lepidic', 'papillary', 'acinar', 'micropapillary', 'solid', 'mucinous'))
luad_meta$Gastric <- module_score["Gastric", luad_meta$region] %>% unlist()
luad_meta$Basal <- module_score["Basal", luad_meta$region] %>% unlist()


ll <- 'Gastric'
p1 <- ggplot(aes(x = luad_subtype, y = Gastric, fill = luad_subtype), data = luad_meta) + geom_violin() + ggtitle(ll)

ll <- "Basal"
p2 <- ggplot(aes(x = luad_subtype, y = Basal, fill = luad_subtype), data = luad_meta) + geom_violin() + ggtitle(ll)

cowplot::plot_grid(p1, p2, ncol = 2)
