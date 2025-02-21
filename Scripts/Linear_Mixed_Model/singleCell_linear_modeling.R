
library(celda)
library(lme4)
library(dplyr)
library(insight)

### 1. Load dataset
op <- '/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))

# data_sce$celda_k <- altExp(data_sce)$celda_cell_cluster
# data_sce$Histology[data_sce$Histology == 'NCSLC'] <- 'NSCLC'


### 2. Quantify module score
y <- factorizeMatrix(data_sce, type = "counts")
cell_m_counts <- y[[1]][['cell']]; cell_m_counts <- t(cell_m_counts)
colSum_sizeFactor <- colSums(t(cell_m_counts)) / 10000 
normalized_m_counts <- sweep(cell_m_counts, 1, colSum_sizeFactor, "/")
module_activity <- log2(normalized_m_counts + 1)



full_meta <- colData(data_sce)
full_meta$stage_bin <- ifelse(full_meta$Stage %in% c('GGO', 'I'), 'Early', 'Late')
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells)

#### 3. perform linear mixed model

#### 3.1 Generated corrected module score
k <- unique(full_meta$celda_k) %>% as.character()
meta <- full_meta[data_sce$celda_k %in% k, c('study', 'celda_k', 'Histology','stage_bin','sample')] ### k of luad cells , 'sizeFactor'

anova_lm <- list()
anova_res <- list()
for (l in colnames(module_activity)) {
  print(l)
  d <- data.frame('module' = module_activity[rownames(meta), l],
                  meta[, c('study','stage_bin', 'celda_k', 'sample', 'Histology')])
  
  m <- lme4::lmer(module ~ stage_bin + Histology + (1|study/sample/celda_k), data=d)

  ### get coefficient and random intercept
  # m_sum <- summary(m)$coefficients
  # sampleK_int <- coef(m)[[1]]; sample_int <- coef(m)[[2]]; study_int <- coef(m)[[3]]
  # fixed_int <- m_sum['(Intercept)', 1]
  # design <- model.matrix(m)
  # 
  # random_int <- data.frame('cell' = rownames(design),
  #                          'sample' = paste(d$sample, d$study, sep = ':'), 
  #                          'celda_k' = paste(d$celda_k, d$sample, d$study, sep = ':'))
  # random_int$k_int <- sampleK_int[random_int$celda_k, '(Intercept)']
  # random_int$sample_int <- sample_int[random_int$sample, '(Intercept)']
  # 
  # random_int$sample_effect <- random_int$sample_int - fixed_int
  # random_int$cluster_effect <- random_int$k_int - fixed_int
  # 
  # res_old <- design[, -1] %*% t(study_int[1, -1]) + residuals(m) + random_int$sample_int + random_int$k_int - fixed_int
  # #anova_res[[l]] <- res[, 1]
  
  ### using ranf function to get random intercept
  design <- model.matrix(m)
  m_sum <- summary(m)$coefficients
  random_b0 <- ranef(m) ### random intercept
  random_int <- data.frame('cell' = rownames(design),
                           'study' = d$study,
                           'sample' = paste(d$sample, d$study, sep = ':'), 
                           'celda_k' = paste(d$celda_k, d$sample, d$study, sep = ':'))
  random_int$k_int <- random_b0$`celda_k:(sample:study)`[random_int$celda_k, '(Intercept)']
  random_int$sample_int <- random_b0$`sample:study`[random_int$sample, '(Intercept)']
  random_int$fixed_int <- m_sum[1,1]
  
  res <- (design %*% as.matrix(m_sum[, 1])) + residuals(m) + random_int$sample_int + random_int$k_int #+ random_int$fixed_int
  anova_res[[l]] <- res[,1]
  
  ### get posterior distribution
  #m_sim <- arm::sim(m)
  
  ### Get summary statistics
  af <- car::Anova(m, type = 'III')
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  anova_lm[[l]] <- data.frame('module' = l, 'celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(sample:study)']/all_var), 'sample_VarPer' = as.numeric(vara$var.intercept['sample:study']/all_var), 'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                              'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                              'Residual_VarPer' = vara[['var.residual']]/all_var,
                              'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                              'Stage_pval' = af['stage_bin', 3],
                              'Stage_coef' =  m_sum['stage_binLate', 1],
                              'Histology_pval' = af['Histology', 3],
                              'Histology_coef' = m_sum['HistologyLUSC', 1])

}
rm_res_mat <- do.call('cbind', anova_res) %>% t()
rm_summary <- do.call('rbind', anova_lm)

saveRDS(rm_res_mat, 
        file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/corrected_module_score.rds'))
write.csv(rm_summary,
          file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/model_summary_AllCells.csv'))



#### 3.2 Investigate histology-associated modules
his <- c('LUSC', 'LUAD') ### not NSCLC. We donn't use it in the original LM model.
k <- unique(data_sce$celda_k[data_sce$Histology %in% his]) %>% as.character()

meta <- full_meta[data_sce$celda_k %in% k & data_sce$Histology %in% his, c('study', 'celda_k', 'Histology','stage_bin','sample')] ### k of luad cells , 'sizeFactor'


##### use the corrected modules score to perform linear regression
cr_module <- readRDS('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/corrected_module_score.rds')
module_activity <- t(cr_module)


histology_lm <- list()
for (l in colnames(module_activity)) {
  print(l)
  d <- data.frame('module' = module_activity[rownames(meta), l],
                  meta[, c('study','stage_bin', 'celda_k', 'sample', 'Histology')])
  
  m <- lme4::lmer(module ~ stage_bin + Histology + (1|study/sample/celda_k), data=d)
  m_sum <- summary(m)$coefficients
  
  
  ### Get summary statistics
  af <- car::Anova(m, type = 'III')
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  histology_lm[[l]] <- data.frame('module' = l,
                              'celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(sample:study)']/all_var), 
                              'sample_VarPer' = as.numeric(vara$var.intercept['sample:study']/all_var), 
                              'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                              'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                              'Residual_VarPer' = vara[['var.residual']]/all_var,
                              'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                              'Stage_pval' = af['stage_bin', 3],
                              'Stage_coef' =  m_sum['stage_binLate', 1],
                              'Histology_pval' = af['Histology', 3],
                              'Histology_coef' = m_sum['HistologyLUSC', 1])
  
}
histo_m_summary <- do.call('rbind', histology_lm)
write.csv(histo_m_summary, 
          file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/model_summary_Histology_CorrectedModule.csv'))





##### 2.2 luad only to identify stage-associated modules
his <- c('LUAD')
hist_k_summary <- table(full_meta$Histology, full_meta$celda_k) %>% as.data.frame()
ggplot(aes(x = Var2, y = Freq, fill = Var1), data = hist_k_summary) + geom_bar(position = 'stack', stat = "identity")

k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57)) ### perLUAD in 63 is 87.5, 12.5 NSCLC
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')

##### use the corrected modules score to perform linear regression
cr_module <- readRDS('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/corrected_module_score.rds')
module_activity <- t(cr_module)


meta <- full_meta[data_sce$celda_k %in% k & data_sce$Histology %in% his, c('study', 'celda_k', 'Histology','stage_bin','sample')] ### k of luad cells , 'sizeFactor'

stage_lm <- list()
for (l in colnames(module_activity)) {
  print(l)
  d <- data.frame('module' = module_activity[rownames(meta), l],
                  meta[, c('study','stage_bin', 'celda_k', 'sample', 'Histology')])
  
  m <- lme4::lmer(module ~ stage_bin + (1|study/sample/celda_k), data=d)
  m_sum <- summary(m)$coefficients
  
  
  ### Get summary statistics
  af <- car::Anova(m, type = 'III')
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  stage_lm[[l]] <- data.frame('module' = l,
                                  'celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(sample:study)']/all_var), 
                                  'sample_VarPer' = as.numeric(vara$var.intercept['sample:study']/all_var), 
                                  'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                                  'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                                  'Residual_VarPer' = vara[['var.residual']]/all_var,
                                  'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                                  'Stage_pval' = af['stage_bin', 3],
                                  'Stage_coef' =  m_sum['stage_binLate', 1])
  
}
stage_m_summary <- do.call('rbind', stage_lm)
write.csv(stage_m_summary, 
          file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Mixed_model_output/model_summary_Stage_CorrectedModule.csv.csv'))




##### 2.3 lusc only to identify stage-associated modules
his <- c('LUSC')
hist_k_summary <- table(full_meta$Histology, full_meta$celda_k) %>% as.data.frame()
ggplot(aes(x = Var2, y = Freq, fill = Var1), data = hist_k_summary) + geom_bar(position = 'stack', stat = "identity")

k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(53)) ### perLUAD in 63 is 87.5, 12.5 NSCLC
k <- sort(as.numeric(k))


module_activity <- log2(normalized_m_counts + 1)
meta <- full_meta[data_sce$celda_k %in% k & data_sce$Histology %in% his, c('study', 'celda_k', 'Histology','Stage','sample')] ### k of luad cells , 'sizeFactor'
meta$stage_bin <- ifelse(meta$Stage == 'III/IV', 'Late', 'Early')

### remove some samples with too few cells
skip_s <- names(table(meta$sample))[table(meta$sample) < 10]
meta <- meta %>% as.data.frame() %>% filter(!sample %in% skip_s)

stage_lm <- list()
for (l in colnames(module_activity)) {
  print(l)
  d <- data.frame('module' = module_activity[rownames(meta), l],
                  meta[, c('study','stage_bin', 'celda_k', 'sample', 'Histology')])
  
  m <- lme4::lmer(module ~ stage_bin + (1|sample), data=d) #study/sample/celda_k
  m_sum <- summary(m)$coefficients
  
  
  ### Get summary statistics
  af <- car::Anova(m, type = 'III')
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  stage_lm[[l]] <- data.frame('module' = l,
                              # 'celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(sample:study)']/all_var), 
                              'sample_VarPer' = as.numeric(vara$var.intercept['sample']/all_var),  #sample:study
                              # 'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                              'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                              'Residual_VarPer' = vara[['var.residual']]/all_var,
                              'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                              'Stage_pval' = af['stage_bin', 3],
                              'Stage_coef' =  m_sum['stage_binLate', 1])
  
}
stage_m_summary <- do.call('rbind', stage_lm)
stage_m_summary$Stage_fdr <- p.adjust(stage_m_summary$Stage_pval, 'fdr')
write.csv(stage_m_summary, 
          file.path('/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Supplementary_LME_summary_LUSC_Stage.csv'))


