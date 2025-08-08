library(celda)
library(singleCellTK)
library(data.table)
library(dplyr)
library(ggplot2)
library(eulerr)

setwd(here::here("./Scripts/Figure"))

### 1. Load dataset
## load sce object
op <- '../../Data'
data_sce <- readRDS(file.path(op, 'celda_Module_CellCluster_Reorder_08.rds'))


x <- factorizeMatrix(data_sce, type = "proportion")
cell_probs <- x[[1]][['cell']]; cell_probs <- t(cell_probs)
module_prob <- cell_probs


### Define lineage module
#### NKX2-1 
#### AT2
#### AT1#
#### Club
#### SCGB3A2

lineage_modules <- c(paste0('L',c(
  10:12, ##LUAD
  1:4,  ## Surfactant
  5:6,  ## AT2
  7:9,  ## AT1
  14, ## Sox2
  19, ## MUC5AC
  20:22, ## Club
  26, ##SCGB3A2
  28 ##SCGB3A1
))
  
)
target_modules <- c(lineage_modules)


#### select the LUAD cells
his <- 'LUAD'
full_meta <- colData(data_sce)
full_meta$stage_bin <- ifelse(full_meta$Stage %in% c('GGO', 'I'), 'Early', 'Late')
k_summary <- full_meta %>% as.data.frame() %>% group_by(celda_k) %>% summarise(nCells = n(),
                                                                               nLUAD = sum(Histology == 'LUAD'),
                                                                               nLUSC = sum(Histology == 'LUSC'),
                                                                               perLUAD = nLUAD / nCells,
                                                                               perLUSC = nLUSC / nCells)


k <- k_summary$celda_k[k_summary[[paste0('per', his)]] > 0.9]
# add some NSCLC cluster that are similar to LUAD in transcriptomic data
k <- union(k, c(12, 29, 34, 36, 41, 45, 49, 57))
k <- sort(as.numeric(k))
his <- c('LUAD', 'NSCLC')
select_s <- unique(full_meta$sample[full_meta$celda_k %in% k & full_meta$Histology %in% his])
select_c <- rownames(full_meta)[full_meta$sample %in% select_s]





### 2. Define module on/off

ip <- '../../Mixed_model_output/'
luad_lm <- fread(file.path(ip, 'model_summary_Stage.csv'))
skip_m <- luad_lm$module[luad_lm$study_VarPer > 0.3]

cr_module <- readRDS(file.path(ip, 'corrected_module_score.rds'))
luad_res <- SingleCellExperiment(assay = list(corrected_module = cr_module[!rownames(cr_module) %in% skip_m, select_c]), #
                                 rowData = NULL,
                                 colData = colData(data_sce)[select_c,])


cf <- 0.005
assay(luad_res, 'module_on') <- 1* t(module_prob[colnames(luad_res), rownames(luad_res)] > cf) #target_module
lm_on <- t(assay(luad_res, 'module_on')) %>% as.data.frame()


# #### checking whether cf = 0.01 cutoff is too high for different module score
# mod_prob_meta <- module_prob[colnames(luad_res), paste0('L', c(1:4, 9, 26))] %>% as.data.frame()
# mod_prob_meta <- reshape2::melt(mod_prob_meta)
# colnames(mod_prob_meta) <- c('Module', 'Probability')
# ggplot(aes(x = Probability, color = Module), data = mod_prob_meta) + geom_density() + 
#   geom_vline(xintercept = cf, color = 'red3') + xlim(0, 0.15) + ylim(0, 500)

# plot(density(module_prob[colnames(luad_res), 'L1'])) + abline(v = cf, col = 'red')
# plot(density(module_prob[colnames(luad_res), 'L2'])) + abline(v = cf, col = 'red')
# plot(density(module_prob[colnames(luad_res), 'L3'])) + abline(v = cf, col = 'red')
# plot(density(module_prob[colnames(luad_res), 'L4'])) + abline(v = cf, col = 'red')
# plot(density(module_prob[colnames(luad_res), 'L9'])) + abline(v = cf, col = 'red')
# plot(density(module_prob[colnames(luad_res), 'L26'])) + abline(v = cf, col = 'red')



### 3. Calculate percentage of module that are on for each lineage subgroup

lm_on <- lm_on %>% mutate(
  #adenoLineage = 1*(rowSums(lm_on[, paste0('L', 10:12)]) != 0),
  SFTPC = 1*(rowSums(lm_on[, paste0('L', 1), drop=F]) != 0),
  `SFTPA1/2` = 1*(rowSums(lm_on[, paste0('L', 2), drop=F]) != 0),
  #`SFTPB/D` = 1*(rowSums(lm_on[, paste0('L', 3:4), drop=F]) != 0),
  `SFTPB` = 1*(rowSums(lm_on[, paste0('L', 3), drop=F]) != 0),
  `SFTPD` = 1*(rowSums(lm_on[, paste0('L', 4), drop=F]) != 0),
  # `AGR3` = 1*(rowSums(lm_on[, paste0('L', 8), drop=F]) != 0),
  `AGER` = 1*(rowSums(lm_on[, paste0('L', 9), drop=F]) != 0),
  
  #AT2 = 1*(rowSums(lm_on[, paste0('L', 5:6)]) != 0),
  #AT1 = 1*(rowSums(lm_on[, paste0('L', 7:9)]) != 0),
  #squamousLineage = 1*(rowSums(lm_on[, paste0('L', 14), drop=F]) != 0),
  #mucinous = 1*(rowSums(lm_on[, paste0('L', 19), drop=F]) != 0),
  #Club = 1*(rowSums(lm_on[, paste0('L', 20:22)]) != 0),
  SCGB3A2 = 1*(rowSums(lm_on[, paste0('L', 26), drop=F]) != 0),
  # SCGB3A1 = 1*(rowSums(lm_on[, paste0('L', 28), drop=F]) != 0)
)


m_x <- lapply(lm_on, function(l) {
  rownames(lm_on)[l == 1]
})

m <- grep('^L',  colnames(lm_on), invert = T, value=T)

md <- m_x[m]
set.seed(2020)
fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,
fill_col <- setNames(RColorBrewer::brewer.pal(length(m), "Set3"), m)


pdf("../..//Figure/Figure4/Figure4_NSCLC_lineageVenn_probbased_cf0.005.pdf", width = 4, height = 4)
plot(fit,
     quantities = list(type = c("percent"), cex = 0.5), #TRUE #, "counts", "percent"
     fill = fill_col, #RColorBrewer::brewer.pal(length(m), "Set3"),
     lty = 1:length(m),
     labels = list(cex = 0.5),
     #labels = NULL,
     legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"), 
                   labels = rownames(fit$ellipses), 
                   cex = 0.5))
dev.off()



# check what intersections are missing
miss_op <- data.frame('OP' = names(fit$original.values), 
                      'Number' = fit$original.values)
miss_op$OP_Per <- miss_op$Number / sum(miss_op$Number) * 100


miss_op <- miss_op %>% as.data.frame()
write.csv(miss_op, 
          "../..//Table/Supplementary_NSCLC_VennDiagram.csv")




### quantify # of module marker express in each cell each sample
marker_on <- lm_on[, c('SFTPA1/2', 'SFTPB', 'SFTPD', 'AGER', 'SCGB3A2')]
marker_on$num_marker <- rowSums(marker_on)
marker_on <- cbind(marker_on, 
                   colData(data_sce)[rownames(marker_on), c("new_sample", "Stage", "study", "celda_k")])
marker_on <- marker_on %>% mutate(
  Stage = case_when(
    Stage %in% c("GGO", "I") ~ "GGO/I", 
    Stage %in% c("II", "III/IV") ~ "II/III/IV"
  )
)


marker_on_summary <- marker_on %>% group_by(new_sample) %>% 
  summarise(avg_num_marker = mean(num_marker), Stage = unique(Stage))

pdf("../..//Figure/Figure4/Figure3_scRNA_Number_Marker_modules.pdf", width = 4, height = 4)
print(
  ggplot(aes(x = Stage, y = avg_num_marker, fill = Stage), data = marker_on_summary) + 
    geom_boxplot() + theme_bw() + 
    geom_jitter() + ylab("Average number of marker modules")
)
dev.off()



### using LME model 

m <- lme4::lmer(num_marker ~ Stage + (1|study/new_sample/celda_k), data=marker_on)

### Get summary statistics
af <- car::Anova(m, type = 'III')
m_sum <- summary(m)$coefficients

#### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
vara <- insight::get_variance(m, 'all', verbose = T)
random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum


m_summary <- data.frame('celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(new_sample:study)']/all_var), 'sample_VarPer' = as.numeric(vara$var.intercept['new_sample:study']/all_var), 'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                            'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                            'Residual_VarPer' = vara[['var.residual']]/all_var,
                            'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                            'Stage_pval' = af['Stage', 3],
                            # 'StageI_coef' =  m_sum['StageI', 1], 
                            # 'StageII_coef' =  m_sum['StageII', 1], 
                            'LateStage_coef' =  m_sum['StageII/III/IV', 1]
                            )


library(jtools)
m_summary_plot <- plot_summs(m)
# summ(m, confint = TRUE, digit = 3)


#### testing individual modules on
marker_coef <- list()
marker_models <- list()
for (marker in c('SFTPA1/2', 'SFTPB', 'SFTPD', 'AGER', 'SCGB3A2')) {
  marker_on <- lm_on[, marker, drop = F]
  marker_on$num_marker <- rowSums(marker_on)
  marker_on <- cbind(marker_on, 
                     colData(data_sce)[rownames(marker_on), c("new_sample", "Stage", "study", "celda_k")])
  
  marker_on <- marker_on %>% mutate(
    Stage = case_when(
      Stage %in% c("GGO", "I") ~ "GGO/I", 
      Stage %in% c("II", "III/IV") ~ "II/III/IV"
    )
  )
  
  # marker_on <- marker_on %>% group_by(new_sample) %>% 
  #   summarise(avg_num_marker = mean(num_marker), Stage = unique(Stage))
  
  sub_m <- lme4::lmer(num_marker ~ Stage + (1|study/new_sample/celda_k), data=marker_on)
  marker_models[[marker]] <- sub_m
  ### Get summary statistics
  af <- car::Anova(sub_m, type = 'III')
  m_sum <- summary(sub_m)$coefficients
  
  #### using interclass correlation to calculate % variability explained by covariates: https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  vara <- insight::get_variance(sub_m, 'all', verbose = T)
  random_var <-  if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}
  all_var <- unlist(c(vara[c('var.fixed', 'var.residual')], random_var)) %>% sum
  
  
  m_summary <- data.frame('celda_VarPer' = as.numeric(vara$var.intercept['celda_k:(new_sample:study)']/all_var), 'sample_VarPer' = as.numeric(vara$var.intercept['new_sample:study']/all_var), 'study_VarPer' = as.numeric(vara$var.intercept['study']/all_var),
                          'Random_VarPer'  = if (is.null(vara[['var.random']])) {sum(vara[['var.intercept']])} else {vara[['var.random']]}/all_var,
                          'Residual_VarPer' = vara[['var.residual']]/all_var,
                          'Fixed_VarPer' = vara[['var.fixed']]/all_var,
                          'Stage_pval' = af['Stage', 3],
                          # 'StageI_coef' =  m_sum['StageI', 1], 
                          # 'StageII_coef' =  m_sum['StageII', 1], 
                          # 'StageIII/IV_coef' =  m_sum['StageIII/IV', 1], 
                          'LateStage_coef' =  m_sum['StageII/III/IV', 1],
                          
                          "marker" = marker
  )
  
  marker_coef[[marker]] <- m_summary
}
marker_coef <- do.call('rbind', marker_coef)

pdf("../..//Figure/Figure4/Figure3_scRNA_Number_Marker_modules_coefficient.pdf", width = 6, height = 6)
print(
  plot_summs(m, marker_models$`SFTPA1/2`, marker_models$SFTPB, marker_models$SFTPD, marker_models$AGER, marker_models$SCGB3A2, 
             model.names = c("AllModules", c('SFTPA1/2', 'SFTPB', 'SFTPD', 'AGER', 'SCGB3A2')))
)
dev.off()


### 4. Calculate percentage of NAPSA, NKX2-1, and linage modules

# lm_on <- lm_on %>% mutate(
#   Lineage = 1*(rowSums(lm_on[, paste0('L', c(1:4, 26, 28))]) != 0),
#   RNASE1 = 1*(rowSums(lm_on[, paste0('L', 11), drop=F]) != 0),
#   `NAPSA` = 1*(rowSums(lm_on[, paste0('L', 10), drop=F]) != 0),
#   `NKX2-1` = 1*(rowSums(lm_on[, paste0('L', 12), drop=F]) != 0),
# )
# 
# 
# m_x <- lapply(lm_on, function(l) {
#   rownames(lm_on)[l == 1]
# })
# 
# library(eulerr)
# m <- c('NAPSA', 'RNASE1', 'NKX2-1', 'Lineage')
# 
# md <- m_x[m]
# fit <-euler(md, quantities = TRUE,  shape = "ellipse", loss  = 'square', loss_aggregator = 'max') #,
# 
# 
# plot(fit,
#      quantities = list(type = c("counts", "percent"), cex = 1), #TRUE #, "counts", "percent"
#      fill = RColorBrewer::brewer.pal(length(m), "Set3"),
#      lty = 1:length(m),
#      labels = list(cex = 1.5),
#      #labels = NULL,
#      legend = list(cols = RColorBrewer::brewer.pal(length(m), "Set2"), 
#                    labels = rownames(fit$ellipses)))
# 
# 
# 
# 
# # check what intersections are missing
# miss_op <- data.frame('OP' = names(fit$fitted.values)[fit$fitted.values == 0],
#                       'NumOP' = fit$original.values[fit$fitted.values == 0])
# miss_op$OP_Per <- miss_op$NumOP / length(Reduce('union', m_x[m])) * 100
# missOP <- miss_op %>% filter(OP_Per != 0) %>% arrange(OP_Per)
# sum(missOP$NumOP); sum(missOP$OP_Per)

