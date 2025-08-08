#### generate an example inferCNV heatmap
library(infercnv)

fp <- '/restricted/projectnb/camplab/home/rzhong/projects/luad_preprocess/1_inferCNV/Chen_Zhang_2020/Chen_Zhang_2020_NSCLC-10'
cnv_obj <- readRDS(file.path(fp, 'run.final.infercnv_obj'))
plot_cnv(cnv_obj, 
         out_dir = '/restricted/projectnb/camplab/projects/20230801_LungMetaSC/Figure/Figure1', 
         title = "Chen_Zhang_2020_NSCLC-10", 
         output_format = 'pdf', 
         cluster_by_groups=FALSE, #TRUE
         cluster_references=FALSE,
         k_obs_groups = 10)
