module_summary <- function(module_score, module_prob, cf, clusters, his, data_sce, module) {
  sample_module <- lapply(clusters, function(k) { #include all cells within that cluster, if it's relatively homogeneous
    cid <- colnames(data_sce)[data_sce$celda_k == k ] ##& data_sce$Histology %in% his
    
    d <- data.frame('cluster' = k, 
                    'per.expr' = colSums(as.matrix(module_prob[cid, module])> cf ) /length(cid), 
                    'avg.expr' = colMeans(as.matrix(module_score[cid, module])))
    d$module <- module
    d
  })
  
  sample_module <- do.call('rbind', sample_module)
  
  
  module_expr <- sample_module %>% dplyr::select(-per.expr) %>% 
    pivot_wider(names_from = cluster, values_from = avg.expr) %>% as.data.frame()
  rownames(module_expr) <- module_expr$module
  module_expr <- module_expr[,-1] %>% as.matrix()
  
  module_per <- sample_module %>% dplyr::select(-avg.expr) %>% 
    pivot_wider(names_from = cluster, values_from = per.expr) %>% as.data.frame()
  rownames(module_per) <- module_per$module
  module_per <- module_per[,-1] %>% as.matrix()
  
  return(list(sample_module, module_expr, module_per))
}


dot_plot_param <- function(col_scale, dot_size_label, cf = 0.005) {
  col_fun = circlize::colorRamp2(col_scale,  c('snow', 'yellow', 'red2')) #viridis, c(1,10, 20) #brewer.pal(n = 3, name = "YlOrRd") , rocket(20)[c(20, 10, 1)]
  
  layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= (pindex(module_per, i, j)) * unit(2, "mm"), #percent_size
                gp = gpar(fill = col_fun(pindex(module_expr, i, j)), col = NA))}
  
  
  lgd_list = list(
    Legend( labels = dot_size_label, title = paste0("Percent Cells with \nModule Detected > ", cf*100, "%"),
            graphics = list(
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.1 * unit(2.5, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2.5, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2.5, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2.5, "mm"),
                                               gp = gpar(fill = "black")))#,
            # function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
            #                                  gp = gpar(fill = "black")))
    ))
  
  
  
  #return(list(col_fun, cell_fun, layer_fun, lgd_list))
  return(list(col_fun, layer_fun, lgd_list))
}



eigenGene <- function(expr, module, scale = T) {
  
  g <- module[module %in% rownames(expr)]
  
  if (length(g) < 2) {return()}
  
  mod_expr <- expr[g, ] %>% t()
  pca <- prcomp(mod_expr, scale = scale, center = T)
  eg <- pca$x[, 1]
  names(eg) <- rownames(mod_expr)
  return(eg)
}

feature_summary <- function(feature, data_sce, module_expr, k, his, color_brewer = 'YlOrRd') {
  cluster_table <- table(data_sce$celda_k[data_sce$celda_k %in% k & data_sce$Histology %in% his], 
                         data_sce[[feature]][data_sce$celda_k %in% k & data_sce$Histology %in% his])
  cluster_table <- cluster_table / rowSums(cluster_table)
  cluster_table <- cluster_table[colnames(module_expr), ]
  
  color <- brewer.pal(n = ncol(cluster_table), name = color_brewer)
  return(list(cluster_table, color))
}


hm_summary <- function(module_score, module_prob, cf, clusters, his, data_sce, module, clu_col = 'cluster') {
  sample_module <- lapply(clusters, function(k) { #include all cells within that cluster, if it's relatively homogeneous
    cid <- colnames(data_sce)[data_sce[[clu_col]] == k ] ##& data_sce$Histology %in% his
    
    d <- data.frame('cluster' = k, 
                    'per.expr' = colSums(as.matrix(module_prob[cid, module])> cf ) /length(cid), 
                    'avg.expr' = colMeans(as.matrix(module_score[cid, module])))
    #d$module <- rownames(d)
    d$module <- module
    d
  })
  
  sample_module <- do.call('rbind', sample_module)
  
  
  module_expr <- sample_module %>% dplyr::select(-per.expr) %>% 
    pivot_wider(names_from = cluster, values_from = avg.expr) %>% as.data.frame()
  rownames(module_expr) <- module_expr$module
  module_expr <- module_expr[,-1] %>% as.matrix()
  
  module_per <- sample_module %>% dplyr::select(-avg.expr) %>% 
    pivot_wider(names_from = cluster, values_from = per.expr) %>% as.data.frame()
  rownames(module_per) <- module_per$module
  module_per <- module_per[,-1] %>% as.matrix()
  
  return(list(sample_module, module_expr, module_per))
} 