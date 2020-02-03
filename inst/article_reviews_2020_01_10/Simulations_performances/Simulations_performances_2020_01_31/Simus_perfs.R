rm(list = ls())

setwd("/home/hulot/Documents/packages_R/mergeTrees/inst/article_reviews_2020_01_10/Simulations_performances/Simulations_performances_2020_01_31/")

library(mergeTrees)
library(aricode)
library(svd)
library(reshape2)
library(tidyverse)
library(parallel)

source("util.R")

# ---- expand grid pour param : ----

#' @param n1 nombre de vecteurs dans la table.

# n_tables_p = c(5, 10)
# n_tables_info_p = c(3, 5, 10)
# n_grp_p = c(5, 10)
# grp_size_p = c(50)
# n1_p = c(10, 50, 100) 
# n_grp_per_data_p = c(2, 4)

n_tables_p = 10
n_tables_info_p = c(5, 7, 10)
n_grp_p = 5
grp_size_p = c(50)
n1_p = c(100, 500)
n_grp_per_data_p = c(2, 3, 4)

grid_param = data.frame(expand.grid(n_tables_p, n_grp_p, grp_size_p,
                                    n1_p, n_grp_per_data_p, n_tables_info_p))
colnames(grid_param) = c("Tables", "Grp", "Grp_sizes", "n1", "ngrpData", "ntablesinfo")
grid_param = grid_param[which(grid_param$ntablesinfo<grid_param$Tables),]

dim(grid_param)

# sep_par = c(seq(0.2, 1.5, by = 0.1))
sep_par = c(0.4, 0.8, 1, 1.5, 2)
neval = 2

# ---- Simu donnees groupes : Fonction de run sur data ----

Run_func =  function(dataSets, n_grp, grp_size, reference){
  dataSets = lapply(dataSets, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1]) # sait-on jamais.
  
  hc_methods = list()
  hc_methods$AD = averagedClustering(dataSets)
  hc_methods$DC = directClustering(dataSets)
  hc_methods$MC = mergeTrees(lapply(dataSets, FUN = function(dat) hclust(dist(dat, method = "euclidean"), method = "ward.D2")))
  
  svd_dat = svd(do.call("cbind", dataSets), nv = 0)
  
  k_svd = 2
  dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  
  hc_methods$spAD_2 = averagedClustering(dat_univar)
  hc_methods$spMC_2 = mergeTreesWard1d(dat_univar)
  hc_methods$spDC_2 = directClustering(dat_univar)
  
  # k_svd = 4
  # dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  # 
  # hc_methods$spAD_4 = averagedClustering(dat_univar)
  # hc_methods$spMC_4 = mergeTreesWard1d(dat_univar)
  # hc_methods$spDC_4 = directClustering(dat_univar)
  # 
  # k_svd = 8
  # dat_univar= as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  # 
  # hc_methods$spAD_8 = averagedClustering(dat_univar)
  # hc_methods$spMC_8 = mergeTreesWard1d(dat_univar)
  # hc_methods$spDC_8 = directClustering(dat_univar)
  
  res_NID = lapply(hc_methods, FUN = function(hc){
    apply(cutree(hc, k = 1:(n_grp*grp_size)), 2, NID, c2 = reference)
  })
  
  tab_res_NID = do.call("rbind", res_NID)
  
  return(tab_res_NID)
}


# ---- Fonction de simulation/resultats : ----

function_sim = function(n_tables,n1, n_grp, grp_size, separability, n_grp_per_data, n_eval = 10, reference,  ntab_info, ncores = 2){
  
  # --- Scenario1 : ----
  res_sim1 = mclapply(1:n_eval, FUN = function(ev){
    dataSets1 = vector(mode = "list", length = n_tables)
    dataSets1[1:ntab_info] = lapply(1:ntab_info, FUN = function(i){set.seed(i*ev*100); scenario1( n1, n_grp, grp_size, separability)}) 
    if(ntab_info!=n_tables){
      dataSets1[(ntab_info+1):n_tables] = lapply((ntab_info+1):n_tables, FUN = function(i){set.seed(i*ev*100); matrix(rnorm(grp_size * n_grp * (n1), sd = 1),  ncol = (n1), nrow = grp_size * n_grp)})
    }
    Run = Run_func(dataSets1, n_grp, grp_size, reference)
    return(list(dataSets1 = dataSets1, 
               Run = Run))
  }, mc.cores = ncores)

  res_tab1 = do.call("cbind", lapply(res_sim1, FUN = function(x) x$Run))
  gg_tab1 = melt(res_tab1)
  colnames(gg_tab1) = c("Method", "N", "NID")
  
  gg_tab1 = gg_tab1 %>% group_by(Method, N) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  gg_tab1$Separability = separability
  gg_tab1$Scenario = "Scenario 1"
  
  p1 <- ggplot(gg_tab1, aes(N, NID, color = Method)) + 
    geom_line() + geom_point(aes(shape = Method)) + ggtitle(paste0("Scenario1 - sep = ", separability))
  
  # --- Scenario2 : ----
  
  res_sim2 = mclapply(1:n_eval, FUN = function(ev){
    # dataSets2 = lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)})
    dataSets2 = vector(mode = "list", length = n_tables)
    dataSets2[1:ntab_info] = lapply(1:ntab_info, FUN = function(i){set.seed(i*ev*100); scenario2( n1, n_grp, n_grp_per_data, grp_size, separability) }) 
    if(ntab_info!=n_tables){
      dataSets2[(ntab_info+1):n_tables] = lapply((ntab_info+1):n_tables, FUN = function(i){set.seed(i*ev*100); matrix(rnorm(grp_size * n_grp * (n1), sd = 1),  ncol = (n1), nrow = grp_size * n_grp)})
    }
    Run = Run_func(dataSets2, n_grp, grp_size, reference)
    return(list(dataSets2 = dataSets2, 
                Run = Run))
  }, mc.cores = ncores)
  
  res_tab2 = do.call("cbind", lapply(res_sim2, FUN = function(x) x$Run))
  gg_tab2 = melt(res_tab2)
  colnames(gg_tab2) = c("Method", "N", "NID")
  
  gg_tab2 = gg_tab2 %>% group_by(Method, N) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  gg_tab2$Separability = separability
  gg_tab2$Scenario = "Scenario 2"
  
  p2 <- ggplot(gg_tab2, aes(N, NID, color = Method)) + 
    geom_line() + geom_point(aes(shape = Method)) + ggtitle(paste0("Scenario2 - sep = ", separability))
  
  
  return(list(gg_tab1 = gg_tab1, res_sim1 = res_sim1, p1 = p1,
              gg_tab2 = gg_tab2, res_sim2 = res_sim2, p2 = p2))
}

for(i in 1:nrow(grid_param)){
  
  # i = 1
  n_tables = grid_param[i, 1]
  n_grp = grid_param[i, 2]
  grp_size = grid_param[i, 3]
  n_ind = n_grp*grp_size
  # n0 = grid_param[i, 4]
  n1 = grid_param[i, 4]
  n_grp_per_data = grid_param[i, 5]
  n_tab_info = grid_param[i,6]
  reference <- rep(1:n_grp, each = grp_size) 

  res_scen = lapply(sep_par, FUN = function(sep){function_sim(n_tables, n1, n_grp, grp_size, sep, n_grp_per_data, n_eval = neval, reference, ntab_info = n_tab_info)})
  
  res_plot = do.call("rbind", c(lapply(res_scen, FUN = function(res) res$gg_tab1), lapply(res_scen, FUN = function(res) res$gg_tab2)))
  
  res_plot$Spectral = "Non"
  res_plot$Spectral[grep("sp", res_plot$Method)] = "Spectral"
  
  
  library(RColorBrewer)
  # Graphe pas lisible
  nb_cols = length(unique(res_plot$Method))
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb_cols)
  # mycolors <- rainbow(nb_cols)
  
  p_all = ggplot(res_plot,  aes(N, NID, colour = Method)) +
    geom_line() +
    geom_point(aes(shape = Method), size = 1.5) +
    scale_shape_manual(values = c(1:nb_cols)) +
    facet_grid(Scenario ~ Separability) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    # scale_color_discrete(NULL)  +
    guides(shape = guide_legend(nrow = 1, ncol = nb_cols), colour = NULL) +
    xlim(c(0, 25)) +
    ylim(c(0,1)) +
    xlab("Number of groups") +
    ggtitle(paste0("_n1=", n1, "_ngroups=", n_grp, "_grp size=", grp_size,
                   "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data, "_n_info=" , n_tab_info))
  
  print(p_all)
  ggsave(plot = p_all, filename = paste0("_n1=", n1, "_ngroups=", n_grp, "_grp size=", grp_size,
                "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data, "_n_info=" , n_tab_info, ".jpeg"), height = 15, width = 25, units = c("cm"))  
}

