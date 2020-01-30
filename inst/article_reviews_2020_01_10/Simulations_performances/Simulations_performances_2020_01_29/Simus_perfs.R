rm(list = ls())

setwd("/home/hulot/Documents/packages_R/mergeTrees/inst/article_reviews_2020_01_10/Simulations_performances/Simulations_performances_2020_01_29/")

library(mergeTrees)
library(aricode)
library(svd)
library(reshape2)
library(tidyverse)
library(parallel)

source("util.R")

# ---- expand grid pour param : ----

n_tables_p = c(4, 5, 100)
n_grp_p = c(10)
grp_size_p = c(50)
n0_p = c(0)
n1_p = c(1)
n_grp_per_data_p = c(2)

grid_param = data.frame(expand.grid(n_tables_p, n_grp_p, grp_size_p,
                                    n0_p, n1_p, n_grp_per_data_p))
colnames(grid_param) = c("Tables", "Grp", "Grp_sizes", "n0", "n1", "ngrpData")

dim(grid_param)

# ---- Simu donnees groupes : test scenar 1 ----

#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information



Run_func =  function(dataSets, n_grp, grp_size, reference){
  dataSets = lapply(dataSets, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1]) # sait-on jamais.
  
  hc_methods = list()
  hc_methods$AD = averagedClustering(dataSets)
  hc_methods$DC = directClustering(dataSets)
  hc_methods$MC = mergeTrees(lapply(dataSets, FUN = function(dat) hclust(dist(dat, method = "euclidean"), method = "ward.D2")))
  
  svd_dat = svd(do.call("cbind", dataSets), nv = 0)
  
  k_svd = 4
  dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  
  hc_methods$spAD_4 = averagedClustering(dat_univar)
  hc_methods$spMC_4 = mergeTreesWard1d(dat_univar)
  hc_methods$spDC_4 = directClustering(dat_univar)
  
  # Qu ce soit 4 ou duex, les courbes sont souvent confondues. 
  # Le 4 s'en sort mieux dans les autres methodes.
  
  res_NID = lapply(hc_methods, FUN = function(hc){
    apply(cutree(hc, k = 1:(n_grp*grp_size)), 2, NID, c2 = reference)
  })
  
  # tab_scen1 = do.call("rbind", lapply(res_NID, FUN = function(vect) c(min(vect), which.min(vect)[1])))
  tab_res_NID = do.call("rbind", res_NID)
  
  return(tab_res_NID)
}



function_sim = function(n_tables, n0, n1, n_grp, grp_size, separability, n_grp_per_data, n_eval = 10, reference, ncores = 2){
  
  # --- Scenario1 : ----
  res_sim1 = mclapply(1:n_eval, FUN = function(ev){
    dataSets1 = lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario1(n0, n1, n_grp, grp_size, separability)})
    Run_func(dataSets1, n_grp, grp_size, reference)
  }, mc.cores = ncores)
  
  res_tab1 = do.call("cbind", res_sim1)
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
    dataSets2 = lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)})
    Run_func(dataSets2, n_grp, grp_size, reference)
  }, mc.cores = ncores)
  
  res_tab2 = do.call("cbind", res_sim2)
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
  
  n_tables = grid_param[i, 1]
  n_grp = grid_param[i, 2]
  grp_size = grid_param[i, 3]
  n_ind = n_grp*grp_size
  n0 = grid_param[i, 4]
  n1 = grid_param[i, 5]
  n_grp_per_data = grid_param[i, 6]
  reference <- rep(1:n_grp, each = grp_size) 
  
  sep_par = c(0.2, 0.4, 1, 1.5)
  res_scen = lapply(sep_par, FUN = function(sep){function_sim(n_tables, n0, n1, n_grp, grp_size, sep, n_grp_per_data, n_eval = 2, reference)})
  
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
    xlab("Number of groups") +
    ggtitle(paste0("n0=", n0, " - n1=", n1, " - ngroups=", n_grp, " - grp size=", grp_size, "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data))
  
  p_all
  ggsave(paste0("n0=", n0, "_n1=", n1, "_ngroups=", n_grp, "_grp size=", grp_size,
                "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data,".jpeg"), height = 10, width = 25, units = c("cm"))  
}

# # ---- Comparaison vecteurs propres : ----
# 
# # Avec un vecteur par table
# n_tables = 4
# n_grp = 10
# grp_size = 50
# n_ind = n_grp*grp_size
# n0 = 0
# n1 = 1
# n_grp_per_data = 2
# reference <- rep(1:n_grp, each = grp_size) 
# separability = 1
# 
# n_eval = 10
# ncores = 2
# 
# dataSets1_list = mclapply(1:n_eval, FUN = function(ev){
#   lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario1(n0, n1, n_grp, grp_size, separability)})
# }, mc.cores = ncores)
# 
# dataSets2_list = mclapply(1:n_eval, FUN = function(ev){
#   lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)})
# }, mc.cores = ncores)
# 
# dataSets = lapply(dataSets1_list, FUN = function(liste) lapply(liste, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1])) 
# svd_dat =lapply(dataSets, FUN = function(dat_list) svd(do.call("cbind", dat_list), nv = 0))
# k_svd = 4
# dat_univar = lapply(svd_dat, FUN = function(svd_res) as.list(data.frame(svd_res$u[,1:k_svd] %*% diag(svd_res$d[1:k_svd]))))
# 
# dataSets2 = lapply(dataSets2_list, FUN = function(liste) lapply(liste, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1])) 
# svd_dat = lapply(dataSets2, FUN = function(dat_list) svd(do.call("cbind", dat_list), nv = 0))
# k_svd = 4
# dat_univar2 = lapply(svd_dat, FUN = function(svd_res) as.list(data.frame(svd_res$u[,1:k_svd] %*% diag(svd_res$d[1:k_svd]))))
# 
# 
# # Avec deux vecteurs par table
# n_tables = 4
# n_grp = 10
# grp_size = 50
# n_ind = n_grp*grp_size
# n0 = 0
# n1 = 2
# n_grp_per_data = 2
# reference <- rep(1:n_grp, each = grp_size) 
# separability = 1
# 
# n_eval = 10
# ncores = 2
# 
# dataSets1_list2 = mclapply(1:n_eval, FUN = function(ev){
#   lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario1(n0, n1, n_grp, grp_size, separability)})
# }, mc.cores = ncores)
# 
# dataSets2_list2 = mclapply(1:n_eval, FUN = function(ev){
#   lapply(1:n_tables, FUN = function(i){set.seed(i*ev*100); scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)})
# }, mc.cores = ncores)
# 
# dataSets_2 = lapply(dataSets1_list2, FUN = function(liste) lapply(liste, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1])) 
# svd_dat =lapply(dataSets_2, FUN = function(dat_list) svd(do.call("cbind", dat_list), nv = 0))
# k_svd = 4
# dat_univar_2 = lapply(svd_dat, FUN = function(svd_res) as.list(data.frame(svd_res$u[,1:k_svd] %*% diag(svd_res$d[1:k_svd]))))
# 
# dataSets2_2 = lapply(dataSets2_list2, FUN = function(liste) lapply(liste, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1])) 
# svd_dat = lapply(dataSets2_2, FUN = function(dat_list) svd(do.call("cbind", dat_list), nv = 0))
# k_svd = 4
# dat_univar2_2 = lapply(svd_dat, FUN = function(svd_res) as.list(data.frame(svd_res$u[,1:k_svd] %*% diag(svd_res$d[1:k_svd]))))
# 
# 
# # Comparaison des vecteurs
# 
# Diff1 = lapply(1:n_eval, FUN = function(i) lapply(1:k_svd, FUN = function(x) abs(sign(dat_univar[[i]][[x]][1])*dat_univar[[i]][[x]] - sign(dat_univar[[i]][[x]][1])*dat_univar_2[[i]][[x]])))
# Diff2 = lapply(1:n_eval, FUN = function(i) lapply(1:k_svd, FUN = function(x) abs(sign(dat_univar_2[[i]][[x]][1])*dat_univar_2[[i]][[x]] - sign(dat_univar_2[[i]][[x]][1])*dat_univar2_2[[i]][[x]])))
# 
# do.call("rbind", lapply(Diff1, FUN = function(dat) unlist(lapply(dat, sum))))
# do.call("rbind", lapply(Diff2, FUN = function(dat) unlist(lapply(dat, sum))))
