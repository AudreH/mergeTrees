rm(list = ls())

setwd("/home/hulot/Documents/packages_R/mergeTrees/inst/article_reviews_2020_01_10/Simulations_performances_2020_01_16/")

library(mergeTrees)
library(aricode)
library(svd)
library(reshape2)
library(tidyverse)
library(dplyr)
library(parallel)

source("util.R")

# ---- Simu donnees groupes : test scenar 1 ----

#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information

n_tables = 5
n_grp = 10
grp_size = 50
n_ind = n_grp*grp_size
n0 = 500
n1 = 50
n_grp_per_data = 3
reference <- rep(1:n_grp, each = grp_size) 

Run_func =  function(dataSets, n_grp, grp_size, reference){
  # dataSets = lapply(dataSets, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1]) # sait-on jamais.
  
  hc_methods = list()
  hc_methods$AD = averagedClustering(dataSets)
  hc_methods$DC = directClustering(dataSets)
  hc_methods$MC = mergeTrees(lapply(dataSets, FUN = function(dat) hclust(dist(dat, method = "euclidean"), method = "ward.D2")))
  
  svd_dat = svd(do.call("cbind", dataSets), nv = 0)
  
  # k_svd = 2
  # dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  # 
  # hc_methods$spAD_2 = averagedClustering(dat_univar)
  # hc_methods$spMC_2 = mergeTreesWard1d(svd_dat)
  # hc_methods$spDC_2 = directClustering(dat_univar)
  
  # k_svd = 3
  # dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  # 
  # hc_methods$spAD_3 = averagedClustering(dat_univar)
  # hc_methods$spMC_3 = mergeTreesWard1d(svd_dat)
  # hc_methods$spDC_3 = directClustering(dat_univar)
  
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
  gg_tab1$Scenario = "Scenario1"
  
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
  gg_tab2$Scenario = "Scenario2"
  
  p2 <- ggplot(gg_tab2, aes(N, NID, color = Method)) + 
    geom_line() + geom_point(aes(shape = Method)) + ggtitle(paste0("Scenario2 - sep = ", separability))
  
  
  return(list(gg_tab1 = gg_tab1, res_sim1 = res_sim1, p1 = p1,
              gg_tab2 = gg_tab2, res_sim2 = res_sim2, p2 = p2))
}

sep_par = c(0.1, 0.2, 0.4, 1)
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
  geom_point(aes(shape = Method), size = 1) +
  scale_shape_manual(values = c(1:nb_cols)) +
  facet_grid(Scenario ~ Separability) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  # scale_color_discrete(NULL)  +
  guides(shape = guide_legend(nrow = 1, ncol = 6), colour = NULL) +
   xlim(c(0, 25)) +
  xlab("Number of groups")
p_all
ggsave("p_all.pdf", height = 1.5*10, width = 1.5*15, units = c("cm"))


# # Juste le pas spectral
# p_non = ggplot(res_plot[res_plot$Spectral == "Non",],  aes(N, NID, color = Method)) +
#   geom_line() + 
#   geom_point(aes(shape = Method), size = 0.5) +
#   facet_grid(Scenario ~ Separability) +
#   theme_bw() + 
#   theme(legend.position="bottom", legend.box = "horizontal") +
#   guides(color = guide_legend(nrow = 1)) + xlim(c(0,25)) +
#   xlab("Number of groups")
# p_non
# 
# ggsave("p_spectral.pdf", height = 1.5*7.5, width = 2*15, units = c("cm"))
# 
# # Juste le spectral
# p_spectral = ggplot(res_plot[res_plot$Spectral == "Spectral",],  aes(N, NID, color = Method)) +
#   geom_line() + 
#   geom_point(aes(shape = Method), size = 0.5) +
#   facet_grid(Scenario ~ Separability) +
#   theme_bw() + 
#   theme(legend.position="bottom", legend.box = "horizontal") +
#   guides(color = guide_legend(nrow = 1)) + xlim(c(0,25)) +
#   xlab("Number of groups")
# p_spectral
# 
# ggsave("p_spectral.pdf", height = 1.5*7.5, width = 2*15, units = c("cm"))

