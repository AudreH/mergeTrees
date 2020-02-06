rm(list = ls())

library(mergeTrees)
library(aricode)
library(svd)
library(reshape2)
library(tidyverse)  
library(RColorBrewer)
library(parallel)

source("util.R")

# ---- Parameters Search Grid : ----

n_tables_p = c(20)
n_tables_info_p = c(15, 20)
n_grp_p = c(5, 10)
grp_size_p = c(10, 25)
n0_p = c(100) 
n1_p = c(50, 100) 
n_grp_per_data_p = c(2, 3, 4)
sep_par = c(0.2, 0.4, 0.8, 1, 2)
spectral_p = c(2, 4)
sd_p = c(1, 5, 10)

grid_param = data.frame(expand.grid(n_tables_p, n_grp_p, grp_size_p,
                                    n0_p, n1_p, n_grp_per_data_p, n_tables_info_p, spectral_p,
                                    sd_p))
colnames(grid_param) = c("Tables", "Grp", "Grp_sizes", "n0", "n1", "ngrpData", "ntablesinfo", "spectral",
                         "sd_par")
grid_param = grid_param[which(grid_param$ntablesinfo<=grid_param$Tables),]

dim(grid_param)

ncores = 4

# ---- Loop for graphs : ----

for(i in 1:nrow(grid_param)){
  
  n_tables = grid_param[i, 1]
  n_grp = grid_param[i, 2]
  grp_size = grid_param[i, 3]
  n_ind = n_grp*grp_size
  n0 = grid_param[i, 4]
  n1 = grid_param[i, 5]
  n_grp_per_data = grid_param[i, 6]
  n_tab_info = grid_param[i,7]
  reference <- rep(1:n_grp, each = grp_size) 
  spectral = grid_param[i,8]
  sd_par = grid_param[i,9]
  
  res_scen = lapply(sep_par, FUN = function(sep){function_sim(n_tables, n0, n1, n_grp, grp_size, sep, n_grp_per_data,
                                                              n_eval = 2, reference = reference, ntab_info = n_tab_info,
                                                              spectral = spectral, ncores = ncores, sd_par = sd_par)})
  
  res_plot = do.call("rbind", c(lapply(res_scen, FUN = function(res) res$gg_tab1), lapply(res_scen, FUN = function(res) res$gg_tab2)))
  
  res_plot$Spectral = "Non"
  res_plot$Spectral[grep("sp", res_plot$Method)] = "Spectral"

  nb_cols = length(unique(res_plot$Method))
  
  p_all = ggplot(res_plot,  aes(N, NID, colour = Method)) +
    geom_line() +
    geom_point(aes(shape = Method), size = 1.5) +
    scale_shape_manual(values = c(1:nb_cols)) +
    facet_grid(Scenario ~ Separability) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    guides(shape = guide_legend(nrow = 1, ncol = nb_cols), colour = NULL) +
    xlim(c(0, 25)) +
    ylim(c(0,1)) +
    xlab("Number of groups") +
    ggtitle(paste0("n0=", n0, "_n1=", n1, "_ngroups=", n_grp, "_grp size=", grp_size,
                   "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data, "_n_info=" , n_tab_info,
                   "_spectral=", spectral,"_sdpar=", sd_par))
  
  print(p_all)
  ggsave(paste0("n0=", n0, "_n1=", n1, "_ngroups=", n_grp, "_grp size=", grp_size,
                "_ntab=",n_tables, "_ngrp_per_data=", n_grp_per_data, "_n_info=" , n_tab_info,
                "_spectral=", spectral, "_sdpar=", sd_par,
                ".jpeg"), height = 20, width = 30, units = c("cm"))  
}