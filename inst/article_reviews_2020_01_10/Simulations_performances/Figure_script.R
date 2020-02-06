rm(list = ls())

library(mergeTrees)
library(aricode)
library(svd)
library(reshape2)
library(tidyverse)  
library(RColorBrewer)
library(parallel)

source("util.R")

n_tables = 7
n_grp = 5
grp_size = 25
n_ind = n_grp*grp_size
n0 = 100
n1 = 50
n_grp_per_data = 2
n_tab_info = n_tables
reference <- rep(1:n_grp, each = grp_size) 
spectral = 2
sd_par = 1
ncores = 4

sep_par = c(0.2, 1, 2)

res_scen = lapply(sep_par, FUN = function(sep){function_sim(n_tables, n0, n1, n_grp, grp_size, sep, n_grp_per_data,
                                                            n_eval = 5, reference = reference, ntab_info = n_tab_info,
                                                            spectral = spectral, ncores = ncores, sd_par = sd_par)})

res_plot = do.call("rbind", c(lapply(res_scen, FUN = function(res) res$gg_tab1), lapply(res_scen, FUN = function(res) res$gg_tab2)))

res_plot$Spectral = "Non"
res_plot$Spectral[grep("sp", res_plot$Method)] = "Spectral"

library(RColorBrewer)
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
  xlab("Number of groups") 

print(p_all)
ggsave(plot = p_all, filename = "Figure3.pdf", device = "pdf", width = 25, height = 15, units = "cm")