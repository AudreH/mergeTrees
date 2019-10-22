rm(list = ls())
library(mergeTrees)
library(univarclust)
library(aricode) 
library(rsvd)
library(parallel)
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
theme_set(theme_bw())
source("util.R")

width_fig = 20
height_fig = 10

## separability
## 10:  everyboby is working

## SIMULATION PARAMETERS (COMMON SETTINGS)

ncores = 10
n_grp <- 10
grp_size <- 50
n0 <- 500
n1 <- 50
n_grp_per_data <- 2
n_sim <- 10
x_lim = 100

# ncores = 10
# n_grp <- 4
# grp_size <- 25
# n0 <- 1000
# n1 <- 100
# n_grp_per_data <- 3
# n_sim <- 10


# ---- Essai : ----

separability = c(seq(0.1,1,by = 0.1), c(seq(1,10, by = 1)))
# separability = 1

separability = c(.2, 0.4, 1)

Res_separability = lapply(separability, FUN = function(sep_param){
  
  res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
    set.seed(sim_label*1000)
    dataSets  <- scenario1(n0, n1, n_grp, grp_size, sep_param)
    reference <- rep(1:n_grp, each = grp_size) 
    oneRun(sim_label, dataSets, reference)
  }, mc.cores = ncores)
  
  
  # res_scenario1 <- lapply(1:n_sim, function(sim_label) {
  #   set.seed(sim_label*1000)
  #   dataSets  <- scenario1(n0, n1, n_grp, grp_size, sep_param)
  #   reference <- rep(1:n_grp, each = grp_size)
  #   oneRun(sim_label, dataSets, reference)
  # })

  
  res_plot1 <- Reduce("rbind", res_scenario1) %>% 
    group_by(method, nb_grp) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  
  p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
    geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
    ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", sep_param, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))
  
  res_plot1$Separability = sep_param
  res_plot1$Scenario = 1
  
  res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
    set.seed(sim_label*1000)
    dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, sep_param)
    reference <- rep(1:n_grp, each = grp_size) 
    oneRun(sim_label, dataSets, reference)
  }, mc.cores = ncores)
  
  res_plot2 <- Reduce("rbind", res_scenario2) %>% 
    group_by(method, nb_grp) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  
  p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
    geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE) +
    ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", sep_param, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))
  
  ggsave(p1, filename = paste0("p1_", sep_param,".pdf"), width = width_fig, height = height_fig, units = "cm")
  ggsave(p2, filename = paste0("p2_", sep_param,".pdf"), width = width_fig, height = height_fig, units = "cm")
  
  res_plot2$Separability = sep_param
  res_plot2$Scenario = 2
  
  return(list(res_plot1  = res_plot1, res_plot2 = res_plot2))
})

# save(Res_separability, file = "Res_separability.RData")
# load("Res_separability.RData")

# ---- Res plot All ----

res_plot = do.call("rbind", c(lapply(Res_separability, FUN = function(res) res$res_plot1),
                                 lapply(Res_separability, FUN = function(res) res$res_plot2)))
colnames(res_plot)[which(colnames(res_plot)=="method")] = "Method"

# Changement noms de methodes pour coller avec le reste de l'article
Method_vect = as.character(res_plot$Method)

res_plot$Method = factor(Method_vect)

p_all = ggplot(res_plot,  aes(nb_grp, NID, color = Method)) +
  geom_line() + geom_point(aes(shape = Method)) +
  facet_grid(Scenario ~ Separability) +
  scale_color_viridis(discrete = TRUE) +
  # scale_x_continuous(trans = 'log', breaks = c(1,2,4,10,25,50,100)) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1)) + xlim(c(0,50)) +
  xlab("Number of groups")
p_all
ggsave(plot = p_all, filename = "Simulations_plots_all.pdf", device = "pdf", width = 100, height = 15, units = "cm")


# ---- Res plot Sauf SpAD et spDC ----
# 
# res_plot = do.call("rbind", c(lapply(Res_separability, FUN = function(res) res$res_plot1),
#                               lapply(Res_separability, FUN = function(res) res$res_plot2)))
# colnames(res_plot)[which(colnames(res_plot)=="method")] = "Method"
# 
# # Changement noms de methodes pour coller avec le reste de l'article
# Method_vect = as.character(res_plot$Method)
# res_plot2 = res_plot[-which(res_plot$Method%in%c("SpAD", "SpDC")),]
# 
# res_plot2$Method = factor(Method_vect)
# 
# p_all = ggplot(res_plot2,  aes(nb_grp, NID, color = Method)) +
#   geom_line() + geom_point(aes(shape = Method)) +
#   facet_grid(Scenario ~ Separability) +
#   scale_color_viridis(discrete = TRUE) +
#   scale_x_continuous(trans = 'log', breaks = c(1,2,4,10,25,50,100)) +
#   theme(legend.position="bottom", legend.box = "horizontal") +
#   guides(color = guide_legend(nrow = 1)) +
#   xlab("Number of groups")
# p_all
# ggsave(plot = p_all, filename = "Simulations_plots_all2.pdf", device = "pdf", width = 100, height = 15, units = "cm")
# 
