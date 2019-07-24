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
theme_set(theme_bw())
source("../Simulations/util.R")

width_fig = 20
height_fig = 10

## separability
## 10:  everyboby is working

## SIMULATION PARAMETERS (COMMON SETTINGS)

ncores = 2
n_grp <- 4
grp_size <- 25
n0 <- 1000
n1 <- 100
n_grp_per_data <- 3
n_sim <- 10
# k_svd = max(ceiling((n0+n1)/10), 2)
k_svd = 5
x_lim = 50 

scen_1 = scenario1(n0, n1, n_grp, grp_size, separability)
scen_2 = scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)

# ---- Separability : 0.5 ----

## SIMULATION PARAMETERS
separability <- 0.5

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
# 
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
# 
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
# 
# FAVORABLE À MERGETREES
# 
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE) +
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

# Affichage

grid.arrange(p1, p2, nrow = 2)
ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")



# ---- Separability : 1 ----

## SIMULATION PARAMETERS

separability <- 1

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
# 
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
# 
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
# 
# FAVORABLE À MERGETREES
# 
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE) +
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

grid.arrange(p1, p2, nrow = 2)

ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")


# ---- Separability : 5 ----

## SIMULATION PARAMETERS
separability <- 5

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
#
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
#
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>%
  group_by(method, nb_grp) %>%
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>%
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
#
# FAVORABLE À MERGETREES
#
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>%
  group_by(method, nb_grp) %>%
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>%
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

# Affichage

grid.arrange(p1, p2, nrow = 2)

ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")

# ---- Separability : 10 ----

## SIMULATION PARAMETERS
separability <- 10

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
#
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
#
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>%
  group_by(method, nb_grp) %>%
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>%
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
#
# FAVORABLE À MERGETREES
#
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>%
  group_by(method, nb_grp) %>%
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>%
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

# Affichage

grid.arrange(p1, p2, nrow = 2)

ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")

