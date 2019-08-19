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

setwd("/home/hulot/Documents/these_documents/articles/Tree_aggregation_method_1d/Simulations/")
source("util.R")

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

# ---- Separability : 0.5 ----

## SIMULATION PARAMETERS
separability <- 0.5

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
# 
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
# 
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
# 
# FAVORABLE À MERGETREES
# 
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE) +
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

# Affichage

grid.arrange(p1, p2, nrow = 2)
ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")

p1.1 <- res_plot1
p1.1$Separability = 0.5
p1.1$Scenario = 1

p1.2 <- res_plot2
p1.2$Separability = 0.5
p1.2$Scenario = 2

# ---- Separability : 1 ----

## SIMULATION PARAMETERS

separability <- 1

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
# 
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
# 
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
# 
# FAVORABLE À MERGETREES
# 
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE) +
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

grid.arrange(p1, p2, nrow = 2)

ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")


p2.1 <- res_plot1
p2.1$Separability = 1
p2.1$Scenario = 1

p2.2 <- res_plot2
p2.2$Separability = 1
p2.2$Scenario = 2

# ---- Separability : 5 ----

## SIMULATION PARAMETERS
separability <- 5

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
#
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
#
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>%
  group_by(method, nb_grp) %>%
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
#
# FAVORABLE À MERGETREES
#
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>%
  group_by(method, nb_grp) %>%
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
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
  set.seed(sim_label*1000)
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot1 <- Reduce("rbind", res_scenario1) %>%
  group_by(method, nb_grp) %>%
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
#
# FAVORABLE À MERGETREES
#
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  set.seed(sim_label*1000)
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size)
  oneRun(sim_label, dataSets, reference, k = k_svd)
}, mc.cores = ncores)

res_plot2 <- Reduce("rbind", res_scenario2) %>%
  group_by(method, nb_grp) %>%
  summarize(NID_sd = sd(NID), NID = mean(NID)) %>% 
  filter(method %in% c("One", "MTW", "DC", "AD", "rScMTW", "rScAD", "rScDC"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) +
  geom_line() + geom_point(aes(shape = method)) + xlim(0, x_lim) + scale_color_viridis(discrete = TRUE)+
  ggtitle(paste0("n_grp=", n_grp, " ; grp_size=", grp_size, " ; sep=", separability, " ; n_grp_per_data=", n_grp_per_data, " ; n0=",  n0, " ; n1=", n1, " ; n_sim=", n_sim))

ggsave(p1, filename = paste0("../Simulations/p1_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")
ggsave(p2, filename = paste0("../Simulations/p2_", separability,".pdf"), width = width_fig, height = height_fig, units = "cm")

# Affichage

grid.arrange(p1, p2, nrow = 2)

ggsave(grid.arrange(p1, p2, nrow = 2), filename = paste0("../Simulations/", separability,".pdf"), width = width_fig, height = 2*height_fig, units = "cm")

# ---- GGPLOT : ----

# p3.1 <- res_plot1
# p3.1$Separability = 10
# p3.1$Scenario = 1
# 
# p3.2 <- res_plot2
# p3.2$Separability = 10
# p3.2$Scenario = 2
# 
# save(p1.1, p1.2, p2.1, p2.2, p3.1, p3.2, file = "plots_data.RData")
setwd("/home/hulot/Documents/these_documents/articles/Tree_aggregation_method_1d/Simulations/")
load("plots_data.RData")


res_plot = do.call("rbind", list(p1.1, p1.2, p2.1, p2.2, p3.1, p3.2))
colnames(res_plot)[which(colnames(res_plot)=="method")] = "Method"

# Changement noms de methodes pour coller avec le reste de l'article
Method_vect = as.character(res_plot$Method)
# Method_vect[which(ggplot_all$Method == "ACuni")] = "AD" 
Method_vect[which(res_plot$Method == "MTW")] = "MT"
Method_vect[which(res_plot$Method == "rScDC")] = "SpDC" 
Method_vect[which(res_plot$Method == "rScAD")] = "SpAD" 
Method_vect[which(res_plot$Method == "rScMTW")] = "SpMT" 

res_plot$Method = factor(Method_vect)

p_all = ggplot(res_plot,  aes(nb_grp, NID, color = Method)) +
  geom_line() + geom_point(aes(shape = Method)) +
  facet_grid(Scenario ~ Separability) +
  scale_color_viridis(discrete = TRUE) +
  scale_x_continuous(trans = 'log', breaks = c(1,2,4,10,25,50,100)) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1)) +
  xlab("Number of groups")
p_all
# ggsave(p_all, filename = paste0("Simulations_plots.pdf"), width = 35, height = 20, units = "cm")
# ggsave(p_all, filename = paste0("Simulations_plots.pdf"), width = 18, height = 10, units = "cm")
ggsave(plot = p_all, filename = "Simulations_plots.pdf", device = "pdf", width = 25, height = 15, units = "cm")

# p_all.1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
#   geom_line() + geom_point(aes(shape = method)) + 
#   scale_color_viridis(discrete = TRUE) +
#   scale_x_continuous(trans = "log", breaks = c(1, 2 , 4, 10, 25, 50, 100))
# 
# p3.2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
#   geom_line() + geom_point(aes(shape = method)) +  
#   scale_color_viridis(discrete = TRUE) +
#   scale_x_continuous(trans = "log", breaks = c(1, 2 , 4, 10, 25, 50, 100))




