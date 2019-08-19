rm(list = ls())
# Code repris de celui de Julien du 26/01/19 (mail)

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(ggplot2)
library(rsvd)
library(parallel)
library(viridis)

ncores = 3
n_eval = 10

draw_data <- function(n, p) {
   matrix(rnorm(n*p),n,p) # Nombre d'individus, nombre de variables, c'est tout. 
}

direct_clustering <- function(data) {
  hclust(dist(Reduce("cbind", data)))
}

averaged_clustering <- function(data) {
  hclust(Reduce("+", lapply(data, dist)))
}

# merged tree multivarie : pas interessant ici.
# merged_clustering <- function(data) {
#   mergeTrees::mergeTrees(
#     lapply(
#       lapply(data, dist),
#       hclust
#     )
#   )
# }

mergeTreesWard <- function(data) {
  mergeTrees::mergeTrees(lapply(data, FUN = function(x) {univarclust::ward_1d(x)}))} 

simu <- function(n, p,  n_eval = 10, prop_rsvd = 0.1) {
  data <- draw_data(n, p)
  data_univar = as.list(as.data.frame(data))
  memory_check = Reduce("cbind", data_univar)
  
  microbenchmark(
    
    DC = {M = matrix(NA, ncol = p, nrow = n); direct_clustering(data_univar)},
    SDC = {rSVD <- rsvd(Reduce("cbind", data_univar), k = prop_rsvd*length(data_univar));
    hclust(dist(as.data.frame(rSVD$u %*% diag(rSVD$d)), method = "euclidean"), method = "ward.D2")},
    
    # Methodes univariees
    # ACuni = averaged_clustering(data_univar),
    # MTuni = mergeTreesWard(data_univar),
    
    # Methods univariees spectrales
    # SACuni = {rSVD <- rsvd(Reduce("cbind", data_univar), k = prop_rsvd*length(data_univar));
    # averaged_clustering(as.list(as.data.frame(rSVD$u %*% diag(rSVD$d))))},
    # 
    # SMTuni = {rSVD <- rsvd(Reduce("cbind", data_univar), k = prop_rsvd*length(data_univar));
    # mergeTreesWard(as.list(as.data.frame(rSVD$u %*% diag(rSVD$d))))},
    
    times = n_eval
  )
}

# print(simu(100, 1000, 3))

# ---- Simulations nombre d'individus ----

lists_res_n = list()
# n_possibilities = c(50,100,250,500, 1000, 5000, 10000)
n_possibilities = c(10000)
p_fix = 1000
# q_fix = 3
# n_eval = 10
prop_rsvd = 0.25

lists_res_n = mclapply(n_possibilities, FUN = function(n_ind) simu(n_ind, p_fix, n_eval, prop_rsvd), mc.cores = ncores)
save(lists_res_n, file = "Timings_individuals_10000_DC.RData")
load("Timings_individuals2.RData")
names(lists_res_n) = n_possibilities

lists_res_n
df_list = lapply(lists_res_n, data.frame)
df_list = lapply(1:length(df_list), FUN = function(indice){
  df_list[[indice]]$n_ind= n_possibilities[[indice]]
  return(df_list[[indice]])
})

all.df_n_ind = do.call("rbind", df_list)
ggplot_n_ind = aggregate(all.df_n_ind$time/10^9, by =  list(all.df_n_ind$expr, all.df_n_ind$n_ind), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_n_ind = cbind(ggplot_n_ind[,1:2], ggplot_n_ind$x)
colnames(ggplot_n_ind) = c("Method", "N", "Mean", "SD")

ggplot_ind = ggplot(data = ggplot_n_ind, aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  # facet_grid(Dimensions ~ Proportions, scales = "free") +
  theme_bw() + xlab("Number of individuals") + ylab("Time (s)") +
  # geom_vline(xintercept = nb_groupes.) +
  scale_x_log10() +
  ggtitle("Simulations results - Time - Individuals per trees")
ggplot_ind

# ---- Simulations nombre de variables/arbres ----

lists_res_p = list()
p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000)
# p_possibilities = c(1000, 2500, 5000, 15000)
n_fix = 100
# q_fix = 3
# n_eval = 10
prop_rsvd = 0.25

lists_res_p = mclapply(p_possibilities, FUN = function(n_var) simu(n_fix, n_var, n_eval, prop_rsvd), mc.cores = ncores)
save(lists_res_p, file = "Timings_features.RData")
load("Timings_features.RData")
names(lists_res_p) = p_possibilities

# lists_res_q
df_list = lapply(lists_res_p, data.frame)
df_list = lapply(1:length(df_list), FUN = function(indice){
  df_list[[indice]]$n_features= p_possibilities[[indice]]
  return(df_list[[indice]])
})

all.df_n_features = do.call("rbind", df_list)
ggplot_n_features = aggregate(all.df_n_features$time/10^9, by =  list(all.df_n_features$expr, all.df_n_features$n_features), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_n_features = cbind(ggplot_n_features[,1:2], ggplot_n_features$x)
colnames(ggplot_n_features) = c("Method", "N", "Mean", "SD")

ggplot_var= ggplot(data = ggplot_n_features, aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  # facet_grid(Dimensions ~ Proportions, scales = "free") +
  theme_bw() + xlab("Number of features/trees") + ylab("Time (s)") +
  # geom_vline(xintercept = nb_groupes.) +
  scale_x_log10() + scale_color_viridis(discrete = TRUE) +
  ggtitle("Simulations results - Time - Features / Trees to aggregate")
ggplot_var

# --- All ggplot ----

load("Timings_features.RData")
load("Timings_individuals.RData")

n_possibilities = c(50,100,250,500, 1000, 5000, 10000)
p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000)
n_eval = 10

names(lists_res_n) = n_possibilities

df_list = lapply(lists_res_n, data.frame)
df_list = lapply(1:length(df_list), FUN = function(indice){
  df_list[[indice]]$n_ind= n_possibilities[[indice]]
  return(df_list[[indice]])
})

all.df_n_ind = do.call("rbind", df_list)
ggplot_n_ind = aggregate(all.df_n_ind$time/10^9, by =  list(all.df_n_ind$expr, all.df_n_ind$n_ind), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_n_ind = cbind(ggplot_n_ind[,1:2], ggplot_n_ind$x)
colnames(ggplot_n_ind) = c("Method", "N", "Mean", "SD")

names(lists_res_p) = p_possibilities
df_list = lapply(lists_res_p, data.frame)
df_list = lapply(1:length(df_list), FUN = function(indice){
  df_list[[indice]]$n_features= p_possibilities[[indice]]
  return(df_list[[indice]])
})

all.df_n_features = do.call("rbind", df_list)
ggplot_n_features = aggregate(all.df_n_features$time/10^9, by =  list(all.df_n_features$expr, all.df_n_features$n_features), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_n_features = cbind(ggplot_n_features[,1:2], ggplot_n_features$x)
colnames(ggplot_n_features) = c("Method", "N", "Mean", "SD")


# ggplot_n_data$Compare = rep("Datasets", nrow(ggplot_n_data))
ggplot_n_features$Compare = rep("Features", nrow(ggplot_n_features))
ggplot_n_ind$Compare = rep("Individuals", nrow(ggplot_n_ind))

# ggplot_all = do.call(rbind, list(ggplot_n_data, ggplot_n_features, ggplot_n_ind))
ggplot_all = do.call(rbind, list(ggplot_n_features, ggplot_n_ind))
ggplot_all$Spectral = "Non spectral"
ggplot_all$Spectral[grep("S", ggplot_all$Method)] = "Spectral"


gg = ggplot(data = ggplot_all, aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  facet_grid(Spectral ~ Compare, scales = "free") + scale_x_continuous(trans = "log")+
  theme_bw() + ylab("Time (s)") + xlab("") + scale_color_viridis(discrete = TRUE) +
  ggtitle("Simulations results - Time")
gg

ggsave(plot = gg, filename = "Simulations_results_time.pdf", device = "pdf", width = 2*15, height = 2*14, units = "cm")

gg = ggplot(data = ggplot_all[grep("Individuals", ggplot_all$Compare),], aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  facet_grid(Spectral ~ ., scales = "free") + scale_x_continuous(trans = "log", breaks = n_possibilities)+
  theme_bw() + ylab("Time (s)") + xlab("") + scale_color_viridis(discrete = TRUE) +
  ggtitle("Simulations results - Time - Individuals")
gg

ggsave(plot = gg, filename = "Simulations_ind_time.pdf", device = "pdf", width = 15, height = 14, units = "cm")

gg = ggplot(data = ggplot_all[grep("Features", ggplot_all$Compare),], aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  facet_grid(Spectral ~ ., scales = "free") + scale_x_continuous(trans = "log", breaks = p_possibilities)+
  theme_bw() + ylab("Time (s)") + xlab("") + scale_color_viridis(discrete = TRUE) +
  ggtitle("Simulations results - Time - Features")
gg

ggsave(plot = gg, filename = "Simulations_feat_time.pdf", device = "pdf", width = 15, height = 14, units = "cm")

