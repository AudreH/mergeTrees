rm(list = ls())

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(rsvd)
n_eval = 3

draw_data <- function(n, p) {
  matrix(rnorm(n*p),n,p) 
}

direct_clustering <- function(data) {
  hclust(dist(Reduce("cbind", data)))
}

averaged_clustering <- function(data) {
  hclust(Reduce("+", lapply(data, dist)))
}

mergeTreesWard <- function(data) {
  mergeTrees::mergeTrees(lapply(data, FUN = function(x) {univarclust::ward_1d(x)}))} 

simu <- function(n, p,  n_eval = 3, k_rsvd = 10) {
  data <- draw_data(n, p)
  data_univar = as.list(as.data.frame(data))
  memory_check = Reduce("cbind", data_univar)
  
  microbenchmark(
    ACuni = averaged_clustering(data_univar),
    SACuni = {rSVD <- rsvd(Reduce("cbind", data_univar), k = k_rsvd);
              averaged_clustering(as.list(as.data.frame(rSVD$u %*% diag(rSVD$d))))},
              
    times = n_eval
  )
}

# ---- Simulations nombre de variables/arbres ----

lists_res_p = list()
p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000, 100000, 500000, 10^6)
n_fix = 100
k_rsvd = 3
lists_res_p = list()

for(i in 1:length(p_possibilities)){
  print(p_possibilities[i])
  lists_res_p[[i]] = simu(n_fix, p_possibilities[i], n_eval, k_rsvd)
  names(lists_res_p)[i] = lists_res_p[i]
  save(lists_res_p, file = paste0("Features_Etape_", i, ".RData"))
}

save(lists_res_p, file = "Timings_features_AD.RData")
