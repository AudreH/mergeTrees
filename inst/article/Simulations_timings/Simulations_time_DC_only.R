rm(list = ls())
library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(rsvd)
n_eval = 3

draw_data <- function(n, p) {
  matrix(rnorm(n*p),n,p) # Nombre d'individus, nombre de variables, c'est tout. 
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
    DC = {M = matrix(NA, ncol = p, nrow = n); direct_clustering(data_univar)},
    SDC = {rSVD <- rsvd(Reduce("cbind", data_univar), k = k_rsvd); 
           hclust(dist(as.data.frame(rSVD$u %*% diag(rSVD$d)), method = "euclidean"), method = "ward.D2")},
    times = n_eval
  )
}

# ---- Simulations nombre d'individus ----

lists_res_n = list()
n_possibilities = c(50,100,250,500, 1000, 5000, 10000, 12500, 15000, 20000, 25000, 30000, 40000, 50000, 100000, 500000, 10^6)
p_fix = 100
k_rsvd = 3

lists_res_n = list()

for(i in 1:length(n_possibilities)){
  print(n_possibilities[i])
  lists_res_n[[i]] = simu(n_possibilities[i], p_fix, n_eval, k_rsvd)
  names(lists_res_n)[i] = n_possibilities[i]
  save(lists_res_n, file = paste0("Ind_Etape_", i,"_p_", p_fix, ".RData"))
}

lists_res_n = lapply(n_possibilities, FUN = function(n_var) simu(n_var, p_fix, n_eval, k_rsvd))
save(lists_res_n, file =  paste0("Timings_individuals_DC_p_", p_fix, ".RData"))
