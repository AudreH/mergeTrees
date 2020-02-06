rm(list = ls())

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(rsvd)
n_eval = 3

source("Fonctions.R")

# ---- Simulations nombre d'individus ----

lists_res_n = list()
n_possibilities = c(50,100,250,500, 1000, 5000, 10000, 12500, 15000, 20000, 25000, 30000, 40000, 50000, 100000, 500000, 10^6)
p_fix = 100
q_fix = 3
k_rsvd = 3

lists_res_n = list()

for(i in 1:length(n_possibilities)){
  print(n_possibilities[i])
  lists_res_n[[i]] = simu(n_possibilities[i], p_fix, q_fix, n_eval, k_rsvd)
  names(lists_res_n)[i] = n_possibilities[i]
  save(lists_res_n, file = paste0("Ind_Etape_", i,"_p_", p_fix, "_q_", q_fix, ".RData"))
}
