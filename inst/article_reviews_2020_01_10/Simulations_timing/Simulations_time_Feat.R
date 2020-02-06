rm(list = ls())

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(rsvd)
n_eval = 3

source("Fonctions.R")

# ---- Simulations nombre de variables/arbres ----

lists_res_p = list()
p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000, 100000, 500000, 10^6)
n_fix = 100
q_fix = 3
k_rsvd = 3
lists_res_p = list()

for(i in 1:length(p_possibilities)){
  print(p_possibilities[i])
  lists_res_p[[i]] = simu(n_fix, p_possibilities[i], q_fix, n_eval, k_rsvd)
  names(lists_res_p)[i] = p_possibilities[i]
  save(lists_res_p, file = paste0("Features_Etape_", i,"_n_", n_fix, "_q_", q_fix,  ".RData"))
}