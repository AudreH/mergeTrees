rm(list = ls())

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(rsvd)
n_eval = 3

source("Fonctions.R")

# ---- Simulations nombre de tables ----

lists_res_n = list()
q_possibilities = c(50,100,250,500, 1000, 5000, 10000, 12500, 15000, 20000, 25000, 30000, 40000, 50000, 100000, 500000, 10^6)
p_fix = 100
n_fix = 100
k_rsvd = 3

lists_res_q = list()

for(i in 1:length(q_possibilities)){
  print(q_possibilities[i])
  lists_res_q[[i]] = simu(n_fix, p_fix, q_possibilities[[i]], n_eval, k_rsvd)
  names(lists_res_q)[i] = q_possibilities[i]
  save(lists_res_q, file = paste0("Tables_Etape_", i,"_n_", n_fix, "_p_", p_fix, ".RData"))
}

