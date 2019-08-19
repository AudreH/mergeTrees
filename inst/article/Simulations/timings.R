rm(list = ls())
library(Rmergetrees)
library(fusedanova)
library(aricode) 
library(rsvd)
library(parallel)
library(dplyr)
library(ggplot2)
library(viridis)
theme_set(theme_bw())
source("util.R")

## separability
## 10:  everyboby is working

## SIMULATION PARAMETERS
separability <- 2

## 
n_grp <- 20
grp_size <- 50

##
n0 <- 30
n1 <- 50
n_grp_per_data <- 2

n_sim <- 5

res_growing_features <- c()
for (factor in c(1:10)) {
  
  res_growing_features <- 
    rbind(res_growing_features, 
          Reduce("rbind",
          mclapply(1:n_sim, function(sim_label) {

  dataSets <- scenario1(n0 * factor, n1 * factor, n_grp, grp_size, separability)
  time_rSVD <- system.time(
    {
      rSVD <- rsvd(do.call("cbind", dataSets), k = 10) 
      dataSets_rspectral <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))  
    }
  )[3]

  # Spectral Merge trees
  time_SMT <- system.time(SMT_res <- mergeTreesWard1d(dataSets_rspectral))[3] + time_rSVD

  # Spectral Merge trees
  time_MT <- system.time(MT_res <- mergeTreesWard1d(dataSets))[3]

  
  data.frame(time = c(time_MT, time_SMT, time_rSVD), method = c("mergeTree", "specMergeTree", "SVD"), sample_size = length(dataSets[[1]]), features = length(dataSets), sim_label = sim_label)

  
  }, mc.cores = 10)))

}


res_growing_samplesize <- c()
for (factor in c(1:10)) {
  
  res_growing_samplesize <- 
    rbind(res_growing_samplesize, 
          Reduce("rbind",
          mclapply(1:n_sim, function(sim_label) {

  dataSets <- scenario1(n0, n1, n_grp , grp_size * factor, separability)
  time_rSVD <- system.time(
    {
      rSVD <- rsvd(do.call("cbind", dataSets), k = 10) 
      dataSets_rspectral <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))  
    }
  )[3]

  # Spectral Merge trees
  time_SMT <- system.time(SMT_res <- mergeTreesWard1d(dataSets_rspectral))[3] + time_rSVD

  # Spectral Merge trees
  time_MT <- system.time(MT_res <- mergeTreesWard1d(dataSets))[3]

  
  data.frame(time = c(time_MT, time_SMT, time_rSVD), method = c("mergeTree", "specMergeTree", "SVD"), sample_size = length(dataSets[[1]]), features = length(dataSets), sim_label = sim_label)

  
  }, mc.cores = 10)))

}

p_sample <- ggplot(res_growing_samplesize, aes(x = sample_size, y = time, shape = method)) + geom_point()
p_feature <-ggplot(res_growing_features  , aes(x = features, y = time, shape = method)) + geom_point()
