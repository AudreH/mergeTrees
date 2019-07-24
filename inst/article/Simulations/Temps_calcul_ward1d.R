rm(list = ls())
# Code repris de celui de Julien du 26/01/19 (mail)
# Investigation temps de calcul étrange de ward1D en comparaison à celui des autres méthodes

library(mergeTrees)
library(univarclust)
library(microbenchmark)
library(ggplot2)
library(rsvd)
library(parallel)
library(viridis)
library(profvis)
library(Rcpp)



# ---- Fonctions : ----

draw_data <- function(n, p) {
  matrix(rnorm(n*p),n,p) # Nombre d'individus, nombre de variables, c'est tout. 
}

direct_clustering <- function(data) {
  hclust(dist(Reduce("cbind", data)))
}

averaged_clustering <- function(data, hc.list, standardize = FALSE) {
  n <- unique(sapply(hc.list, function(hc) length(hc$order)))
  stopifnot(is.integer(n))
  p <- length(hc.list)
  if (standardize) {
    hc.list <- lapply(hc.list, function(hc) {
      hc$height <- hc$height/max(hc$height)
      hc
    })
  }
  labels <- Reduce(intersect, lapply(hc.list, function(hc) hc$labels))
  if (!is.null(labels)) {
    stopifnot(length(labels) == n)
  }
  hclust(Reduce("+", lapply(data, dist)))
}

mergeTreesWard <- function(data) {
  mergeTrees::mergeTrees(lapply(data, FUN = function(x) {univarclust::ward_1d(x)}))} 


# Deuxieme version beaucoup plus longue que la premiere. Reste a savoir ce qui prend du temps dans le mergeTrees
# mergeTreesWard2 <- function(data){
#   mergeTrees(lapply(data, FUN = function(x) hclust(dist(x, method = "euclidean"), method = "ward.D2")))
# }


mergeTreesWard4 = function(data){
  lapply(data, FUN = function(x) {univarclust::ward_1d(x)})
}

# ---- Settings data test : ----

n_ = 500
p_ = 2500

data_sets = draw_data(n_, p_)
data_univar = as.list(as.data.frame(data_sets))
length(data_univar)

n_eval = 10

hc_list = lapply(data_univar, FUN = function(dat) hclust(dist(dat), method = "ward.D2"))

# ---- Tests : profvis ----

profvis({		
  ACuni = averaged_clustering(data_univar, hc.list = hc_list)
})

profvis({		
  MTuni = mergeTreesWard(data_univar)
})

profvis({		
  MTuni = mergeTreesWard4(data_univar)
})

# --- Tests : microbenchmark : ----

test_m = microbenchmark(
  ACuni = averaged_clustering(data_univar, hc.list = hc_list),
  MTuni = mergeTreesWard(data_univar),
  MTuni3 = mergeTreesWard4(data_univar),
  times = n_eval
)

test_m
