#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information
#' @param n_grp  nombre de groupe
#' @param grp_size  taille des groupes
#' @param separability coefficient de séparabilité des groupes
scenario1 <- function(n0, n1, n_grp, grp_size, separability){
  
  get_one_data <- function(i) {
    noise <- rnorm(grp_size * n_grp, sd = 1)
    mu <- rep(i * (1:n_grp), each = grp_size)
    scale(mu + noise, TRUE, FALSE)
  }
  
  grp_levels <- rep(c(0, separability), c(n0, n1))
  
  do.call("cbind", lapply(grp_levels, get_one_data))
}

#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information
#' @param n_grp  nombre de groupe
#' @param n_grp_per_data  nombre de groupes propres à chaque dataset
#' @param grp_size  taille des groupes
#' @param separability coefficient de séparabilité des groupes
scenario2 <- function(n0, n1, n_grp, n_grp_per_data, grp_size, separability){
  
  get_one_data <- function(separability_) {
    
    noise <- rnorm(grp_size * n_grp, sd = 1)
    
    ## draw n.dif group per feature
    grp_means <- sample(1:n_grp, n_grp_per_data)
    
    ## how many times per feature (at least one)
    grp_means_sizes <- as.vector(rmultinom(1, size = n_grp - n_grp_per_data, prob = rep(1/n_grp_per_data, n_grp_per_data)) + 1)
    
    ## group means in all data sets    
    means <- sample(rep(grp_means, grp_means_sizes))
    
    mu <- rep(separability_ * means, each = grp_size)
    scale(mu + noise, TRUE, FALSE)
  }
  
  grp_levels <- rep(c(0, separability), c(n0, n1))
  
  do.call("cbind", lapply(grp_levels, get_one_data))
}

# ------ FUNCTIONS ----- 

directClustering <- function(dataSets) {
  hclust(dist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
}

averagedClustering <- function(dataSets) {
  AD <- Reduce("+", lapply(dataSets, dist, method = "euclidean")) / length(dataSets)
  hclust(AD, method = "ward.D2")
}

mergeTreesWard1d <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    # hclust(dist(x, method = "euclidean"), method = "ward.D2")
        univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
} 
