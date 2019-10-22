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
  
  lapply(grp_levels, get_one_data)
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
  
  lapply(grp_levels, get_one_data)
}

directClustering <- function(dataSets) {
  hclust(dist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
}

averagedClustering <- function(dataSets) {
  AD <- Reduce("+", lapply(dataSets, dist, method = "euclidean")) / length(dataSets)
  hclust(AD, method = "ward.D2")
}

mergeTreesWard1d <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    hclust(dist(x, method = "euclidean"), method = "ward.D2")
    #     univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
} 

## run approaches
oneRun <- function(sim_label, dataSets, reference, score = NID){
  
  ## Merge Tree
  time_MTW <- system.time(MTW_res <- mergeTreesWard1d(dataSets))[3]
  
  # Direct Clustering
  time_DC <- system.time(DC_res <- directClustering(dataSets))[3]
  
  # Average Distance
  time_AD <- system.time(AD_res <- averagedClustering(dataSets))[3]
  
  
  nb_ind <- length(dataSets[[1]])
  
  ## spectral versions with random SVD
  time_rSVD2 <- system.time( {rSVD <- rsvd(do.call("cbind", dataSets), k = 2) ; dataSets_rspectral2 <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))})[3]
  time_rSVD4 <- system.time( {rSVD <- rsvd(do.call("cbind", dataSets), k = 4) ; dataSets_rspectral4 <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))})[3]

  # # Spectral Direct Clustering
  # time_rScDC <- system.time(rScDC_res <- directClustering(dataSets_rspectral10))[3]
  
  # # Spectral Average Distance Clustering
  # time_rScAD <- system.time(rScAD_res <- averagedClustering(dataSets_rspectral10))[3]
  
  # Spectral WARD Merge trees
  time_rScMTW2 <- system.time(rScMTW_res2 <- mergeTreesWard1d(dataSets_rspectral2))[3] + time_rSVD2
  time_rScMTW4 <- system.time(rScMTW_res4 <- mergeTreesWard1d(dataSets_rspectral4))[3] + time_rSVD4

  # NID cutree
  NID_MTW      <- apply(cutree(MTW_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_DC       <- apply(cutree(DC_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_AD       <- apply(cutree(AD_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_rScMTW2  <- apply(cutree(rScMTW_res2, 1:nb_ind), 2, score, c2 = reference)
  NID_rScMTW4  <- apply(cutree(rScMTW_res4, 1:nb_ind), 2, score, c2 = reference)

  do.call(rbind, 
          list(
            MTW  =     data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="MT"      , NID=NID_MTW      , time = time_MTW),
            DC  =      data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="DC"      , NID=NID_DC       , time = time_DC),
            AD  =      data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="AD"      , NID=NID_AD       , time = time_AD),
            rScMTW2 =  data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SpMT_2"  , NID=NID_rScMTW2  , time = time_rScMTW2),
            rScMTW4 =  data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SpMT_4"  , NID=NID_rScMTW4  , time = time_rScMTW4)
          ))
  
}
