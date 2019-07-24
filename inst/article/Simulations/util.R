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

# mergeTreesClustering <- function(dataSets) {
#   hc_list <- lapply(dataSets, FUN = function(x) {
#     hclust(dist(x, method = "euclidean"), method = "ward.D2")
#   })
#   mergeTrees::mergeTrees(hc_list)
# } 

mergeTreesWard1d <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    hclust(dist(x, method = "euclidean"), method = "ward.D2")
    #     univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
} 

# oneDimClustering <- function(dataSet) {
#   hclust(dist(dataSet, method = "euclidean"), method = "ward.D2")
# }

## run approaches
oneRun <- function(sim_label, dataSets, reference, score = NID, k = 10){
  
  # ## Merge Tree
  # time_MT <- system.time(MT_res <- mergeTreesClustering(dataSets))[3]

  ## Merge Tree
  time_MTW <- system.time(MTW_res <- mergeTreesWard1d(dataSets))[3]
    
  ## One dimension
  # time_OD <- system.time(OD_res <- oneDimClustering(dataSets[[length(dataSets)]]))[3]
  
  # Direct Clustering
  time_DC <- system.time(DC_res <- directClustering(dataSets))[3]

  # Average Distance
  time_AD <- system.time(AD_res <- averagedClustering(dataSets))[3]

  # # Fused-ANOVA Merge trees
  # time_FA <- system.time(FA_res <- fusedanova(Reduce("cbind", dataSets), gamma = rep(gamma, length(dataSets))))[3]

  # ## Some spectral versions
  # time_SVD <- system.time(
  #   {
  #     SVD <- svd(do.call("cbind", dataSets))
  #     dataSets_spectral <- as.list(as.data.frame(SVD$u %*% diag(SVD$d)))
  #   }
  # )[3]

  # # Spectral Direct Clustering
  # time_SDC <- system.time(SDC_res <- directClustering(dataSets_spectral))[3] + time_SVD

  # # Spectral Merge trees
  # time_SMT <- system.time(SMT_res <- mergeTreesClustering(dataSets_spectral))[3] + time_SVD

  # # Spectral Ward Merge trees
  # time_SMTW <- system.time(SMTW_res <- mergeTreesWard1d(dataSets_spectral))[3] + time_SVD
  
  # # Spectral Average distance Clustering
  # time_SAD <- system.time(SAD_res <- averagedClustering(dataSets_spectral))[3] + time_SVD
      
  # # Spectral Fused-ANOVA Merge trees
  # time_SFA <- system.time(SFA_res <- fusedanova(Reduce("cbind", dataSets_spectral), gamma = rep(gamma, length(dataSets))))[3] + time_SVD
    
  nb_ind <- length(dataSets[[1]])

  ## spectral versions with random SVD
  time_rSVD <- system.time(
    {
      rSVD <- rsvd(do.call("cbind", dataSets), k = k)
      dataSets_rspectral <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))
    }
  )[3]
  # 
  # # Spectral Direct Clustering
  time_rScDC <- system.time(rScDC_res <- directClustering(dataSets_rspectral))[3]

  # # Spectral Average Distance Clustering
  time_rScAD <- system.time(rScAD_res <- averagedClustering(dataSets_rspectral))[3]
  
  # Spectral Merge trees
  # time_rSMT <- system.time(rSMT_res <- mergeTreesClustering(dataSets_rspectral))[3] + time_rSVD
  
  # Spectral WARD Merge trees
  time_rScMTW <- system.time(rScMTW_res <- mergeTreesWard1d(dataSets_rspectral))[3] + time_rSVD
  #   
  # # Spectral Fused-ANOVA Merge trees
  # time_rSFA <- system.time(rSFA_res <- fusedanova(Reduce("cbind", dataSets_rspectral), gamma = rep(gamma, k)))[3]
  # NID_1    <- apply(cutree(OD_res  , 1:nb_ind), 2, score, c2 = reference)
  # NID_MT   <- apply(cutree(MT_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_MTW  <- apply(cutree(MTW_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_DC   <- apply(cutree(DC_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_AD   <- apply(cutree(AD_res  , 1:nb_ind), 2, score, c2 = reference)
  # NID_FA   <- apply(cutree(FA_res  , 1:nb_ind), 2, score, c2 = reference)
  # NID_SMT  <- apply(cutree(SMT_res , 1:nb_ind), 2, score, c2 = reference)
  # NID_SMTW <- apply(cutree(SMTW_res , 1:nb_ind), 2, score, c2 = reference)
  # NID_SDC  <- apply(cutree(SDC_res , 1:nb_ind), 2, score, c2 = reference)
  # NID_SAD   <- apply(cutree(AD_res  , 1:nb_ind), 2, score, c2 = reference)
  # NID_SFA  <- apply(cutree(SFA_res , 1:nb_ind), 2, score, c2 = reference)
  # NID_rScMT <- apply(cutree(rSMT_res, 1:nb_ind), 2, score, c2 = reference)
  NID_rScMTW <- apply(cutree(rScMTW_res, 1:nb_ind), 2, score, c2 = reference)
  NID_rScDC <- apply(cutree(rScDC_res, 1:nb_ind), 2, score, c2 = reference)
  NID_rScAD <- apply(cutree(rScAD_res, 1:nb_ind), 2, score, c2 = reference)
  # NID_rSFA <- apply(cutree(rSFA_res, 1:nb_ind), 2, score, c2 = reference)
  
  do.call(rbind, 
          list(
            # OD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="One"  , NID=NID_1    , time = time_OD),
               # MT  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="MT"   , NID=NID_MT   , time = time_MT),
               MTW  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="MTW"   , NID=NID_MTW   , time = time_MTW),
               DC  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="DC"   , NID=NID_DC   , time = time_DC),
               AD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="AD"   , NID=NID_AD   , time = time_AD),
               # FA  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="FA"   , NID=NID_FA   , time = time_FA),
               # SMT = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SMT"  , NID=NID_SMT  , time = time_SMT),
               # SDC = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SDC"  , NID=NID_SDC  , time = time_SDC),
               # SAD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SAD"   , NID=NID_SAD   , time = time_SAD)
               # SFA = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SFA"  , NID=NID_SFA  , time = time_SFA),
               # rSMT = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="rSMT" , NID=NID_rSMT, time = time_rSMT),
               rScDC = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="rScDC" , NID=NID_rScDC, time = time_rScDC),
               # rSFA = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="rSFA" , NID=NID_rSFA, time = time_rSFA),
               rScMTW = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="rScMTW"  , NID=NID_rScMTW  , time = time_rScMTW),
               rScAD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="rScAD"   , NID=NID_rScAD   , time = time_rScAD)
          ))

}
