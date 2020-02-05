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
  
  noise <- rnorm(grp_size * n_grp, sd = 1)
  
  ## draw n.dif group per feature
  grp_means <- sample(1:n_grp, n_grp_per_data)
  
  ## how many times per feature (at least one)
  grp_means_sizes <- as.vector(rmultinom(1, size = n_grp - n_grp_per_data, prob = rep(1/n_grp_per_data, n_grp_per_data)) + 1)
  
  ## group means in all data sets    
  means <- sample(rep(grp_means, grp_means_sizes))
  get_one_data <- function(separability_) {
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

# ---- Simu donnees groupes : Fonction de run sur data ----

Run_func =  function(dataSets, n_grp, grp_size, reference, spectral = 2){
  dataSets = lapply(dataSets, FUN = function(dat) dat/svd(dat, nu = 0, nv = 0)$d[1]) # sait-on jamais.
  
  hc_methods = list()
  hc_methods$AD = averagedClustering(dataSets)
  hc_methods$DC = directClustering(dataSets)
  hc_methods$MC = mergeTrees(lapply(dataSets, FUN = function(dat) hclust(dist(dat, method = "euclidean"), method = "ward.D2")))
  
  svd_dat = svd(do.call("cbind", dataSets), nv = 0)
  
  k_svd = spectral
  dat_univar = as.list(data.frame(svd_dat$u[,1:k_svd] %*% diag(svd_dat$d[1:k_svd])))
  
  if(spectral == 2){
    hc_methods$spAD_2 = averagedClustering(dat_univar)
    hc_methods$spMC_2 = mergeTreesWard1d(dat_univar)
    hc_methods$spDC_2 = directClustering(dat_univar)
  }
  if(spectral == 4){
    hc_methods$spAD_4 = averagedClustering(dat_univar)
    hc_methods$spMC_4 = mergeTreesWard1d(dat_univar)
    hc_methods$spDC_4 = directClustering(dat_univar)
  }
  
  res_NID = lapply(hc_methods, FUN = function(hc){
    apply(cutree(hc, k = 1:(n_grp*grp_size)), 2, NID, c2 = reference)
  })

  tab_res_NID = do.call("rbind", res_NID)
  
  return(tab_res_NID)
}


# ---- Fonction de simulation/resultats : ----

function_sim = function(n_tables, n0, n1, n_grp, grp_size, separability, n_grp_per_data, 
                        n_eval = 10, reference,  ntab_info, ncores = 2, spectral = 2, sd_par = 1){
  
  # --- Scenario1 : ----
  res_sim1 = mclapply(1:n_eval, FUN = function(ev){
    dataSets1 = vector(mode = "list", length = n_tables)
    dataSets1[1:ntab_info] = lapply(1:ntab_info, FUN = function(i){set.seed(i*ev*100); scenario1(n0, n1, n_grp, grp_size, separability)}) 
    if(ntab_info!=n_tables){
      dataSets1[(ntab_info+1):n_tables] = lapply((ntab_info+1):n_tables, FUN = function(i){set.seed(i*ev*100); matrix(rnorm(grp_size * n_grp * (n0+n1), sd = sd_par), 
                                                                                                                      ncol = (n0+n1), nrow = grp_size * n_grp)})
    }
    Run_func(dataSets1, n_grp, grp_size, reference, spectral = spectral)
  }, mc.cores = ncores)
  
  res_tab1 = do.call("cbind", res_sim1)
  gg_tab1 = melt(res_tab1)
  colnames(gg_tab1) = c("Method", "N", "NID")
  
  gg_tab1 = gg_tab1 %>% group_by(Method, N) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  gg_tab1$Separability = separability
  gg_tab1$Scenario = "Scenario 1"
  
  p1 <- ggplot(gg_tab1, aes(N, NID, color = Method)) + 
    geom_line() + geom_point(aes(shape = Method)) + ggtitle(paste0("Scenario1 - sep = ", separability))
  
  # --- Scenario2 : ----
  
  res_sim2 = mclapply(1:n_eval, FUN = function(ev){
    dataSets2 = vector(mode = "list", length = n_tables)
    dataSets2[1:ntab_info] = lapply(1:ntab_info, FUN = function(i){set.seed(i*ev*100); scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability) }) 
    if(ntab_info!=n_tables){
      dataSets2[(ntab_info+1):n_tables] = lapply((ntab_info+1):n_tables, FUN = function(i){set.seed(i*ev*100); matrix(rnorm(grp_size * n_grp * (n0+n1), sd = sd_par), 
                                                                                                                      ncol = (n0+n1), nrow = grp_size * n_grp)})
    }
    Run_func(dataSets2, n_grp, grp_size, reference, spectral = spectral)
  }, mc.cores = ncores)
  
  res_tab2 = do.call("cbind", res_sim2)
  gg_tab2 = melt(res_tab2)
  colnames(gg_tab2) = c("Method", "N", "NID")
  
  gg_tab2 = gg_tab2 %>% group_by(Method, N) %>% 
    summarize(NID_sd = sd(NID), NID = mean(NID)) 
  gg_tab2$Separability = separability
  gg_tab2$Scenario = "Scenario 2"
  
  p2 <- ggplot(gg_tab2, aes(N, NID, color = Method)) + 
    geom_line() + geom_point(aes(shape = Method)) + ggtitle(paste0("Scenario2 - sep = ", separability))
  
  
  return(list(gg_tab1 = gg_tab1, res_sim1 = res_sim1, p1 = p1,
              gg_tab2 = gg_tab2, res_sim2 = res_sim2, p2 = p2))
}

