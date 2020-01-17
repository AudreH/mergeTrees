

draw_data <- function(n, p, q) {
  lapply(1:q, FUN = function(x){matrix(rnorm(n*p),n,p)})
}


concat_mat = function(data){
  p = ncol(data[[1]])
  Dat = matrix(NA, ncol = p*length(data), nrow = nrow(data[[1]]))
  for(i in 0:(length(data)-1)){Dat[,(i*p+1):((i+1)*p)] = data[[i+1]]}
  return(Dat)
}

# ---- Servant pour des tables : ----

direct_clustering <- function(data) {
  hclust(dist(concat_mat(data)))
}

averaged_clustering <- function(data) {
  hclust(Reduce("+", lapply(data, dist)))
}

# ---- Avec decomp spectrale :  ----

spectral_DC = function(data, k_svd){
  Dat = concat_mat(data)
  
  dat_svd = rsvd(Dat, k = k_svd)
  hclust(dist(dat_svd$u %*% dat_svd$d))
}

spectral_AD = function(data, k_svd){
  Dat = concat_mat(data)
  
  dat_svd = rsvd(Dat, k = k_svd)
  hclust(Reduce("+", lapply(as.list(as.data.frame(dat_svd$u %*% dat_svd$d)), dist)))
}

spectral_MC = function(data, k_svd){
  Dat = concat_mat(data)
  
  dat_svd = rsvd(Dat, k = k_svd)
  mergeTrees::mergeTrees(lapply(as.list(as.data.frame(dat_svd$u %*% dat_svd$d)), univarclust::ward_1d))
}

# ---- Fonction de simu temps : ----
simu <- function(n, p, q,  n_eval = 3, k_rsvd = 10) {
  data <- draw_data(n, p, q)
  memory_check = Reduce("cbind", data)
  
  microbenchmark(
    DC = {M = matrix(NA, ncol = p*q, nrow = n); direct_clustering(data)},
    AC = averaged_clustering(data),
    MC = mergeTrees::mergeTrees(lapply(data, FUN = function(dat) hclust(dist(dat)))),
    spDC = spectral_DC(data, k_rsvd),
    spAD = spectral_AD(data, k_rsvd),
    spMT = spectral_MC(data, k_rsvd),
    times = n_eval
  )
}