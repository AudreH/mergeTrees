### CELLTYPE DIVISION PAR première valeur singulière


library(tidyverse)
library(devtools)
library(mergeTrees)
library(rsvd)

library(univarclust)
library(aricode)

# ---- SETTINGS ----------------------------------------------------------

dist_arg = "euclidean"
linkage_arg = "ward.D2"
center_arg = TRUE
scale_arg = FALSE
k_svd = 3

setwd("")
wd_var = paste0("k_svd_", k_svd)
dir.create(wd_var, showWarnings = FALSE, recursive = TRUE)

# ---- FONCTIONS ----------------------------------------------------------

directClustering <- function(dataSets) {
  hclust(dist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
}

averagedClustering <- function(dataSets) {
  hclust(Reduce("+", lapply(dataSets, dist, method = "euclidean")) / length(dataSets), method = "ward.D2")
}

mergeTreesWard <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
} 

# ------ Source : ----


# ------ Prepartion des donnees ------------------------------------------------

load("data_tablesformat.RData")

dim(dat1)
dim(dat2)

dat1 = dat1[which(apply(dat1, 1, var)!=0),]
dat2 = dat2[which(apply(dat2, 1, var)!=0),]
dim(dat1)
dim(dat2)

dataSets = list("Gene Expression" = t(dat1), "Methylation" = t(dat2))
dataSets = lapply(dataSets, scale, center = TRUE, scale = FALSE)
dataSets = lapply(dataSets, FUN = function(dat){dat = dat/svd(dat,  nu = 0, nv = 0)$d[1]})

lapply(dataSets, dim)

dat3 = matrix(NA, ncol = ncol(dataSets[[1]])+ncol(dataSets[[2]]), nrow = nrow(dataSets[[1]]))
dat3[,1:ncol(dataSets[[1]])] = dataSets[[1]]
dat3[,(ncol(dataSets[[1]])+1):(ncol(dat3))] = dataSets[[2]]

dim(dat3)

dist1 = dist(dataSets[[1]], method = dist_arg)
dist2 = dist(dataSets[[2]], method = dist_arg)

rSVD <- lapply(dataSets, rsvd, k = k_svd)
rSVD_dataSets = lapply(rSVD, FUN = function(svd_res) svd_res$u%*% diag(svd_res$d))
names(rSVD_dataSets) = paste0(names(dataSets), "_sp")

rSVD_dat3 = rsvd(dat3, k = k_svd)

# ---- Dossier des résultats : -----

setwd(wd_var)

# ------ Arbres sur tables séparées : -----------------------------------------------

hc_list = list("Gene expression" = hclust(dist1, method = linkage_arg),
               "Methylation" = hclust(dist2, method = linkage_arg))

hc_list_sp = list("Gene expression sp" = hclust(dist(rSVD_dataSets[[1]]), method = linkage_arg),
                  "Methylation sp" = hclust(dist(rSVD_dataSets[[2]]), method = linkage_arg))

# ------ Méthodes multivariees : -----------------------------------------------

hc_list_methods = list()

hc_list_methods$AD = hclust((dist1+dist2)/2, method = linkage_arg)
hc_list_methods$DC = hclust(dist(dat3, method = dist_arg), method = linkage_arg)
hc_list_methods$MT = mergeTrees(hc_list)

save(hc_list_methods, file = "hc_list_methods.RData")

# ------ Méthodes spectrales sur tables séparées : -----------------------------------------------

dist_sp_list = lapply(rSVD_dataSets, FUN = function(dat) dist(dat, method = dist_arg))
hc_sp_list = lapply(dist_sp_list, FUN = function(dist_mat) hclust(dist_mat, method = linkage_arg))

hc_list_methods$SDC = directClustering(rSVD_dataSets)
hc_list_methods$SAD = averagedClustering(rSVD_dataSets)
hc_list_methods$SMT = mergeTrees(hc_sp_list)

save(hc_list_methods, file = "hc_list_methods.RData")

# ------ Méthodes univariées : -----------------------------------------------

dataSets_univar = as.list(data.frame(dat3))

hc_list_methods$ADuni = averagedClustering(dataSets_univar)
hc_list_methods$MTuni = mergeTreesWard(dataSets_univar)

# ------ Méthodes spectrales univariées : --------------------------------------------------------

# --- Sur tables concaténées

dataSets_spectral <- as.data.frame(rSVD_dat3$u %*% diag(rSVD_dat3$d))

hc_list_methods$ScDC = hclust(dist(dataSets_spectral, method = "euclidean"), method = "ward.D2")
hc_list_methods$ScADuni = averagedClustering(as.list(dataSets_spectral))
hc_list_methods$ScMTuni = mergeTreesWard(as.list(dataSets_spectral))


save(hc_list_methods, file = "hc_list_methods.RData")
save.image(file='Session5.RData')
load("Session5.RData")

# ------------ EXPLORATION DES RESULTATS : --------------------------

cutree_index_max = 20

get_nid <- function(clustering, reference) {
  clusterings <- cutree(clustering, seq.int(1:cutree_index_max)) %>% as.data.frame() %>% as.list()
  nid <- map_dbl(clusterings, ~NID(., reference))  
  nid
}

# ------- Tables séparées ------ 

hc_list = c(hc_list, hc_list_sp)

save(hc_list, file = "hc_list_data.RData")

nids_celltype_data <- map_df(hc_list, get_nid, phen1$`cell type:ch1`) %>%
  add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "DataSet", value = "NID", -ngroup) %>%
  add_column(clinical = "CellType")

nids_df = as.data.frame(nids_celltype_data)
mat_res_best_nids = matrix(NA, ncol = 2, nrow = length(hc_list))
compteur_methode = 1
for(method in unique(nids_df$DataSet)){
  subset_df = nids_df[nids_df$DataSet==method,]
  mat_res_best_nids[compteur_methode, 1:2] = c(subset_df$ngroup[which.min(subset_df$NID)], subset_df$NID[which.min(subset_df$NID)])
  compteur_methode = compteur_methode + 1
}  
colnames(mat_res_best_nids) =  c("Nb Groupes", "NID")
rownames(mat_res_best_nids) = unique(nids_df$DataSet)
mat_res_data = mat_res_best_nids

# -------- Méthodes : -----

nids_celltype_methods <- map_df(hc_list_methods, get_nid, phen1$`cell type:ch1`) %>%
  add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "Method", value = "NID", -ngroup) %>%
  add_column(clinical = "CellType")


nids_df = as.data.frame(nids_celltype_methods)
mat_res_best_nids = matrix(NA, ncol = 2, nrow = length(unique(nids_df$Method)))
compteur_methode = 1
for(method in unique(nids_df$Method)){
  subset_df = nids_df[nids_df$Method==method,]
  mat_res_best_nids[compteur_methode, 1:2] = c(subset_df$ngroup[which.min(subset_df$NID)], subset_df$NID[which.min(subset_df$NID)])
  compteur_methode = compteur_methode + 1
}  
colnames(mat_res_best_nids) =  c("Nb Groupes", "NID")
rownames(mat_res_best_nids) = unique(nids_df$Method)
mat_res_methods = mat_res_best_nids

save(mat_res_data, mat_res_methods, file = "mat_res.RData")

write.table(mat_res_data, file = "mat_res_data.txt", quote = FALSE, sep = "\t")
write.table(mat_res_methods, file = "mat_res_methods.txt", quote = FALSE, sep = "\t")
