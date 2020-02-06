library(tidyverse)
library(devtools)
library(mergeTrees)

library(parallelDist)
library(parallel)

library(univarclust)
library(aricode)

# ---- SETTINGS -----------------------------------------------------------

dist_arg = "euclidean"
linkage_arg = "ward.D2"
center_arg = TRUE
scale_arg = FALSE

k_svd = 3

# ---- FONCTIONS CONSENSUS TREES ------------------------------------------

directClustering <- function(dataSets) {
  hclust(parDist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
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

# ------------FUNCTIONS FOR RESULTS : -------------------------------------

cutree_index_max = 20

get_nid <- function(clustering, reference) {
  clusterings <- cutree(clustering, seq.int(1:cutree_index_max)) %>% as.data.frame() %>% as.list()
  nid <- map_dbl(clusterings, ~NID(., reference))  
  nid
}

get_ari <- function(clustering, reference) {
  clusterings <- cutree(clustering, seq.int(1:cutree_index_max)) %>% as.data.frame() %>% as.list()
  ari <- map_dbl(clusterings, ~ARI(., reference))  
  ari
}

# ------ DATA -------------------------------------------------------------

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

dist1 = parDist(dataSets[[1]], method = dist_arg)
dist2 = parDist(dataSets[[2]], method = dist_arg)


SVD <- lapply(dataSets, svd, nv = 0) 
SVD_dat3 = svd(dat3, nv = 0) 

# save.image("Session_base.RData")
# load("Session_base.RData")

# ---- Data Reconstruction : -----

SVD_dataSets = lapply(SVD, FUN = function(svd_res) svd_res$u[,1:k_svd]%*% diag(svd_res$d[1:k_svd]))
names(SVD_dataSets) = paste0(names(dataSets), "_sp")

# ------ Arbres sur tables séparées : -------------------------------------

hc_list = list("Gene expression" = hclust(dist1, method = linkage_arg),
               "Methylation" = hclust(dist2, method = linkage_arg))

hc_list_sp = list("Gene expression sp" = hclust(dist(SVD_dataSets[[1]]), method = linkage_arg),
                  "Methylation sp" = hclust(dist(SVD_dataSets[[2]]), method = linkage_arg))

# ------ Non spectral : ---------------------------------------------------

hc_list_methods = list()

hc_list_methods$AD = hclust((dist1+dist2)/2, method = linkage_arg)
hc_list_methods$DC = hclust(dist(dat3, method = dist_arg), method = linkage_arg)
hc_list_methods$MT = mergeTrees(hc_list)

# --- Spectral : ----------------------------------------------------------

dataSets_spectral <- as.data.frame(SVD_dat3$u[,1:k_svd] %*% diag(SVD_dat3$d[1:k_svd]))

hc_list_methods$ScDC = hclust(dist(dataSets_spectral, method = "euclidean"), method = "ward.D2")
hc_list_methods$ScADuni = averagedClustering(as.list(dataSets_spectral))
hc_list_methods$ScMTuni = mergeTreesWard(as.list(dataSets_spectral))

# ------- Consensus Trees list : ------------------------------------------

hc_list = c(hc_list, hc_list_sp)

# ------- NIDS Individual tables  -----------------------------------------

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

# -------- NIDS Consensus trees : -----------------------------------------

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

# ------- ARIS Individual tables ------------------------------------------

aris_celltype_data <- map_df(hc_list, get_ari, phen1$`cell type:ch1`) %>%
  add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "DataSet", value = "ARI", -ngroup) %>%
  add_column(clinical = "CellType")

aris_df = as.data.frame(aris_celltype_data)
mat_res_best_aris = matrix(NA, ncol = 2, nrow = length(hc_list))
compteur_methode = 1
for(method in unique(aris_df$DataSet)){
  subset_df = aris_df[aris_df$DataSet==method,]
  mat_res_best_aris[compteur_methode, 1:2] = c(subset_df$ngroup[which.max(subset_df$ARI)], subset_df$ARI[which.max(subset_df$ARI)])
  compteur_methode = compteur_methode + 1
}  
colnames(mat_res_best_aris) =  c("Nb Groupes", "ARI")
rownames(mat_res_best_aris) = unique(aris_df$DataSet)
mat_res_data = mat_res_best_aris

# -------- ARIS Consensus trees : -----------------------------------------

aris_celltype_methods <- map_df(hc_list_methods, get_ari, phen1$`cell type:ch1`) %>%
  add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "Method", value = "ARI", -ngroup) %>%
  add_column(clinical = "CellType")


aris_df = as.data.frame(aris_celltype_methods)
mat_res_best_aris = matrix(NA, ncol = 2, nrow = length(unique(aris_df$Method)))
compteur_methode = 1
for(method in unique(aris_df$Method)){
  subset_df = aris_df[aris_df$Method==method,]
  mat_res_best_aris[compteur_methode, 1:2] = c(subset_df$ngroup[which.max(subset_df$ARI)], subset_df$ARI[which.max(subset_df$ARI)])
  compteur_methode = compteur_methode + 1
}  
colnames(mat_res_best_aris) =  c("Nb Groupes", "ARI")
rownames(mat_res_best_aris) = unique(aris_df$Method)
mat_res_methods = mat_res_best_aris

save(mat_res_data, mat_res_methods, file = "mat_res.RData")

write.table(mat_res_data, file = "mat_res_data_ari.txt", quote = FALSE, sep = "\t")
write.table(mat_res_methods, file = "mat_res_methods_ari.txt", quote = FALSE, sep = "\t")

