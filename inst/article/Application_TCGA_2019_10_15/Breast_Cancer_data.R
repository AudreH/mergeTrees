rm(list = ls())

library(aricode)
library(devtools)
library(mergeTrees)
library(RColorBrewer)
library(dendextend)
library(tidyverse)
library(viridis)
library(rsvd)
library(svd)

setwd("/home/hulot/Documents/packages_R/mergeTrees/inst/article/Application_TCGA_2019_10_15/")

# Settings
dist_arg = "euclidean" 
linkage_arg = "ward.D2"

## Scale/center settings
center_arg = TRUE
scale_arg = FALSE

## SVD settings
k_svd = 5

# Fonctions

directClustering <- function(dataSets) {
  hclust(dist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
}

averagedClustering <- function(dataSets) {
  AD <- Reduce("+", lapply(dataSets, dist, method = "euclidean")) / length(dataSets)
  hclust(AD, method = "ward.D2")
}

mergeTreesWard <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
}

# ---- Chargement des données ---- 
load("tcga_brca_data.RData")
load("clinical.RData")

clinic1 = clinic1[order(clinic1$bcr_patient_barcode),]

dataSets = list(
  "methyl" = zmethyl,
  "mirna" = zmirna,
  "protein" = zprotein,
  "rna" = log2(zrna+1))

dataSets = lapply(dataSets, FUN = function(dat) dat[order(rownames(dat)),])

dataSets_0 = dataSets = lapply(dataSets, scale, center = center_arg, scale = scale_arg)

dataSets = lapply(dataSets, FUN = function(dat){
  dat = dat/svd(dat,  nu = 0, nv = 0)$d[1] # division par premiere valeur singuliere
})

# Traitement des données cliniques
clinical = clinic1
rownames(clinical) = clinical$bcr_patient_barcode
clinical = clinical[,-which(colnames(clinical)=="bcr_patient_barcode")]


# ---- Application des méthodes -----
hc_list_methods = list() 

## Méthodes multivariées
hc_list = lapply(dataSets, FUN = function(x) hclust(dist(x, method = dist_arg), method = linkage_arg))

hc_list_methods$AD = averagedClustering(dataSets)
hc_list_methods$DC = directClustering(lapply(dataSets, scale, center = center_arg, scale = scale_arg))
hc_list_methods$MT = mergeTrees(hc_list)

## Méthodes spectrales multivariées

### Spectral sur tables séparées

rSVD <- lapply(dataSets, rsvd, k = k_svd)
rSVD_dataSets = lapply(rSVD, FUN = function(svd_res) svd_res$u%*% diag(svd_res$d))

dist_sp_list = lapply(rSVD_dataSets, FUN = function(dat) dist(dat, method = dist_arg))
hc_sp_list = lapply(dist_sp_list, FUN = function(dist_mat) hclust(dist_mat, method = linkage_arg))

hc_list_methods$SDC = directClustering(rSVD_dataSets)
hc_list_methods$SAD = averagedClustering(rSVD_dataSets)
hc_list_methods$SMT = mergeTrees(hc_sp_list)

## Méthodes univariees

dataSets_univar = as.list(data.frame(Reduce("cbind", dataSets)))

hc_list_methods$ADuni = averagedClustering(dataSets_univar)
hc_list_methods$MTuni = mergeTreesWard(dataSets_univar)

## Méthodes spectrale univariees

### Spectral sur tables séparées

dataSets_sp <- lapply(lapply(dataSets, rsvd, k = k_svd), FUN = function(svd_res) svd_res$u%*% diag(svd_res$d))
data_sp_univar <- as.list(as.data.frame(do.call("cbind", dataSets_sp)))

hc_list_methods$SADuni = averagedClustering(data_sp_univar)
hc_list_methods$SMTuni = mergeTreesWard(data_sp_univar)

### Spectral sur tables concaténées 

rSVD <- rsvd(do.call("cbind", dataSets), k = k_svd)
dataSets_spectral <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))

hc_list_methods$ScDC = hclust(dist(as.data.frame(rSVD$u %*% diag(rSVD$d)), method = "euclidean"), method = "ward.D2")
hc_list_methods$ScADuni = averagedClustering(dataSets_spectral)
hc_list_methods$ScMTuni = mergeTreesWard(dataSets_spectral)

## Sectral sur tables de base

rSVD <- lapply(dataSets, rsvd, k = k_svd)
rSVD_dataSets = lapply(rSVD, FUN = function(svd_res) svd_res$u%*% diag(svd_res$d))
names(rSVD_dataSets) = paste0(names(rSVD_dataSets), "sp")

hc_list = c(hc_list, lapply(rSVD_dataSets, FUN = function(x) hclust(dist(x, method = dist_arg), method = linkage_arg)))

# ---- Comparaison des arbres des tables de base ----

NID_compare = function(tree_1, tree_2, cut_index_max = NULL){
  if(is.null(cut_index_max)) cut_index_max = length(tree_1$order)
  unlist(lapply(2:cut_index_max, FUN = function(cut_index) NID(cutree(tree_1, k = cut_index), 
                                                               cutree(tree_2, k = cut_index))))
}

mat_NID_compare = sapply(hc_list, function(x) sapply(hc_list, function(y) min(NID_compare(x,y, cut_index_max = 20))))

#  NID exploration fonction @ Julien

## Parametre pour niveau de coupure maximum

cutree_index_max = 104

## Fonction get_nid
get_nid <- function(clustering, reference) {
  clusterings <- cutree(clustering, seq.int(1:cutree_index_max)) %>% as.data.frame() %>% as.list()
  nid <- map_dbl(clusterings, ~NID(., reference))  
  nid
}

## Fonction figure 
plot_method = function(arbres_liste){
  nids_ER_status <- map_df(arbres_liste, get_nid, clinical$ER_status) %>%
    add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "DataSet", value = "NID", -ngroup) %>%
    add_column(clinical = "ER Status")
  nids_PR_status <- map_df(arbres_liste, get_nid, clinical$PR_status) %>%
    add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "DataSet", value = "NID", -ngroup) %>%
    add_column(clinical = "PR Status")
  nids_subtype <- map_df(arbres_liste, get_nid, clinical$subtype) %>%
    add_column(ngroup = seq.int(1:cutree_index_max)) %>% gather(key = "DataSet", value = "NID", -ngroup) %>%
    add_column(clinical = "Subtype")
  nids <- rbind(nids_ER_status, nids_PR_status, nids_subtype)
  
  nids %>% group_by(DataSet) %>% 
    ggplot(aes(x = ngroup, y = NID, color = DataSet))  + geom_line() + facet_grid(.~clinical) + theme_bw() + theme(legend.position="bottom")+
    scale_color_viridis(discrete = TRUE) -> plot_data
  
  return(list(nids = nids, plot_data = plot_data))
}

res = plot_method(hc_list)
print(res$plot_data)

### Meilleur NID possible et nombre de groupes associé


nids_df = as.data.frame(res$nids)
mat_res_best_nids = matrix(NA, ncol = 6, nrow = length(hc_list))
compteur_clinique = 1
for(clinique in unique(nids_df$clinical)){
  compteur_methode = 1
  for(dataset in unique(nids_df$DataSet)){
    subset_df = nids_df[nids_df$DataSet==dataset & nids_df$clinical==clinique,]
    mat_res_best_nids[compteur_methode, compteur_clinique:(compteur_clinique+1)] = c(subset_df$ngroup[which.min(subset_df$NID)], subset_df$NID[which.min(subset_df$NID)])
    compteur_methode = compteur_methode + 1
  }  
  compteur_clinique = compteur_clinique+2
}
rownames(mat_res_best_nids) = names(hc_list)
mat_res_sp_data = mat_res_best_nids

round(mat_res_sp_data, 2)

## Méthodes d'agrégation 

### Figures (code Julien)

#### Toutes les méthodes ensemble

res = plot_method(hc_list_methods)
print(res$plot_data)

nids_df = data.frame(res$nids)
nids_df$Method = nids_df$DataSet
mat_res_best_nids = matrix(NA, ncol = 6, nrow = length(hc_list_methods))
compteur_clinique = 1
for(clinique in unique(nids_df$clinical)){
  compteur_methode = 1
  for(method in unique(nids_df$Method)){
    subset_df = nids_df[nids_df$Method==method & nids_df$clinical==clinique,]
    mat_res_best_nids[compteur_methode, compteur_clinique:(compteur_clinique+1)] = c(subset_df$ngroup[which.min(subset_df$NID)], subset_df$NID[which.min(subset_df$NID)])
    compteur_methode = compteur_methode + 1
  }  
  compteur_clinique = compteur_clinique+2
}
rownames(mat_res_best_nids) = unique(nids_df$Method)
mat_res_methods = mat_res_best_nids

round(mat_res_methods, 2)

