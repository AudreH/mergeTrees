# rm(list = ls())
<<<<<<< HEAD
library(devtools)
devtools::install_github("AudreH/Rmergetrees")
# use_rcpp()
# document()
=======
# devtools::install_github("AudreH/Rmergetrees")
>>>>>>> 146e485936c8aa40fe9d4fc30f917dc095d67656
# build()
# install()
library(Rmergetrees)
library(aricode)

hc_1 <- hclust(dist(iris[, 1:4], "euclidean"), method = "ward.D2")
hc_2 <- hclust(dist(iris[, 1:4], "euclidean"), method = "complete")

tree.list = list(hc_1, hc_2)

hc_merged <- mergeTrees(list(hc_1, hc_2))

par(mfrow = c(1,3))
plot(hc_1, main= "Ward")
plot(hc_2, main= "complete")
plot(hc_merged, main= "merged")


par(mfrow = c(1,1))
plot(apply(cutree(hc_1, 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(hc_2, 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged, 1:10), 2, NID, iris$Species), col = "red")

# VERIFICATION
par(mfrow = c(1,2))
plot(hc_1)
plot(mergeTrees(list(hc_1)))
#
# ---- Standardization : max height : ----

tree.list.stdMax = lapply(tree.list, FUN = function(x){
  x$height = x$height/max(x$height) # tous compris entre 0 et 1
  return(x)
})

hc_merged.stdMax <- mergeTrees(tree.list.stdMax)

par(mfrow = c(1,3))
plot(tree.list.stdMax[[1]], main= "Ward - Std Max height")
plot(tree.list.stdMax[[2]], main= "complete - Std Max height")
plot(hc_merged.stdMax, main= "merged")

par(mfrow = c(1,1))
plot(apply(cutree(tree.list.stdMax[[1]], 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(tree.list.stdMax[[2]], 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged.stdMax, 1:10), 2, NID, iris$Species), col = "red")

# ---- Standardization : sum of heights : ----

# tree.list.stdSum = lapply(tree.list, FUN = function(x){
#   x$height = x$height/sum(x$height) # tous compris entre 0 et 1
#   return(x)
# })

tree.list.stdSum = lapply(tree.list, FUN = function(x){
  x$height = x$height/(x$height[length(x$height)] + x$height[length(x$height)-1])# tous compris entre 0 et 1
  return(x)
})

hc_merged.stdSum <- mergeTrees(tree.list.stdSum)

par(mfrow = c(1,3))
plot(tree.list.stdSum[[1]], main= "Ward - Std Max height")
plot(tree.list.stdSum[[2]], main= "complete - Std Max height")
plot(hc_merged.stdSum, main= "merged")

par(mfrow = c(1,1))
plot(apply(cutree(tree.list.stdSum[[1]], 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(tree.list.stdSum[[2]], 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged.stdSum, 1:10), 2, NID, iris$Species), col = "red")

# ---- Standardization : cumsum of heights : (bad idea) ----

tree.list.stdCumSum = lapply(tree.list, FUN = function(x){
  y = x$height[x$height!=0]
  # if(y[1]<1) y = y/y[1]
  y = y/cumsum(y)
  x$height[x$height!=0] = y
  return(x)
})

hc_merged.stdCumSum <- mergeTrees(tree.list.stdCumSum)

par(mfrow = c(1,3))
plot(tree.list.stdCumSum[[1]], main= "Ward - Std Max height")
plot(tree.list.stdCumSum[[2]], main= "complete - Std Max height")
plot(hc_merged.stdCumSum, main= "merged")

par(mfrow = c(1,1))
plot(apply(cutree(tree.list.stdCumSum[[1]], 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(tree.list.stdCumSum[[2]], 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged.stdCumSum, 1:10), 2, NID, iris$Species), col = "red")

# ---- Standardization : height(1)-height(2) : ----

tree.list.stdDiff = lapply(tree.list, FUN = function(x){
  x$height = x$height/((x$height[length(x$height)]- x$height[length(x$height)-1])*2)
  # avec ca on penalise plus les bons arbres !!
  return(x)
})

hc_merged.stdDiff <- mergeTrees(tree.list.stdDiff)

par(mfrow = c(1,3))
plot(tree.list.stdDiff[[1]], main= "Ward - Std Max height")
plot(tree.list.stdDiff[[2]], main= "complete - Std Max height")

plot(hc_merged.stdDiff, main= "merged")

par(mfrow = c(1,1))
plot(apply(cutree(tree.list.stdDiff[[1]], 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(tree.list.stdDiff[[2]], 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged.stdDiff, 1:10), 2, NID, iris$Species), col = "red")

# ---- Standardization : height(2)/height(A) : ----

tree.list.stdMult = lapply(tree.list, FUN = function(x){
  x$height = x$height*((x$height[length(x$height)-1]/x$height[length(x$height)]))
  # avec ca on penalise plus les bons arbres !!
  return(x)
})

hc_merged.stdMult <- mergeTrees(tree.list.stdMult)

par(mfrow = c(1,3))
plot(tree.list.stdMult[[1]], main= "Ward - Std Max height")
plot(tree.list.stdMult[[2]], main= "complete - Std Max height")

plot(hc_merged.stdMult, main= "merged")

par(mfrow = c(1,1))
plot(apply(cutree(tree.list.stdMult[[1]], 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(tree.list.stdMult[[2]], 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged.stdMult, 1:10), 2, NID, iris$Species), col = "red")


#---- Verification: compute only one tree : ----

tree.list = list(hc_1)
hc_merged_1 = mergeTrees(tree.list)
par(mfrow = c(1,2))
plot(hc_1, main= "Ward")
plot(hc_merged_1, main= "merged")

setequal(hc_merged_1$height, hc_1$height) # TRUE
hc_merged_1$merge[1:10,] # ok
hc_1$merge[1:10,] # ok

# ---- Essai autres methods de consensus tree : package ape ----

library(ape)
library(phytools)
# install.packages("phylogram")
library(phylogram)
# source("https://bioconductor.org/biocLite.R")
# biocLite("DECIPHER")
library(DECIPHER)
ape_consensus_p1 = consensus(as.phylo(hc_1), as.phylo(hc_2), p = 1, check.labels = TRUE) # strict consensus
ape_consensus_p05 = consensus(as.phylo(hc_1), as.phylo(hc_2), p = 0.5, check.labels = TRUE) # majority rule (p = 50%) consensus

dend_ape_consensus_p05 = as.dendrogram.phylo(ape_consensus_p05)
dend_ape_consensus_p05
plot(dend_ape_consensus_p05)

ape_consensus_p05_2 = as.phylo(dend_ape_consensus_p05)
tS_05 = treeSlice(ape_consensus_p05_2, slice = 3)
tS_05
class(tS_05)
tS_05[[1]]
tS_05[[2]]
tS_05 = treeSlice(ape_consensus_p05_2, slice = 4)
tS_05

tS_list = lapply(2:5, FUN = function(h) treeSlice(ape_consensus_p05_2, slice = h))
par(mfrow = c(1,1))
plot(ape_consensus_p05_2)

library(dendextend)
cut_dend_05 = cut(dend_ape_consensus_p05, h = 1)
get_leaves_attr(as.dendrogram(cut_dend_05$upper), "label")
get_leaves_attr(cut_dend_05$lower, "label")

cut_dend_05$lower
class(cut_dend_05$lower)

dend_ape_consensus_p1 = as.dendrogram.phylo(ape_consensus_p1)
dend_ape_consensus_p1
plot(dend_ape_consensus_p1)
ape_consensus_p1_2 = as.phylo(dend_ape_consensus_p1)

#### OUIIIIIIII
sapply(dendextend_cut_lower_fun(dend_ape_consensus_p05,6, labels), length)

# cutree(dend_ape_consensus_p05, k = 3)

# as.hclust(dend_ape_consensus_p05)

# as.hclust(as.dendrogram(ape_consensus_p05))

# plot(ape_consensus_p05)
# cutree(as.dendrogram(ape_consensus_p05), k = 3)

# tf <- tempfile()
# write.tree(ape_consensus_p05, file = "ape_consensus_p05.tree", append = FALSE, digits = 10, tree.names = FALSE)
# dend_consensus_p05 <- ReadDendrogram("ape_consensus_p05.tree", internalLabels = FALSE)

# Manually add branches lengths to consensus trees (to transform them into hclust objects)
# ape_consensus_p05$edge.length = as.phylo(hc_1)$edge.length
# ape_consensus_p05 <- force.ultrametric(ape_consensus_p05, method=c("nnls"))
# ape_consensus_p1 <- force.ultrametric(ape_consensus_p1, method=c("nnls"))
# ape_consensus_p1 = force.binary()
# ape_consensus_p1$edge.length = seq(from = 1, to = nrow(ape_consensus_p1$edge), by = 1)
# is.ultrametric(ape_consensus_p05)
#
par(mfrow = c(2,2))
plot(hc_1, main= "Ward")
plot(hc_2, main= "complete")
plot(dend_ape_consensus_p1, main= "Strict consensus - ape")
plot(dend_ape_consensus_p05, main= "Majority rule consensus - ape")
# affreux ?

library(aricode)
par(mfrow = c(1,1))
plot(apply(cutree(hc_1, k = 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(hc_2, k = 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(ape_consensus_p1, k = 1:10), 2, NID, iris$Species), col = "red")
lines(apply(cutree(ape_consensus_p05, k = 1:10), 2, NID, iris$Species), col = "magenta")



<<<<<<< HEAD
# ---- Test timing / limites fonction d'agregation : ----

# nb_ind  = 1000
# nb_var = 1000
#
# dataSet1 = matrix(rnorm(nb_var*nb_ind, 0, 1), ncol = nb_var, nrow = nb_ind) # individuals = rows
# dataSet2 = matrix(rnorm(nb_var*nb_ind, 0, 2), ncol = nb_var, nrow = nb_ind) # individuals = rows
#
# hc_1 = hclust(dist(dataSet1, method = "euclidean"), method = "ward.D2")
# hc_2 = hclust(dist(dataSet2, method = "euclidean"), method = "ward.D2")
#
# res = mergeTrees(list(hc_1, hc_2))
#
# par(mfrow = c(1,3))
# plot(hc_1)
# plot(hc_2)
# plot(res)
#
# library(profr)
# res_prof = profr({mergeTrees(list(hc_1, hc_2))},  interval = 0.02)


# nb_ind  = seq(1000,5000, by = 1000)
# nb_var = 1000
#
# for(i in 1:length(nb_ind)){
#   print(nb_ind[i])
#   dataSet1 = matrix(rnorm(nb_var*nb_ind[i], 0, 1), ncol = nb_var, nrow = nb_ind[i]) # individuals = rows
#   dataSet2 = matrix(rnorm(nb_var*nb_ind[i], 0, 2), ncol = nb_var, nrow = nb_ind[i]) # individuals = rows
#
#   hc_1 = hclust(dist(dataSet1, method = "euclidean"), method = "ward.D2")
#   hc_2 = hclust(dist(dataSet2, method = "euclidean"), method = "ward.D2")
#
#   res = mergeTrees(list(hc_1, hc_2))
# }

# ---- Labels issues ----

nb_genes = 10
nb_indiv = 3
nb_trees = 2

# par(mfrow = c(1,3))
list.trees = lapply(1:nb_trees, FUN = function(x){
  set.seed(x)
  data_mat = matrix(rnorm(nb_genes*nb_indiv), ncol = nb_indiv, nrow = nb_genes)
  colnames(data_mat) = paste0("labels_", 1:ncol(data_mat))
  rownames(data_mat) = paste0("G_", 1:nrow(data_mat))
  data_mat = data_mat[sample(1:nrow(data_mat), nrow(data_mat), replace = FALSE),
                      sample(1:ncol(data_mat), ncol(data_mat), replace = FALSE)]
  # plot(hclust(dist(data_mat)))
  return(hclust(dist(data_mat)))
})

h1 = list.trees[[1]]
h2 = list.trees[[2]]

h1$labels
h1$order

h2$labels
h2$order

new_order = match(h1$labels, h2$labels)
h2$labels = h2$labels[new_order]
# h2$order = h2$order[new_order]
h2$order = new_order


h2$order
h2$labels

h1$labels
h2$labels

h1 = list.trees[[1]]
h2 = list.trees[[2]]
mergedTree_noStd = mergeTrees(hc.list = list.trees, standardize = FALSE)
par(mfrow = c(1,3))
plot(h1, ylim = c(0, max(h1$height, h2$height)))
plot(h2, ylim = c(0, max(h1$height, h2$height)))
plot(mergedTree_noStd)
=======
ape_consensus_p1 <- consensus(as.phylo(hc_1), as.phylo(hc_2), p = 1, check.labels = TRUE) # strict consensus
ape_consensus_p1$edge.length <- rep(1, nrow(ape_consensus_p1$edge))
lapply(1:6, function(h) phytools::treeSlice(as.phylo(as.dendrogram.phylo(ape_consensus_p1)), slice = h))


plot(phytools::force.ultrametric(ape_consensus_p1))
>>>>>>> 146e485936c8aa40fe9d4fc30f917dc095d67656
