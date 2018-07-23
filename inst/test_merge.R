# rm(list = ls())
library(devtools)
devtools::install_github("AudreH/Rmergetrees")
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

tree.list.stdSum = lapply(tree.list, FUN = function(x){
  x$height = x$height/sum(x$height) # tous compris entre 0 et 1
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
# source("https://bioconductor.org/biocLite.R")
# biocLite("DECIPHER")
library(DECIPHER)
ape_consensus_p1 = consensus(as.phylo(hc_1), as.phylo(hc_2), p = 1, check.labels = TRUE) # strict consensus
ape_consensus_p05 = consensus(as.phylo(hc_1), as.phylo(hc_2), p = 0.5, check.labels = TRUE) # majority rule (p = 50%) consensus
tf <- tempfile()
write.tree(ape_consensus_p05, file = "ape_consensus_p05.tree", append = FALSE, digits = 10, tree.names = FALSE)
dend_consensus_p05 <- ReadDendrogram("ape_consensus_p05.tree", internalLabels = FALSE)

# Manually add branches lengths to consensus trees (to transform them into hclust objects)
ape_consensus_p05$edge.length = as.phylo(hc_1)$edge.length
ape_consensus_p05 <- force.ultrametric(ape_consensus_p05, method=c("nnls"))
ape_consensus_p1 <- force.ultrametric(ape_consensus_p1, method=c("nnls"))
# ape_consensus_p1 = force.binary()
# ape_consensus_p1$edge.length = seq(from = 1, to = nrow(ape_consensus_p1$edge), by = 1)
# is.ultrametric(ape_consensus_p05)
#
par(mfrow = c(2,2))
plot(hc_1, main= "Ward")
plot(hc_2, main= "complete")
plot(as.hclust(ape_consensus_p1), main= "Strict consensus - ape")
plot(as.dendrogram(ape_consensus_p05), main= "Majority rule consensus - ape")
# affreux ?

library(aricode)
par(mfrow = c(1,1))
plot(apply(cutree(hc_1, k = 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(hc_2, k = 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(ape_consensus_p1, k = 1:10), 2, NID, iris$Species), col = "red")
lines(apply(cutree(ape_consensus_p05, k = 1:10), 2, NID, iris$Species), col = "magenta")



ape_consensus_p1 <- consensus(as.phylo(hc_1), as.phylo(hc_2), p = 1, check.labels = TRUE) # strict consensus
ape_consensus_p1$edge.length <- rep(1, nrow(ape_consensus_p1$edge))
cut(as.dendrogram(ape_consensus_p1), 2)

lapply(1:6, function(h) phytools::treeSlice(as.phylo(as.dendrogram.phylo(ape_consensus_p1)), slice = h))


plot(phytools::force.ultrametric(ape_consensus_p1))


ape_consensus_p1 <- consensus(as.phylo(hc_1), as.phylo(hc_2), p = 1, check.labels = TRUE) # strict consensus
D <- as.dendrogram.phylo(ape_consensus_p1)

cut(D, h=1)$upper
