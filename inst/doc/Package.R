## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ---- include = FALSE----------------------------------------------------
library(knitr)

## ------------------------------------------------------------------------
library(Rmergetrees)

## ------------------------------------------------------------------------
data("iris")
head(iris)

## ------------------------------------------------------------------------
set.seed(2000)
iris = iris[sample(1:nrow(iris), 20, replace = FALSE),]
rownames(iris) = 1:nrow(iris)

## ------------------------------------------------------------------------
data_list = list("Sepal" = iris[,grep("Sepal", colnames(iris))], "Petal" = iris[,grep("Petal", colnames(iris))])

## ------------------------------------------------------------------------
library(dendextend) # labels_colors for highlighting groups in a dendrogram

## ---- fig.align = "center", fig.width = 7, fig.height = 4----------------
group_color = c("setosa" = "blue", "versicolor" = "black", "virginica" = "red")
par(mfrow = c(1,2), mar = c(2,2,2,0))
res = lapply(1:length(data_list), FUN = function(i){
  hc = hclust(dist(data_list[[i]], method = "euclidean"), method = "ward.D2")
  dend = as.dendrogram(hc)
  labels_colors(dend) = group_color[iris$Species][order.dendrogram(dend)]
  plot(dend, main = names(data_list)[i])
})

## ------------------------------------------------------------------------
DC = consensusTree(lis = data_list, method = "directClustering", data_st = TRUE, distance_method = "euclidean", linkage_method = "ward.D2", verbose = TRUE)
AD = consensusTree(lis = data_list, method = "averageDistance", dist_st = TRUE, distance_method = "euclidean", linkage_method = "ward.D2", verbose = TRUE)
MT = consensusTree(lis = data_list, method = "mergeTrees", dist_st = TRUE, distance_method = "euclidean", linkage_method = "ward.D2", verbose = TRUE)

## ---- fig.align = "center", fig.width = 9, fig.height = 4----------------
group_color = c("setosa" = "blue", "versicolor" = "black", "virginica" = "red")
list_res = list("DirectClustering" = DC, "AverageDistance" = AD, "MergeTrees" = MT)
par(mfrow = c(1,3), mar = c(2,2,2,0))
res = lapply(1:length(list_res), FUN = function(i){
  hc = list_res[[i]]
  dend = as.dendrogram(hc)
  labels_colors(dend) = group_color[iris$Species][order.dendrogram(dend)]
  plot(dend, main = names(list_res)[i])
})

## ---- out.width = "80%", out.height = "50%", fig.align = "center"--------
knitr::include_graphics("Tree_aggregation_description.png")

## ------------------------------------------------------------------------
devtools::install_github("jchiquet/fusedanova")

## ------------------------------------------------------------------------
library(fusedanova)

## ------------------------------------------------------------------------
sessionInfo()

