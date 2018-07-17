# rm(list = ls())
library(devtools)
# devtools::install_github("AudreH/Rmergetrees")
build()
install()
library(Rmergetrees)

hc_1 <- hclust(dist(iris[, 1:4], "euclidean"), method = "ward.D2")
hc_2 <- hclust(dist(iris[, 1:4], "euclidean"), method = "complete")

tree.list = list(hc_1, hc_2)

hc_merged <- mergeTrees(list(hc_1, hc_2))

par(mfrow = c(1,3))
plot(hc_1, main= "Ward")
plot(hc_2, main= "complete")
plot(hc_merged, main= "merged")

library(aricode)
par(mfrow = c(1,1))
plot(apply(cutree(hc_1, 1:10), 2, NID, iris$Species), type = "l", ylim = c(0,1))
lines(apply(cutree(hc_2, 1:10), 2, NID, iris$Species), col = "blue")
lines(apply(cutree(hc_merged, 1:10), 2, NID, iris$Species), col = "red")
