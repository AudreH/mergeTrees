rm(list=ls())
library(Rmergetrees)

hc_1 <- hclust(dist(USArrests[1:5, 1:2]), method = "ward.D2")
hc_2 <- hclust(dist(USArrests[5:1, 3:4]), method = "ward.D2")
hc_3 <- hclust(dist(USArrests[sample(1:nrow(USArrests)), c(1,3)]), method = "ward.D2")

merged <- mergeTrees(list(hc_1, hc_2))
