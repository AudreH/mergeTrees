library(Rmergetrees)

n <- 1000
p <- 100
X <- matrix(rnorm(n * p), n, p)

hc <- hclust(dist(X), method = "ward.D2")

hcToPath(hc)

as.fusionTree(hc)

library(microbenchmark)

res <- microbenchmark(hcToPath(hc), as.fusionTree(hc))
ggplot2::autoplot(res)
