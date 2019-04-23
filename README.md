Rmergetrees
================

A package for fastly merging tree-like objects.

## Install

``` r
devtools::install_github("AudreH/Rmergetrees")
```

```{r, echo = TRUE}
nb_genes = 50
nb_indiv = 3
nb_trees = 10

list.trees = lapply(1:nb_trees, FUN = function(x){
  set.seed(x)
  data_mat = matrix(rnorm(nb_genes*nb_indiv), ncol = nb_indiv, nrow = nb_genes)
  return(hclust(dist(data_mat)))
})
```

```{r, echo = TRUE}
mergedTree = merge.trees(hc.list = list.trees)
plot(mergedTree_noStd, main = "noStd")
```
