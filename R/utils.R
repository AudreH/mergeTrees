reorder_hc <- function(hc) {
  tmp <- order(hc$labels)
  lexico_order <- (1:length(hc$labels))[order(tmp)]
  hc$labels <- hc$labels[tmp]
  hc$order <- lexico_order[hc$order]
  hc$merge[hc$merge < 0] <- -lexico_order[-hc$merge[hc$merge < 0]]
  hc
}

##' Turn an hclust object into a fusion path.
##'
##' @param hc_obj and hclust object
##' @return a list with three element: path of fusion, vector of height, ordering.
as.fusionTree <- function(hc_obj) {

  fusionTree <- as_fusionTree(hc_obj$merge, hc_obj$order)
  res <- list(path   = fusionTree[nrow(fusionTree):1,],
              height = rev(hc_obj$height),
              order  = hc_obj$order)
  res
}
