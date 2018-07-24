#' Merge a set of hclust objet into a single tree
#'
#' @param hc.list a list with a least one hclust object to be merge in a single consensus tree
#' @param standardize a boolean indicating wether the heights of the different trees should be normalized before merged
#' @export
mergeTrees = function(hc.list, standardize = FALSE){

  n <- unique(sapply(hc.list, function(hc) length(hc$order)))
  stopifnot(is.integer(n))
  p <- length(hc.list)

  #############################################
  # ----- Standardization : -------------------
  #############################################
  # fix tree comparisons issues.
  if (standardize) {
    lapply(hc.list, function(hc){
      hc$height <<- hc$height/max(hc$height) # tous compris entre 0 et 1
    })
  }

  #############################################
  # ----- Labels : ---------------------------
  #############################################

  # In case the trees have different labels, no merging possible
  labels <- Reduce(intersect, lapply(hc.list, function(hc) hc$labels))
  if (!is.null(labels)) {
    stopifnot(length(labels) == n)
    hc.list <- lapply(hc.list, reorder_hc)
  }

  #############################################
  # ----- Reconstitution paths : -------------
  #############################################
  rules   <- lapply(hc.list, as.fusionTree)
  heights <- unlist(lapply(rules, function(l) l$height))
  rules_order <- cbind(rep(1:(n - 1), p), rep(1:p, each = n - 1))
  rules_order <- rules_order[order(heights, decreasing = TRUE), ]

  aggreg <- pruneSplits(listSetRules = rules, orderRules = rules_order, n)
  lambdaRules <- sapply(aggreg$index_rule[-1], function(i) {
    rules[[rules_order[i,2]]]$height[rules_order[i,1]]
  })

  #############################################
  # ----- Reconstitution hclust : ------------
  #############################################

  # merging matrix
  merging <- getMergeMatrix(
    aggreg$group,
    aggreg$parent,
    order(aggreg$group, decreasing = TRUE) - 1)
  mergeMatrix <- merging$merge
  mergeMatrix[mergeMatrix <= n] <- -mergeMatrix[mergeMatrix <= n]
  mergeMatrix[mergeMatrix >  n] <-  mergeMatrix[mergeMatrix >  n] - n

  # Topological order of the dendrogram
  Order <- export_order(mergeMatrix, merging$node_size)

  ## Output
  consensusTree <- structure(
    list(merge = mergeMatrix,
        height = rev(lambdaRules),
        order  = Order,
        labels = hc.list[[1]]$labels,
        method = "consensus",
        dist.method = NA),
    class = "hclust")

  consensusTree
}

