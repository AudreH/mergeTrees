#' Merge a set of hclust objet into a single tree
#'
#' The consensus tree is built from top to bottom. Trees from hc.list argument are transformed into a list of possible splits. Each split is characterized by its heights in the tree and
#' the two clusters it is creating. If the trees are non-binary trees (i.e. create more than two clusters in one split, or equivalently, if they merge more than two clusters), the method creates as many splits as
#' there are clusters created.
#' All the splits from the trees are then ordered by decreasing height. The method follow the list to find the active splits, and create a consensus tree.
#' The method can be summarized with this property: "At height h, if elements i and j are groupes in the same clusters in all the trees, then they are in the same cluster in the consensus tree", or, equivalently:
#' "At height h, if elements i and j are not in the same clusters in at least one of the trees, then they are not in the same cluster in the consensus tree".
#' @name mergeTrees
#' @param hc.list a list with at least one hclust object to be merge in a single consensus tree. No other tree format is supported. The trees should have the same leaves labels, otherwise the function will stop.
#' @param standardize a boolean indicating wether the heights of the different trees should be normalized before merged. Normalization is done by divinding the heights of fusion by their maximum in each tree.
#' @return A list of class hclust, being the consensus tree, with the following components: height, merge, method, order, and, if any, labels. For more information about these components, please see hclust function help page.
#' @author Audrey Hulot, \email{audrey.hulot@@inra.fr}, Julien Chiquet, Guillem Rigaill
#' @examples
#'   library(mergeTrees)
#'   M1 = matrix(c(0,2,3,4,2,0,1,5,3,1,0,7,4,5,7,0), ncol = 4, nrow = 4)
#'   M2 = matrix(c(0,1,5,6,1,0,7,9,5,7,0,2,6,9,2,0), ncol = 4, nrow = 4)
#'   h1 = hclust(as.dist(M1))
#'   h2 = hclust(as.dist(M2))
#'   MT = mergeTrees(list(h1, h2))
#'   oldpar <- par(mfrow = c(1,3))
#'   plot(h1)
#'   plot(h2)
#'   plot(MT)
#'   par(oldpar)
#' @export
mergeTrees = function(hc.list, standardize = FALSE){

  n <- unique(sapply(hc.list, function(hc) length(hc$order)))
  stopifnot(is.integer(n))
  p <- length(hc.list)

  #############################################
  # ----- Standardization : -------------------
  #############################################
  # fix tree comparisons issues.
  if (standardize) { # each tree is standardized so the heights are all between 0 and 1
    hc.list <- lapply(hc.list, function(hc){
      # hc$height <<- hc$height/max(hc$height)  # 23/04/19 Probleme avec cette ligne
      hc$height <- hc$height/max(hc$height)
      hc
    })
  }

  #############################################
  # ----- Labels : ---------------------------
  #############################################

  # In case the trees have different labels, no merging possible
  labels <- Reduce(intersect, lapply(hc.list, function(hc) hc$labels))
  if (!is.null(labels)) {
    stopifnot(length(labels) == n)
    # hc.list <- lapply(hc.list, reorder_hc)
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

