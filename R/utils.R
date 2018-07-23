reorder_hc <- function(hc) {
  tmp <- order(hc$labels)
  lexico_order <- (1:length(hc$labels))[order(tmp)]
  hc$labels <- hc$labels[tmp]
  hc$order <- lexico_order[hc$order]
  hc$merge[hc$merge < 0] <- -lexico_order[-hc$merge[hc$merge < 0]]
  hc
}

# reorder_hc_list <- function(hc_list) {
#
#   lapply()
#   # reference: les labels/order of the first tree in the list
#   ref_label <- hc_list[[1]]$labels
#
#   if (!is.null(ref_label)){
#     hc_list <- lapply(hc_list, FUN = function(hc){
#       if (!identical(hc$labels, ref_label)) {
#         new_order <- match(ref_label, hc$labels)
#         hc$labels <- hc$labels[new_order]
# ### NOT SURE YET OF THAT PART...
#         ## hc$order <- hc$order[new_order]
#         hc$order <- new_order
#       }
#       hc
#     })
#   }
#   hc_list
# }

##' @export
as.fusionTree <- function(hc_obj) {

  fusionTree <- as_fusionTree(hc_obj$merge, hc_obj$order)
  res <- list(path   = fusionTree[nrow(fusionTree):1,],
              height = rev(hc_obj$height),
              order  = hc_obj$order)
  res
}
