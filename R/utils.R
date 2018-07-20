reorder_hc_list <- function(hc_list) {

  # reference: les labels/order of the first tree in the list
  ref_label <- hc_list[[1]]$labels

  if (!is.null(ref_label)){
    hc_list <- lapply(hc_list, FUN = function(hc){
      if (!identical(hc$labels, ref_label)) {
        new_order <- match(ref_label, hc$labels)
        hc$labels <- hc$labels[new_order]
### NOT SURE YET OF THAT PART...
        ## hc$order <- hc$order[new_order]
        hc$order <- new_order
      }
      hc
    })
  }
  hc_list
}
