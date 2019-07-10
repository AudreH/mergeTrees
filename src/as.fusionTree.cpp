#include <Rcpp.h>
using namespace Rcpp;

IntegerVector get_fusion_desc(int child,
                               const IntegerVector& topoOrder,
                               const IntegerVector& idown,
                               const IntegerVector& iup) {

    IntegerVector fusion_desc (2) ;

    if (child < 0) {
      fusion_desc.fill(topoOrder[-child - 1]) ;
    } else {
      fusion_desc[0] = idown [child - 1] ;
      fusion_desc[1] = iup   [child - 1] ;
    }

    return fusion_desc ;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix as_fusionTree (
    IntegerMatrix merge,
    IntegerVector order) {

  int n = order.size() ;
  IntegerVector idown  (n-1) ;
  IntegerVector iup    (n-1) ;
  IntegerVector isplit (n-1) ;

  IntegerVector topoOrder (n) ;
  // topological index of each element
  for (int i = 0; i < n; i ++){
    topoOrder[order[i] - 1] =  i ;
  }

  for (int i = 0; i < (n-1); i ++) {

    // description of the fusion (indices of downer/upper elements in the fusion + splitting element)
    IntegerVector child1 =  get_fusion_desc (merge(i, 0), topoOrder, idown, iup) ;
    IntegerVector child2 =  get_fusion_desc (merge(i, 1), topoOrder, idown, iup) ;

    idown[i]  = std::min(child1[0], child2[0]) ;
    isplit[i] = std::min(child1[1], child2[1]) ;
    iup[i]    = std::max(child1[1], child2[1]) ;

  }

  IntegerMatrix fusion (n - 1, 3) ;
  fusion(_, 0) = idown  + 1;
  fusion(_, 1) = isplit + 1;
  fusion(_, 2) = iup    + 1;
  colnames(fusion) = CharacterVector::create("index_down", "index_split", "index_up");

  return fusion ;

}

//
// /*** R
// library(mergeTrees)
// hc_1 <- hclust(dist(USArrests[1:5, 1:2]), method = "ward.D2")
// as_fusionTree(hc_1$merge, hc_1$order)
// */
