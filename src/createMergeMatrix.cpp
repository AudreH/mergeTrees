#include <Rcpp.h>
#include<string.h>
#include "match_func_int.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix createMergeMatrix(int n,
                                IntegerMatrix prune_res){

  // Initialization: rows of prune_res
  int r_element = 0;
  int r_group = 1;
  int r_parent = 2;
  int r_child = 3;

  // Initialization out
  IntegerMatrix mergeMatrix(n-1,2);

  // Initialization for loop
  int i = 0 ;
  int current_element; // Element to consider
  int next_element = n+1; // Next element to create (new cluster)
  int parent_group; // Parent group of i
  int which_merge; // column index of element to merge with current element
  int old_group; // old current group of the element
  int new_group; // new current group of the element

  // IntegerVector merge_element; // information on which_merge

  // BEGIN LOOP
  while(i<(n-2)){ // merge matrix has n-1 rows

    // Reinitialized every step
    IntegerVector vectorToMatch1(n) ;
    IntegerVector vectorToMatch2(n) ;
    IntegerVector which_child;


    current_element = prune_res(r_element,i);
    parent_group = prune_res(r_parent,i);

    for(int j = 0; j<n; j++){vectorToMatch1(j) = prune_res(r_group,j) ;}
    which_merge = match_func_int(parent_group, vectorToMatch1) ;

    mergeMatrix(i, 0) = current_element ;
    mergeMatrix(i, 1) = prune_res(r_element, which_merge) ;
    prune_res(r_element, which_merge) = next_element ;

    // Replace values by -2*n so they are never find in match functions
    prune_res(r_element,i) = -2*n;
    prune_res(r_group, i) = -2*n;
    prune_res(r_parent, i) = -2*n;
    prune_res(r_child, i) = -2*n;

    next_element = next_element+1 ;
    i = i+1;
  } // END LOOP

  // Fusion of the two last columns.
  mergeMatrix(i,0) = prune_res(r_element,(n-2)) ;
  mergeMatrix(i,1) = prune_res(r_element,(n-1)) ;

  return(mergeMatrix);
}

/*** R
## Test for the function
# load("/home/hulot/Documents/packages_R/Rmergetrees/prune_res.RData")
# mergeMatrix = createMergeMatrix(ncol(matrice_aide2), matrice_aide2)
# mergeMatrix[1:10,]
# test = match_func_vectorInt(17, matrice_aide2[3,]) # Nombre d'enfants engendres par 17
# test
# test2 = match_func_vectorInt(17, matrice_aide2[2,]) # Logiquement un seul : nombre d'element dans le groupe 17
# test2
*/
