#include <Rcpp.h>
#include "match_func_int.h"
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerMatrix createMergeMatrix(int n,
                                IntegerMatrix prune_res){

  // Initialization: rows of prune_res
  int r_element = 0;
  int r_group = 1;
  int r_parent = 2;
  int r_child = 3;

  // Initialization out
  IntegerVector mergeMatrix_col1;
  IntegerVector mergeMatrix_col2;

  IntegerMatrix mergeMatrix(n-1,2);

  // Initialization for loop
  int i = 0 ;
  int current_element; // Element to consider
  int next_element = n+1; // Next element to create (new cluster)
  int parent_group; // Parent group of i
  int which_merge; // column index of element to merge with current element
  int old_group; // old current group of the element
  int new_group; // new current group of the element
  IntegerVector merge_element; // information on which_merge


  while(i<(n-2)){ // merge matrix has n-1 rows

    current_element = prune_res(r_element,i);
    parent_group = prune_res(r_parent,i);

    which_merge = match_func_int(parent_group, prune_res(r_group,_)) ;
    merge_element = prune_res(_,which_merge) ;
    old_group = prune_res(r_group, 1) ;
    new_group = prune_res(r_group, which_merge) ;


    mergeMatrix(i, 1) = current_element ;
    mergeMatrix(i, 2) = merge_element(r_element);

    prune_res(r_element, which_merge) = next_element ;
    prune_res(r_group, which_merge) = merge_element(r_group) ;
    prune_res(r_parent, which_merge) = merge_element(r_parent) ;
    prune_res(r_child, which_merge) = merge_element(r_child) ;

    prune_res(r_element,i) = -2*n;
    prune_res(r_group, i) = -2*n;
    prune_res(r_parent, i) = -2*n;
    prune_res(r_child, i) = -2*n;

    next_element = next_element+1 ;
    i = i+1;
  }

  // Fusion of the two last columns.
  mergeMatrix(i,1) = prune_res(r_element,(n-2)) ;
  mergeMatrix(i,2) = prune_res(r_element,(n-1)) ;

  return(mergeMatrix);
}
// CreationMatriceMerge <- function(n, matrice_aide){
//   l_element = 1 ; l_groupeAct = 2; l_parent = 3 ; l_enfant = 4 ;
//
//   merge_mat = matrix(NA, ncol = 2, nrow = n-1)
//     iterations = 1
//   elementSuivant = n+1
//
//   while(iterations<(n-1)){
//     elementCourant = matrice_aide[l_element,1]
//     groupe_parent = matrice_aide[l_parent, 1]
//
//     which_merge = which(matrice_aide[l_groupeAct,]==groupe_parent)
//     merge_element = matrice_aide[,which_merge]
//
//     groupe_act_old = matrice_aide[l_groupeAct,1]
//     groupe_act_new = matrice_aide[l_groupeAct, which_merge]
//
//     matrice_aide[l_enfant,which(matrice_aide[l_enfant,]==groupe_act_old)] = groupe_act_new
//     merge_mat[iterations, 1] = elementCourant
//     merge_mat[iterations, 2] = merge_element[l_element]
//     matrice_aide[l_element, which_merge] = elementSuivant
//     matrice_aide[l_groupeAct, which_merge] = merge_element[l_groupeAct]
//     matrice_aide[l_parent, which_merge] = merge_element[l_parent]
//     matrice_aide[l_enfant, which_merge] = merge_element[l_enfant]
//     matrice_aide = matrice_aide[,-1]
//
//     elementSuivant = elementSuivant+1
//     iterations = iterations+1
//   }
//
// # fusion des deux dernieres colonnes (il est cense n'en rester que deux)
//   merge_mat[iterations,] = c(matrice_aide[l_element,1], matrice_aide[l_element,2])
//
//     merge_mat[merge_mat<=n] <- -merge_mat[merge_mat<=n] # les elements commencent a 1 en R pas en Cpp
//     merge_mat[merge_mat>n] <- merge_mat[merge_mat>n]-n
//
//     return(merge_mat)
// }
