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
  // IntegerVector mergeMatrix_col1;
  // IntegerVector mergeMatrix_col2;

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

    // Reinitialized every step.
    IntegerVector vectorToMatch1(n) ;
    IntegerVector vectorToMatch2(n) ;
    IntegerVector which_child;


    current_element = prune_res(r_element,i);  //     elementCourant = matrice_aide[l_element,1]
    // std::cout << current_element << " - ";
    parent_group = prune_res(r_parent,i);  //     groupe_parent = matrice_aide[l_parent, 1]
    // std::cout << parent_group << '\n';

    for(int j = 0; j<n; j++){vectorToMatch1(j) = prune_res(r_group,j) ;}
    which_merge = match_func_int(parent_group, vectorToMatch1) ; //     which_merge = which(matrice_aide[l_groupeAct,]==groupe_parent)

    // //   merge_element = matrice_aide[,which_merge] // retire version cpp
    // old_group = prune_res(r_group, i) ;     //     groupe_act_old = matrice_aide[l_groupeAct,1]
    // // std::cout << old_group << " - ";
    // new_group = prune_res(r_group, which_merge) ;    //     groupe_act_new = matrice_aide[l_groupeAct, which_merge]
    // // std::cout << new_group << " - ";
    //
    // for(int j = 0; j<n; j++){vectorToMatch2(j) = prune_res(r_child,j) ;}
    // which_child = match_func_vectorInt(old_group, vectorToMatch2);  // ok jusque la

    // std::cout << i <<  " - " ; // Colonne consideree
    // std::cout << which_merge << " - "; // Colonne a merger
    // std::cout << prune_res(r_element, which_merge) << " - " ; // Element a merger
    // std::cout << which_child << " - "; // Elements issus du meme groupe
    // std::cout << which_child.length() << '\n' ;

    // for(int j=0; j<which_child.length(); j++){ // pas encore teste
    //     prune_res(r_child, which_child(j)) = new_group ;
    // } //     matrice_aide[l_enfant,which(matrice_aide[l_enfant,]==groupe_act_old)] = groupe_act_new

    mergeMatrix(i, 0) = current_element ;   //     merge_mat[iterations, 1] = elementCourant
    mergeMatrix(i, 1) = prune_res(r_element, which_merge) ;   //     merge_mat[iterations, 2] = merge_element[l_element]
    //
    prune_res(r_element, which_merge) = next_element ;     //     matrice_aide[l_element, which_merge] = elementSuivant
    // //     matrice_aide[l_groupeAct, which_merge] = merge_element[l_groupeAct] // inutile
    // //     matrice_aide[l_parent, which_merge] = merge_element[l_parent] // inutile
    // //     matrice_aide[l_enfant, which_merge] = merge_element[l_enfant] //inutile
    //
    // //     matrice_aide = matrice_aide[,-1] // Ne peut pas se faire comme ca en cpp
    prune_res(r_element,i) = -2*n;
    prune_res(r_group, i) = -2*n;
    prune_res(r_parent, i) = -2*n;
    prune_res(r_child, i) = -2*n;
    //
    next_element = next_element+1 ; //     elementSuivant = elementSuivant+1
    i = i+1;     //     iterations = iterations+1
  }
  //
  // // Fusion of the two last columns.
  mergeMatrix(i,0) = prune_res(r_element,(n-2)) ;
  mergeMatrix(i,1) = prune_res(r_element,(n-1)) ;

  return(mergeMatrix);
}

/*** R
# load("/home/hulot/Documents/packages_R/Rmergetrees/prune_res.RData")
# mergeMatrix = createMergeMatrix(ncol(matrice_aide2), matrice_aide2)
# mergeMatrix[1:10,]
# test = match_func_vectorInt(17, matrice_aide2[3,]) # Nombre d'enfants engendres par 17
# test
# test2 = match_func_vectorInt(17, matrice_aide2[2,]) # Logiquement un seul : nombre d'element dans le groupe 17
# test2
*/

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
