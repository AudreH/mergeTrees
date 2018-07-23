#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pruneSplits(List listSetRules, IntegerMatrix orderRules, int n) {

  int p = listSetRules.length() ;
  IntegerVector groups(n, 3);
  // INITIALIZATION //
  // for groups
  IntegerVector groupsNbIndividuals(n, 0);
  IntegerVector groupsNbCalled(n, 0);
  IntegerVector groupsParent(n, 0);
  IntegerVector groupsChildCurrent(n, 0);
  IntegerVector groupsIndexRules(n, 0);
  LogicalVector groupsToVisit(n, false);
  LogicalVector groupsToSplit(n, false);


  // set first group 0
  groupsNbIndividuals[0] = n;

  // for samples
  IntegerVector currentGroup(n, 0);


  // END INITIALIZATION //


  // CORE OF THE ALGO //
  int nbGroup=1;
  int  i_rule=0;
  int nbMoves=0;
  while( (nbGroup < n) & (i_rule < (n-1)*p) ){

    // 0) Initialize groups...
    int index = orderRules(i_rule, 0);
    int p = orderRules(i_rule, 1);

    List setRules = listSetRules[p-1];
    IntegerVector order_current = setRules["order"];
    IntegerMatrix rules_current = setRules["path"];

    // order in R -1 to get order in Cpp
    int j_low   = rules_current(index-1, 0);
    int j_split = rules_current(index-1, 1);
    int j_high  = rules_current(index-1, 2);

    // 1) DEFINE SMALLEST SUBSET [j_low, j_high]
    if( (j_split - j_low + 1) <= (j_high - j_split) ) {
      j_high = j_split;
    } else {
      j_low  = j_split + 1;
    }
    // END 1

    // 2) COUNT CALLS PER SUBSET
    for(int j=j_low; j <= j_high; j++){

      int groupOfJ = currentGroup[ order_current[j-1] -1 ];
      groupsNbCalled[groupOfJ] = groupsNbCalled[groupOfJ]+1;
      groupsToVisit[groupOfJ]  = true;
    }
    // END 2


    // 3) CHECK NB CALLS AND DEFINE NEW GROUPS (loop over subset)
    for(int j=j_low; j <= j_high; j++){

      int groupOfJ = currentGroup[ order_current[j-1] -1];

      if(groupsToVisit[groupOfJ] ){ // split (first time we check this group)
        groupsToVisit[groupOfJ] = false;
        if(groupsNbIndividuals[groupOfJ] != groupsNbCalled[groupOfJ]){ // do we need to split
          groupsToSplit[groupOfJ] = true;

          // set new child
          groupsParent[nbGroup]        = groupOfJ;
          groupsIndexRules[nbGroup]    = i_rule;
          groupsNbIndividuals[nbGroup] = groupsNbCalled[groupOfJ];

          // reset parent
          groupsChildCurrent[groupOfJ]  = nbGroup;
          groupsNbIndividuals[groupOfJ] = groupsNbIndividuals[groupOfJ] - groupsNbCalled[groupOfJ];
          groupsNbCalled[groupOfJ] = 0;

          // next group number
          nbGroup++;
        }
        else { // second we check this group no-need to do anything
          groupsToSplit[groupOfJ]  = false;

          // reset parent
          groupsNbCalled[groupOfJ] = 0;
        }
      }

    }
    // END 3

    // 4 SET TO NEW GROUPS (loop over subset)
    for(int j=j_low; j <= j_high; j++){

      int groupOfJ = currentGroup[ order_current[j-1] -1];
      if(groupsToSplit[groupOfJ] == true){
        currentGroup[ order_current[j-1] -1] = groupsChildCurrent[groupOfJ];
      }
    }
    // END 4


    i_rule++;
  }

  // END CORE //

  // finding group order
  IntegerVector sorted = clone(currentGroup).sort();

  return List::create(
    Rcpp::Named("parent")      = groupsParent[currentGroup],
    Rcpp::Named("group")       = currentGroup,
    Rcpp::Named("index_rule" ) = groupsIndexRules + 1
  );

}

