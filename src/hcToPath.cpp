#include <Rcpp.h>
// #include <vector>
// #include <iostream>
// #include <string>
using namespace Rcpp;


// // Function match for Cpp (character)
// // https://github.com/hadley/adv-r/blob/master/extras/cpp/match.cpp

// Fonction match faite dans l'optique ou de toute facon on trouvera les elements de "x" dans ceux d'"y"
// et une seule et unique fois.

// [[Rcpp::export]]
int match_func_int(int x, IntegerVector y){
  int out; // integer de la position de x dans y
  int position = 0 ; 
  int j = 0;
  while(position == 0){
    if(x == y(j)){position = 1; out = j ;}
    j++;
  }
  return out ;
}

// [[Rcpp::export]]
List hcToPath_cpp(List successives_steps,
                  IntegerMatrix merge1, 
                  IntegerVector order1, // needs to be all negative !!!! convert before apply this function
                  int n
                 ){
  
  // printf("characters\n") ;
  
  // merge is ordered according to height already.

  // Initialization: vectors to build the response matrix
  IntegerVector minVector;
  IntegerVector splitVector;
  IntegerVector maxVector;
  // IntegerVector heightVector(n, 0); // height already known...

  int i = 0 ;

  while(i<(n-1)){ // merge a n-1 lignes, donc de 0 a (n-2) ?
  
    if((merge1(i,0) < 0) & (merge1(i,1) < 0)){
      // Case 1: two single elements to merge
      // printf("CASE 1") ;
      
      int el1 = match_func_int(merge1(i,0), order1) ;
      int el2 = match_func_int(merge1(i,1), order1) ;
      int position = 0;
      if(el1>el2) position = 1 ;

      NumericVector vect_res = NumericVector::create(el1, el2);
      successives_steps.push_back(vect_res) ;

      int minus ;
      int split ;
      int maximum ;

      if(position == 1){//  max = el1, min = el2

        minus = el2 ;
        split = el2 ;
        maximum = el1 ;
        
      }else{

        minus = el1 ;
        split = el1 ;
        maximum = el2 ;

      }
      
      minVector.push_back(minus) ;
      splitVector.push_back(split) ;
      maxVector.push_back(maximum) ;

    } // End if case 1
    else if((merge1(i,0) < 0 ) & (merge1(i,1) > 0)){
      // Second case: first is a singleton, second is a cluster
      // printf("CASE 2") ;
      
      int el_singleton = match_func_int(merge1(i,0), order1) ; // element singleton
      NumericVector el_cluster = successives_steps[merge1(i,1)-1] ;

      int position = 0;
      if(el_singleton>max(el_cluster)) position = 1 ;

      int minus ;
      int split ;
      int maximum ;

      if(position == 1){ // singleton > elements in the cluster
        NumericVector vect_res = el_cluster;
        vect_res.insert(vect_res.end(), el_singleton) ;
        successives_steps.push_back(vect_res) ;

        minus = min(el_cluster) ;
        split = max(el_cluster) ;
        maximum = el_singleton ;

      }else{ // singleton < elements in the cluster
        NumericVector vect_res = el_cluster;
        vect_res.insert(vect_res.begin(), el_singleton) ;
        successives_steps.push_back(vect_res) ;

        minus = el_singleton ;
        split = el_singleton ;
        maximum = max(el_cluster) ;

      }

      minVector.push_back(minus) ;
      splitVector.push_back(split) ;
      maxVector.push_back(maximum) ;
    
    }// End else if case 2
    else if((merge1(i,0) > 0) & (merge1(i,1)< 0)){
      // Third case: first is a cluster, second is a singleton
      // printf("CASE 3") ;
      
      int el_singleton = match_func_int(merge1(i,1), order1) ; // element singleton
      NumericVector el_cluster = successives_steps[merge1(i,0)-1] ;

      int position = 0;
      if(el_singleton>max(el_cluster)) position = 1 ;

      int minus ;
      int split ;
      int maximum ;

      if(position == 1){ // singleton plus grand que les elements du cluster
        NumericVector vect_res = el_cluster;
        vect_res.insert(vect_res.end(), el_singleton) ;
        successives_steps.push_back(vect_res) ;

        minus = min(el_cluster) ;
        split = max(el_cluster) ;
        maximum = el_singleton ;

      }else{ // singleton plus petit que les elements du cluster
        NumericVector vect_res = el_cluster;
        vect_res.insert(vect_res.begin(), el_singleton) ;
        successives_steps.push_back(vect_res) ;

        minus = el_singleton ;
        split = el_singleton ;
        maximum = max(el_cluster) ;

      }

      minVector.push_back(minus) ;
      splitVector.push_back(split) ;
      maxVector.push_back(maximum) ;
    }// End else if case 3
    else if((merge1(i,0) > 0) & (merge1(i,1) > 0)){
      // Last case: both are already defined clusters
      // printf("CASE 4") ;
      
      NumericVector el_cluster1 = successives_steps[merge1(i,0)-1] ;
      NumericVector el_cluster2 = successives_steps[merge1(i,1)-1] ;
      
      int position = 0;
      if(max(el_cluster1)<max(el_cluster2)) position = 1 ;

      int minus ;
      int split ;
      int maximum ;

      if(position==1){ // Cluster 1 < cluster 2
        NumericVector vect_res(el_cluster1.length()+el_cluster2.size(),0) ;
        for(int j = 0; j<el_cluster1.length(); j++){
          int el = el_cluster1(j) ;
          vect_res.push_back(el) ;
        }
        for(int j = 0; j<el_cluster2.length(); j++){
          int el = el_cluster2(j) ;
          vect_res.push_back(el) ;
        }
        successives_steps.push_back(vect_res) ;
        minus = min(el_cluster1) ;
        split = max(el_cluster1) ;
        maximum = max(el_cluster2) ;

      }else{// Cluster 1 > cluster 2
        NumericVector vect_res(el_cluster1.length()+el_cluster2.size(),0) ;
        for(int j = 0; j<el_cluster2.length(); j++){
          int el = el_cluster2(j) ;
          vect_res.push_back(el) ;
        }
        for(int j = 0; j<el_cluster1.length(); j++){
          int el = el_cluster1(j) ;
          vect_res.push_back(el) ;
        }
        successives_steps.push_back(vect_res) ;
        minus = min(el_cluster2) ;
        split = max(el_cluster2) ;
        maximum = max(el_cluster1) ;
      }

      minVector.push_back(minus) ;
      splitVector.push_back(split) ;
      maxVector.push_back(maximum) ;
    } // End else if case 4

    i++;
  } // End while

  List out;
  out["minvector"] = minVector;
  out["splitVector"] = splitVector;
  out["maxVector"] = maxVector;
  return(out);
}

