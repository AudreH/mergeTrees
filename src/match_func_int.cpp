#include <Rcpp.h>
using namespace Rcpp;

// Match function: assume there IS a match, and that only the first is interesting to return.
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

// Match function returning a vector: assume there IS at least one match.
// [[Rcpp::export]]
IntegerVector match_func_vectorInt(int x, IntegerVector y){
  IntegerVector out; // integer de la position de x dans y
  for(int i=0; i<y.length(); i++){
    if(x == y(i)){out.push_back(i) ;}
  }
  return out ;
}
