#include <Rcpp.h>
#include <string>
#include <vector>

# define square(x) ((x)*(x))
int match_func_int(int x, Rcpp::IntegerVector y) ;

Rcpp::IntegerVector match_func_vectorInt(int x, Rcpp::IntegerVector y);
