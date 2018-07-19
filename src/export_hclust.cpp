// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

struct pos_node {
  int_fast32_t pos;
  int_fast32_t node;
};

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector export_order(const IntegerMatrix& merge, const IntegerVector& size) {
  // function to recover a correct order of the nodes for the dendrogram (required for hclust object)

  /* Parameters:
    N         : number of data points
  merge     : (N-1)×2 array which specifies the node indices which are merged in each step of the clustering procedure.
  Negative entries -1...-N point to singleton nodes, while positive entries 1...(N-1) point to nodes which are themselves parents of other nodes.
  size : array of node sizes - makes it easier
  order     : output array of size N
  Runtime: Θ(N)
  */

  int N = size.size() + 1 ;
  std::vector<pos_node> queue(N/2);

  int_fast32_t parent;
  int_fast32_t child;
  int_fast32_t pos = 0;
  IntegerVector order = IntegerVector (N);
  queue[0].pos = 0;
  queue[0].node = N-2;
  int_fast32_t idx = 1;

  do {
    --idx;
    pos = queue[idx].pos;
    parent = queue[idx].node;

    // First child
    child = merge(parent, 0);
    if (child<0) { // singleton node, write this into the 'order' array.
      order[pos] = -child;
      ++pos;
    }
    else { /* compound node: put it on top of the queue and decompose it in a later iteration. */
        queue[idx].pos = pos;
        queue[idx].node = child-1; // convert index-1 based to index-0 based
        ++idx;
        pos += size[child-1];
    }
    // Second child
    child = merge(parent,1);
    if (child<0) {
      order[pos] = -child;
    }
    else {
      queue[idx].pos = pos;
      queue[idx].node = child-1;
      ++idx;
    }
  } while (idx>0);

  return (order) ;
};

// [[Rcpp::export]]
Rcpp::IntegerMatrix export_merge(const IntegerVector& parent1, const IntegerVector& parent2) {

  // function to ouput the dendrogram of fusion into the hclust format
  int K = parent1.size() + 1 ;
  IntegerMatrix merge = IntegerMatrix(K - 1, 2);

  for (int k = 0; k < K - 1; k++) {

    int_fast32_t merge1, merge2;

    if (parent1[k] > K) {
      merge1 = parent1[k] - K;
    } else {
      merge1 = -parent1[k];
    }
    if (parent2[k] > K) {
      merge2 = parent2[k] - K;
    } else {
      merge2 = -parent2[k];
    }
    if (merge1 < merge2) {
      merge(k, 0) = merge1;
      merge(k, 1) = merge2;
    }
    else {
      merge(k, 0) = merge2;
      merge(k, 1) = merge1;
    }
  }
  return(merge) ;
};
