// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

struct pos_node {
  int_fast32_t pos;
  int_fast32_t node;
};

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
Rcpp::List getMergeMatrix(IntegerVector group,
                              IntegerVector parent,
                              IntegerVector order) {

  // Initialization
  int n = order.size() ;
  IntegerMatrix merge(n - 1, 2);

  IntegerVector vanishing = order;
  IntegerVector tmp = parent[order] ;
  IntegerVector aggregating = order[(n-1) - tmp] ;
  IntegerVector label (n) ;
  for (int i = 0; i < n ; i++) {label[i] = i;}

  IntegerVector node_size(n - 1) ;
  IntegerVector group_size(n, 1) ;

  for (int i = 0; i < (n - 1) ; i++) {

    merge(i, 0) = label[vanishing[i]] ;
    merge(i, 1) = label[aggregating[i]] ;

    label[aggregating[i]] = i + n ;

    group_size[aggregating[i]] += group_size[vanishing[i]] ;
    node_size[i] = group_size[aggregating[i]] ;

  }

  return Rcpp::List::create(
    Rcpp::Named("merge") = merge + 1,
    Rcpp::Named("node_size") = node_size
  );
}

