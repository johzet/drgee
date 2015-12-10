#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP tapplysum(SEXP Uin, SEXP ID) {
  NumericMatrix U_tmp(Uin);
  mat U = as<mat>(U_tmp);

  IntegerVector id_tmp(ID);
  uvec id = as<uvec>(id_tmp);

  int n_obs = U.n_rows;
  int n_col = U.n_cols;
  int n_clust = (unique(id)).n_elem;
  // int n_clust(nclust);

  mat Uout = as<mat>(n_clust, n_col);

  // Initialize the cluster identifier flag
  int clust_id = id[0];
  // Initialize the cluster number flag
  int clust_idx = 0;


  for(int k = 0; k < n_obs; ++k){
    if( id[k] != clust_id){
      clust_id = id[k];
      clust_idx += 1;
    }

    Uout.row(clust_idx) = Uout.row(clust_idx) + U.row(k)

  }

  return Uout;

}
