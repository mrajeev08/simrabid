#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sim_trans(NumericVector row_id, NumericVector S,
               NumericVector E, NumericVector I, NumericVector V,
               int nlocs) {

  IntegerVector counts = rep(0, nlocs);
  for (int i = 0; i < nlocs; i++) {
    if (row_id[i] > 0 && row_id[i] <= nlocs)
      counts[row_id[i] - 1]++;
  }

  return counts;
  IntegerVector exps = tabulatecpp(row_id, nlocs);
  arma::uvec id_list = arma::find(exps > 0); // Find indices
  int nexps = row_id.length();

  // track who exposure was allocated to
  arma::uvec contact = rep(0, nexps);
  arma::uvec infected = rep(0, nexps);

  int n = id_list.length();
  for(int i = 0; i < n; ++i) {
    // Sample the states
    int ind = id_list[i];
    // taking out the infectious dog doing the biting (I[i] - 1)
    arma::uvec states = rep(c(1, 2, 3, 4),
                  c(S[ind], E[ind], ifelse(I[ind] > 0, I[ind] - 1, I[ind]),
                    V[ind]));
    int sample_length = states.length();

    if(sample_length > 0) {
      arma::uvec rows = arma::find(row_id >= ind);
      sample_length = ifelse(exps[ind] > sample_length,
                              sample_length, exps[ind]);
      arma::uvec rows_now = rows[1:sample_length];
      arma::uvec sampled_states = sample(states, size = sample_length,
                                 replace = FALSE);
      contact = arma_sub(contact, rows_now, sampled_states);
      // this updates infected to true if contact == 1
      infected = arma_sub_cond(infected, contact, 1, true);
    }
  }
  List L = List::create(v1, v2);
  return L;
}


// https://stackoverflow.com/questions/31001392/rcpp-version-of-tabulate-is-slower-where-is-this-from-how-to-understand
// [[Rcpp::export]]
IntegerVector tabulatecpp(const IntegerVector& x, const unsigned max) {
  IntegerVector counts(max);
  std::size_t n = x.size();
  for (std::size_t i = 0; i < n; i++) {
    if (x[i] > 0 && x[i] <= max)
      counts[x[i] - 1]++;
  }
  return counts;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
// https://gallery.rcpp.org/articles/armadillo-subsetting/index.html
// update a vector at index with replacement
arma::uvec arma_sub(arma::uvec x, arma::uvec pos, arma::uvec vals) {

    x.elem(pos) = vals;  // Subset by element position and set equal to

    return x;
}


// update based on condition in another vector
// return a boolean
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec arma_sub_cond(arma::rowvec x, arma::rowvec y,
                         double val, bool replace) {

  arma::uvec ids = arma::find(y == val); // Find indices

  x.elem(ids).fill(replace);       // Assign value to condition

  return x;
}

