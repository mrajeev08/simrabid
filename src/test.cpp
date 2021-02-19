#include <Rcpp.h>
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
double sumtrue(LogicalVector x) {
  return sum(x);
}

// [[Rcpp::export]]
NumericVector indbooldiff(NumericVector t, LogicalVector within) {
  NumericVector td = t[within];
  NumericVector tt = t[1] - td;
  return tt;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
sumtrue(c(0, 1, 0))
indbooldiff(c(1, 2, 5, 7), c(0, 1, 1, 0))
*/
