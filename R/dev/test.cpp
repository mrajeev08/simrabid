#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerVector tabulate4(const IntegerVector& x, const unsigned max) {
  IntegerVector counts(max);
  auto start = x.begin();
  auto end = x.end();
  for (auto it = start; it != end; it++) {
    if (*(it) > 0 && *(it) <= max)
      counts[*(it) - 1]++;
  }
  return counts;
}
