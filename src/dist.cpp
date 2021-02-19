//dist.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double distmean(NumericVector x, NumericVector y, NumericVector t,
                int twindow){

  int n = x.size();
  double dist = 0;
  int nobs = 0;

  // faster looping in CPP?
  for(int i = 0; i < n; ++i) {

    NumericVector tdiff = t[i] - t; // change this to loop
    LogicalVector within = tdiff > 0 & tdiff < twindow;

    NumericVector xd = x[within]; //calculating number to square in next step
    NumericVector yd = y[within];

    // change this to loop
    // and then sum sugar vs. c
    dist += sum(sqrt (pow(x[i] - xd, 2) + pow(y[i] - yd, 2))); //calculating Euclidean distance
    nobs += sum(within);

  }

  return dist / nobs;

}


