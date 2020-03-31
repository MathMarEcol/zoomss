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

NumericMatrix inner_project_loop(int no_sp, int no_w,
                                 NumericMatrix n, 
                                 NumericMatrix A, 
                                 NumericMatrix B,
                                 NumericMatrix S, 
                                 NumericVector w_min_idx, 
                                 NumericVector w_max_idx) {
  
  for (int i = 0; i < no_sp; i++) {
    for (int j = w_min_idx[i]+1; j < w_max_idx[i]; j++) {
      n(i,j) = (S(i,j) + A(i,j)*n(i,j-1)) / B(i,j);
    }
    n(i,w_max_idx[i]) = 0;
  }
}
