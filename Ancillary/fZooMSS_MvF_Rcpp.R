

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix fZooMSS_MvF_Rcpp(int no_sp, int no_w,
                                         NumericMatrix niter, NumericMatrix Aiter, NumericMatrix Citer,
                                         NumericMatrix Siter, NumericMatrix n, NumericMatrix A, NumericMatrix B, NumericMatrix C,
                                         NumericMatrix S, IntegerVector w_min_idx, IntegerVector w_max_idx) {

    for(int i = 0; i < no_sp; i++){
        niter(i,w_min_idx[i]) = (Siter(i,w_min_idx[i]) + Aiter(i,w_min_idx[i])*niter(i,(w_min_idx[i]-1))) / Citer(i,w_min_idx[i]);
        for(int j = (w_min_idx[i]+1); j < w_max_idx[i]; j++){
            niter(i,j) = (Siter(i,j) + Aiter(i,j)*n(i,j-1)) / Citer(i,j);
            n(i,j - 1) = (S(i,j - 1) + A(i,j - 1) * n(i,j - 2) + B(i,j - 1) * niter(i,j)) / C(i,j - 1);
        }
        n(i,(w_max_idx[i]-1)) = 0;
    }
    return n;
}')
