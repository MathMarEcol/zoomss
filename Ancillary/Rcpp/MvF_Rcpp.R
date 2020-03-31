
library(Rcpp)

MvF_Rcpp <- Rcpp::cppFunction('NumericMatrix MvF_Rcpp(int no_sp, int no_w,
                NumericMatrix n, NumericMatrix A, NumericMatrix B,
                NumericMatrix S, NumericVector w_min_idx, NumericVector w_max_idx) {

    for (int i = 0; i < no_sp; i++) {
        for (int j = w_min_idx[i]+1; j < w_max_idx[i]; j++) {
            n(i,j) = (S(i,j) + A(i,j)*n(i,j-1)) / B(i,j);
        }
        n(i,w_max_idx[i]) = 0;
    }
    return n;
}')
