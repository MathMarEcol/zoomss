fZooMSS_inner_project_loop <- Rcpp::cppFunction('NumericMatrix fZooMSS_inner_project_loop(int no_sp, int no_w,
                NumericMatrix niter, NumericMatrix Aiter, NumericMatrix Citer,
                NumericMatrix Siter, NumericMatrix n, NumericMatrix A, NumericMatrix B, NumericMatrix C,
                NumericMatrix S, NumericVector w_min_idx, NumericVector w_max_idx) {

    for (int i = 0; i < no_sp; i++) {
        for (int j = w_min_idx[i]+1; j < w_max_idx[i]; j++) {
            niter(i,j) = (Siter(i,j) + Aiter(i,j)*niter(i,j-1)) / Citer(i,j);
            if (j > w_min_idx[i]+1) {
            n(i,j - 1) <- (S(i,j - 1) + A(i,j - 1) * n(i,j - 2) + B(i,j - 1) * niter(i,j)) / C(i,j - 1);
            }
        }
        n(i,w_max_idx[i]) = 0;
    }
    return n;
}')