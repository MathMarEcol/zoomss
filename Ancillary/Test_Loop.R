load("Ancillary/Run_out.RData")
source("fZooMSS_inner_project_loop.R")

N_out <- fZooMSS_inner_project_loop(no_sp = ngrps, no_w = ngrid, niter = N.iter,
                                     Aiter = A.iter, Citer = C.iter, Siter = S.iter,
                                     S = S, n = N, A = A, B = B, C = C,
                                     w_min_idx = curr_min_size, w_max_idx = curr_max_size)


load("Ancillary/Run_out.RData") #reload initial data

for(i in 1:ngrps){

  idx_curr <- (curr_min_size[i]+1):curr_max_size[i] ## Set size range index for current group

  for(j in idx_curr){ ## Find the abundance at the next size class with standard MvF
    N.iter[i,j] <- (S.iter[i,j] + A.iter[i,j]*N[i,j-1])/(C.iter[i,j])
    if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
      k <- j - 1
      N[i,k] <- (S[i,k] + A[i,k] * N[i,k-1] + B[i,k] * N.iter[i,k+1]) / C[i,k]
    }

    # MvF without diffusion for last size class
    if(j == idx_curr[length(idx_curr)]){
      N[i,j] <- 0
      # N[i,curr_min_size] <- N.iter[i,curr_min_size] # Keep starting sizes constant
    }
  }
}

identical(N, N_out)
