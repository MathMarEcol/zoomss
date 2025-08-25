#' McKendrick-von Foerster Equation Solver (Base R Implementation)
#'
#' @title Solve McKendrick-von Foerster equation for size-structured populations
#' @description Solves the McKendrick-von Foerster (MvF) partial differential equation
#'   for size-structured population dynamics using a finite difference approach in base R.
#' @details This function implements the numerical solution to the McKendrick-von Foerster
#'   equation, which describes how populations change across size classes over time.
#'   The equation is solved using an upwind finite difference scheme that handles:
#'   - Growth through size classes (advection term)
#'   - Diffusion between adjacent size classes
#'   - Source and sink terms from feeding and mortality
#'
#'   The function processes each functional group separately and applies boundary conditions
#'   appropriate for size-structured models. The last size class is set to zero abundance
#'   to represent maximum size limits.
#'
#'   This is a core computational component of ZooMSS that updates population abundances
#'   at each time step based on growth, mortality, and reproduction processes.
#'
#' @param ngrps Number of functional groups in the model
#' @param curr_min_size Vector of minimum size class indices for each group
#' @param curr_max_size Vector of maximum size class indices for each group
#' @param A_iter Matrix of advection coefficients for current iteration
#' @param C_iter Matrix of diagonal coefficients for current iteration
#' @param Nb_iter Matrix to store updated abundances for current iteration
#' @param S_iter Matrix of source terms for current iteration
#' @param A Matrix of advection coefficients
#' @param B Matrix of diffusion coefficients
#' @param C Matrix of diagonal coefficients
#' @param Nb Matrix of abundances to be updated
#' @param S Matrix of source terms
#'
#' @return Updated abundance matrix (Nb) with new size-class distributions
#'
#' @examples
#' \dontrun{
#' # This function is typically called internally by zoomss_run
#' # Example shows the structure of parameters needed:
#' ngrps <- 9
#' ngrid <- 100
#' curr_min_size <- c(1, 10, 20, 30, 40, 50, 60, 70, 80)
#' curr_max_size <- c(30, 40, 50, 60, 70, 80, 90, 95, 100)
#'
#' # Initialize coefficient matrices
#' A <- matrix(0, nrow = ngrps, ncol = ngrid)
#' B <- matrix(0, nrow = ngrps, ncol = ngrid)
#' C <- matrix(1, nrow = ngrps, ncol = ngrid)
#' S <- matrix(0, nrow = ngrps, ncol = ngrid)
#' Nb <- matrix(0.1, nrow = ngrps, ncol = ngrid)
#'
#' # Run MvF solver
#' updated_abundances <- zoomss_mvf(ngrps, curr_min_size, curr_max_size,
#'                                        A, C, Nb, S, A, B, C, Nb, S)
#' }
#'
#' @noRd
#'
zoomss_mvf <- function(ngrps, curr_min_size, curr_max_size, A_iter, C_iter, Nb_iter, S_iter, A, B, C, Nb, S){

  for(i in 1:ngrps){

    idx_curr <- (curr_min_size[i]+1):curr_max_size[i] ## Set size range index for current group

    for(j in idx_curr){ ## Find the abundance at the next size class with standard MvF
      Nb_iter[i,j] <- (S_iter[i,j] + A_iter[i,j]*Nb[i,j-1])/(C_iter[i,j])

      if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
        k <- j - 1
        Nb[i,k] <- (S[i,k] + A[i,k] * Nb[i,k-1] + B[i,k] * Nb_iter[i,k+1]) / C[i,k]
      }

      # MvF without diffusion for last size class
      if(j == idx_curr[length(idx_curr)]){
        Nb[i,j] <- 0
      }
    }
  }
  return(Nb)
}
