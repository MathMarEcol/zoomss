
r <- rnorm(12*178*178) # GEnerate some random numbers of the correct length for the array
r2 <- rnorm(178) # Random number for the vector


# f1 is the code as we currently run it in ZooMSS.
# Here we create arbitary matrices to test
f1 <- function(r,r2){
  ar <- array(r,c(12, 178, 178)) # Large Array
  ve <- r2 # numeric vector

  sw <- sweep(ar, 3, ve, '*') # n_species x n_sizes x n_sizes # Large Array
  ap <- aperm(sw, c(3,1,2)) # n_sizes x n_species x n_sizes # LArge Array
  cs <- colSums(ap) # n_species x n_sizes # Matrix
  cs
}



# out_f1 <- f1(r,r2)

bench::mark(f1(r,r2),
            f1(r,r2)
            )
