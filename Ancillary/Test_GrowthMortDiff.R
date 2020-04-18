library(bench)
r <- rnorm(12*178*178) # Generate some random numbers of the correct length for the array
r2 <- rnorm(178) # Random numbers for the vector

# f1 is the code as we currently run it in ZooMSS.
# Here we create arbitrary matrices to test
f1 <- function(r, r2) {
  # ar <- array(r,c(12, 178, 178)) # Large Array
  dim(r) <- c(12, 178, 178) # Large Array
  # ve <- r2 # numeric vector
  sw <- sweep(r, 3, r2, '*') # sw will be a Large Array with size 12,178,178
  ap <- aperm(sw, c(3, 1, 2)) # ap is a Large Array with size 178,12,178
  cs <- colSums(ap) # cs is a Matrix of size 12,178
  cs
}

f2 <- function(r, r2) {
  dim(r) <- c(12, 178, 178) # Large Array


  r <- sweep(r, 3, r2, '*', FALSE) # r will be a Large Array with size 12,178,178
  r <- aperm(r, c(3, 1, 2)) # ap is a Large Array with size 178,12,178
  colSums(r) # cs is a Matrix of size 12,178
}

f3 <- function(r, r2) {
  dim(r) <- c(12, 178, 178)
  colSums(
    aperm(
      sweep(r, 3, r2, '*', FALSE),
      c(3, 1, 2)
    )
  )
}

f4 <- function(r, r2) {
  dim(r) <- c(12*178, 178)
  r <- colSums(r2 * t(r))
  dim(r) <- c(12, 178)
  r
}

f5 <- function(r, r2) {
  dim(r) <- c(12, 178, 178) # Large Array
  dim(r) <- c(12*178, 178)
  r <- colSums(r2 * t(r))
  dim(r) <- c(12, 178)
  r
}


# Are they identical?
out_f1 <- f1(r, r2)
out_f2 <- f2(r, r2)
out_f3 <- f3(r, r2)
out_f4 <- f4(r, r2)
out_f5 <- f5(r, r2)

identical(out_f4, out_f5)

identical(out_f1, out_f2) &&
  identical(out_f2, out_f3) &&
  identical(out_f3, out_f4) # Must be TRUE or your answer isn't valid

# Which is faster?
bench::mark(f1(r, r2), f2(r, r2), f3(r, r2), f4(r, r2), iterations = 1000)

rbenchmark::benchmark(f1(r, r2), f2(r, r2), f3(r, r2), f4(r, r2), replications = 1000)

microbenchmark::microbenchmark(f1(r, r2), f2(r, r2), f3(r, r2), f4(r, r2))