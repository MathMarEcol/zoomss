
# Function to calculate the mean of the last 50 % of the model
fZooMSS_AveOutput = function(x){
  ave_x <- colMeans(x[(ceiling(0.5*(dim(x)[1])):dim(x)[1]),,], dims = 1)
  return(ave_x)
}

untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense