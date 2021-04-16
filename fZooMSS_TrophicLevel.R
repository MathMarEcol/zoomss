## TROPHIC LEVEL FUNCTION, takes a 12x15 matrix of diets (predators are rows, prey are columns) from ZooMSS and calculates
## trophic levels of predator groups.

fZooMSS_TrophicLevel <- function(diet_matrix){

  phyto_tl <- 1 # Phyto TL is 1
  start_dynam_tl <- rep(2,12) # Start TL - start at 2 for all zoo and fish groups

  # To be truly generic I need to fix these up with testgroup references #TODO
  curr_phyto_diet <- rowSums(diet_matrix[,1:3]) # Current phyto diet
  curr_dynam_diet <- diet_matrix[,4:15] # Current heterotroph diet

  total_diet <- curr_phyto_diet + rowSums(curr_dynam_diet) # Total consumption, in grams wet weight, by pred group

  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
  curr_dynam_frac <- sweep(curr_dynam_diet, 1, total_diet, '/') # Fraction of diet from each prey group for each pred group

  n <- 1
  eps_diff <- 1

  while(eps_diff > 0.01 & n < 100){ # Gauss-Siedel iterative loop to calculate trophic levels, stops when converged or number of loops reaches 100
    n <- n + 1
    eps <- start_dynam_tl[10]

    calc_dynam_tl = sweep(curr_dynam_frac, 2, start_dynam_tl, '*')
    calc_dynam_tl[which(is.nan(calc_dynam_tl) == TRUE)] = 0 # Get rid of nans - these are entrys where there is no biomass for a given group
    #calc_dynam_tl[which(calc_dynam_tl == Inf)] = 0 # Get rid of infinite values, occurs with asymptotic size bins, because there is no biomass to have a diet in those bins
    start_dynam_tl = 1+phyto_tl*curr_phyto_frac + rowSums(calc_dynam_tl) # Update trophic level matrix

    eps_diff <- abs(eps - start_dynam_tl[10])
  } # End Gauss-Siedel loop

  return(start_dynam_tl)

} # End trophic level function