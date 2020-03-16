## Set up Model Parameter List, this imports "Groups", but also sets parameters that are
## fixed across all groups, or required to run the model
fZooMSS_Params <- function(fileGroups, enviroo, tmax, f_mort){
  
  groups <- fileGroups  # Read in functional group specific parameters from file
  nutrition <- groups$carbon # Extract carbon content (nutr. quality) of each group
  environ <- enviroo # Environmental information
  
  # Set up parameter list
  param <- list(groups = groups,
               environ = environ,
               nutrition = nutrition,
               
               # Model parameters
               ngrps = dim(groups)[1],		# no. of groups
               tmax = tmax,					# no. of years
               dt = environ$dt, 				# time step
               dx = 0.1,         # log10 weight step
               day = 12,          # day length (hours of each day in sun)
               gge_base = 0.25, # baseline gross growth efficiency
               w0 = 10^(min(groups$W0)),		# minimum dynamic size class
               wMax = 10^(max(groups$Wmax)),# maximum dynamic size class
               ZSpre = 1, # senescence mortality prefactor
               ZSexp = 0.3, # senescence mortality exponent
               f_mort = f_mort, # fishing mortality (yr^-1)
               w0_phyto = 10^(-14.5),		# minimum phytoplankton size class (1um)
               wMax_phyto = 10^environ$phyto_max,		# maximum phytoplankton size class
               fixed_filterPPMR = fixed_filterPPMR
  )
  
  param$isave <- 100	# how often to save results every 'isave' time steps
  param$fish_grps <- which(is.na(groups$prop) == TRUE) # Which rows are fish
  param$zoo_grps <- which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  param$num_zoo <- sum(is.na(param$groups$prop) == FALSE) # How many zooplankton
  param$num_fish <- sum(is.na(param$groups$prop) == TRUE) # How many fish
  
  return(param)
}