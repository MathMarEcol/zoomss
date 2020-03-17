## Set up Model Parameter List, this imports "Groups", but also sets parameters that are
## fixed across all Groups, or required to run the model

## Last Updated Tuesday 17th March 2020
##

fZooMSS_Params <- function(Groups, input_params){

  param <- list(Groups = Groups, # Read in functional group specific parameters from file
                nutrition = Groups$Carbon, # Extract carbon content (nutr. quality) of each group
                ngrps = dim(Groups)[1],		# no. of Groups
                dx = 0.1,         # log10 weight step
                day = 12,          # day length (hours of each day in sun)
                gge_base = 0.25, # baseline gross growth efficiency
                w0 = 10^(min(Groups$W0)),		# minimum dynamic size class
                wMax = 10^(max(Groups$Wmax)),# maximum dynamic size class
                ZSpre = 1, # senescence mortality prefactor
                ZSexp = 0.3, # senescence mortality exponent
                f_mort = 0, # fishing mortality (yr^-1)
                w0_phyto = 10^(-14.5),		# minimum phytoplankton size class (1um)
                wMax_phyto = 10^input_params$phyto_max,		# maximum phytoplankton size class
                cc_phyto = 0.1,  # Carbon content of phytoplankton size classes
                fish_on = TRUE,
                isave = 100,	# how often to save results every 'isave' time steps
                zoo_grps = which(Groups$Type == "Zooplankton"), # Which rows are zooplankton
                fish_grps = which(Groups$Type == "Fish"), # Which rows are fish
                num_zoo = sum(Groups$Type == "Zooplankton"), # How many zooplankton
                num_fish = sum(Groups$Type == "Fish") # How many fish
  )
  param <- c(input_params, param) # Join with input_params
  return(param)
}
