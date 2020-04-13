## The Setup function calculates the feeding kernels, dvm (if included), temp effects...
## everything that can be solved before we start iterating through time to solve MvF-D

## Last Updated Tuesday 17th March 2020
##

fZooMSS_Setup <- function(param){

  # Set up dynamic weight grid
  w <- 10^(seq(from = log10(param$w0), to = log10(param$wMax), param$dx))
  ngrid <- length(w) # total number of size classes for zoo and fish

  # Set up phytoplankton size classes
  w_phyto <- 10^(seq(from = log10(param$w0_phyto), to = log10(param$wMax_phyto), param$dx))
  ngridPP <- length(w_phyto) # total number of size classes for phyto

  # Number of time slots to save
  nsave  <- floor(param$tmax/(param$dt*param$isave))

  # Commented out DVM as it is unused at the moment (April 2020)
  # # Diel Vertical Migration - change availability of phyto to zoo
  # # and zoo to fish based on slope of phytoplankton (calculated as
  # # proportion of day searching for food and available for predation)
  # dvm_max <- param$day/24 # maximum proportion of day spent migrating
  # ESD_sizes <- 2*(3/(4*pi)*w)^(1/3) # convert g wet weight to ESD (cm)
  # dvm <- dvm_max*(1.02*(ESD_sizes) - 0.02) # size-dependent amount of time spent away from surface
  # dvm[which(w < 10^-5.4)] <- 0 # Microzoo don't migrate (ESD < 0.02cm)
  # dvm[which(w > 10^-0.3)] <- dvm_max # Macrozoo migrate max time (ESD > 10cm)
  # dvm_mat <- matrix(dvm, nrow = ngrid, ncol = ngrid, byrow = TRUE) # Matrix of dvm, nrow = number of pred size classes
  #
  # # This works out the proportion of time an predator of size w will have access to a prey of size w', for all w and w'
  # dvm_mat <- 1 - dvm_mat
  # dvm_mat[lower.tri(dvm_mat)] <- 0
  # dvm_mat <- t(dvm_mat) + dvm_mat
  # diag(dvm_mat) <- diag(dvm_mat)/2

  # Dynamic prey availability matrix: dim1 is predators, dim2 is predator size classes,
  # dim3 is prey groups, dim 4 is prey size classes.
  dynam_theta <- array(1, dim = c(param$ngrps, ngrid, param$ngrps, ngrid))
  dynam_theta <- sweep(dynam_theta, c(2,4), dvm_mat,"*")

  # Phyto availability matrix: rows are predators, columns are their size classes,
  # entries are time spent feeding on phytoplankton for the size class
  phyto_theta <- matrix(1-dvm, nrow = param$ngrps, ncol = ngrid, byrow = TRUE)

  ### REMOVE DVM
  dynam_theta <- array(1, dim = c(param$ngrps, ngrid, param$ngrps, ngrid))
  phyto_theta <- matrix(1, nrow = param$ngrps, ncol = ngrid, byrow = TRUE)
  ###

  carn_grps <- which(param$Groups$FeedType == 'Carnivore')
  phyto_theta[carn_grps,] <- 0 # Carnivorous groups can't eat phyto

  ## Makes the model object, full of constant functions for model
  model <- list(
    param = param,

    # Phytoplankton abundance
    nPP = 10^(param$phyto_int)*(w_phyto^(param$phyto_slope)), # phyto abundance spectrum

    # Grid parameters
    ngrid = ngrid,
    w = w,
    w_phyto = w_phyto,

    # Group parameters storage
    phyto_growthkernel = array(NA, dim = c(param$ngrps, ngrid, ngridPP)), # predation on phytoplankton
    phyto_diffkernel = array(NA, dim = c(param$ngrps, ngrid, ngridPP)), # diffusion from phytoplankton consumption
    phyto_dietkernel = array(NA, dim = c(param$ngrps, ngrid, ngridPP)), # diet from phytoplankton
    dynam_growthkernel = array(NA, dim = c(param$ngrps, ngrid, ngrid)), # predation on zoo and fish
    dynam_diffkernel = array(NA, dim = c(param$ngrps, ngrid, ngrid)), # diffusion from zoo and fish consumption
    dynam_dietkernel = array(NA, dim = c(param$ngrps, ngrid, ngrid)), # diet from zoo and fish
    dynam_mortkernel = array(NA, dim = c(param$ngrps, ngrid, ngrid)), # mortality from predation on dynamic component
    M_sb = matrix(0, nrow = param$ngrps, ncol = ngrid), # senescence mortality
    fish_mort = matrix(0, nrow = param$ngrps, ncol = ngrid), # fishing mortality

    # Assimilation efficiency and temperature effect storage, by group and size class
    assim_eff = matrix(0, nrow = param$ngrps, ncol = ngrid),
    temp_eff = matrix(0, nrow = param$ngrps, ncol = ngrid),

    #### STORAGE FOR DIET KERNELS
    phyto_dietkernel = array(NA, dim = c(param$ngrps, ngrid, ngridPP)),
    dynam_dietkernel = array(NA, dim = c(param$ngrps, ngrid, ngrid)),

    # Output storage
    N = array(NA, dim = c(nsave, param$ngrps, ngrid)), # dynamic abundance spectrum
    # M2 = array(NA, dim = c(nsave, param$ngrps, ngrid)), # Predation mortality
    Z = array(NA, dim = c(nsave, param$ngrps, ngrid)), # Total mortality
    gg = array(NA, dim = c(nsave, param$ngrps, ngrid)), # Growth
    diet = array(NA, dim = c(nsave, c(param$ngrps), c(param$ngrps+3))) # diet
    # Biomass = matrix(NA, nrow = nsave, ncol = param$ngrps), # Biomass of each group
    # Abundance = matrix(NA, nrow = nsave, ncol = param$ngrps), # Abundance of each group
    # Diff = array(NA, dim = c(nsave, param$ngrps, ngrid)) # save diffusion
  )

  # GGE for different groups
  assim_phyto <- (param$Groups$GrossGEscale)*param$cc_phyto # Phytoplankton
  # assim_dynam <- matrix(param$Groups$GrossGEscale*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE) # rows are predators, columns are prey
  model$assim_eff = matrix(param$Groups$GrossGEscale*param$nutrition, nrow = param$ngrps, ncol = length(model$w))

  #### INITIAL DYNAMIC POPULATION ABUNDANCES
  a_dynam <- 10^(param$phyto_int)*(w[1]^(param$phyto_slope+1)) # calculate coefficient for initial dynamic spectrum, so that N(w_phyto) equals
  # N(w_dynam) at w[1]

  # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
  tempN <- matrix(a_dynam*(w)^-1, nrow = param$ngrps, ncol = ngrid, byrow = TRUE)
  props_z <- param$Groups$Prop[param$zoo_grps] # Zooplankton proportions
  tempN[param$zoo_grps,] <- props_z * tempN[param$zoo_grps,] # Set abundances of diff zoo groups based on smallest size class proportions
  tempN[param$fish_grps,] <- (1/param$num_fish) * tempN[param$fish_grps,] # Set abundandances of fish groups based on smallest size class proportions

  # For each group, set densities at w > Winf and w < Wmin to 0
  tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Winf) Winf < wx, Winf = (param$Groups$Wmax)))] <- 0
  tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Wmin) Wmin > wx, Wmin = (param$Groups$W0)))] <- 0
  model$N[1,,] <- tempN

  # Fishing mortality
  model$fish_mort[param$fish_grps, c(w >= 1)] <- param$f_mort

  ### MATRICES FOR LOG TRANSFORM OF EQUATION
  # Predators are rows, phyto prey weights are columns
  gg_log_t_phyto <- ((w^-1) %*% t(w_phyto))/log(10) # Growth
  diff_log_t_phyto <- ((w^-2) %*% t(w_phyto^2))/log(10) # Diffusion
  diet_log_t_phyto <- matrix(w_phyto, nrow = length(w), ncol = length(w_phyto), byrow = TRUE) # Diet/Ingestion

  # Predators are rows, dynam prey weights are columns
  gg_log_t_dynam <- ((w^-1) %*% t(w))/log(10) # Growth
  diff_log_t_dynam <- ((w^-2) %*% t(w^2))/log(10) # Diffusion
  diet_log_t_dynam <- matrix(w, nrow = length(w), ncol = length(w), byrow = TRUE) # Diet/ingestion

  ### PREDATION KERNELS FOR PHYTOPLANKTON SPECTRUM AND DYNAMIC SPECTRUM
  phyto_pred_weight_matrix <- matrix(w, nrow = ngrid, ncol = ngridPP)
  dynam_pred_weight_matrix <- matrix(w, nrow = ngrid, ncol = ngrid)
  phyto_prey_weight_matrix <- matrix(w_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)
  dynam_prey_weight_matrix <- matrix(w, nrow = ngrid, ncol = ngrid, byrow = TRUE)

  ## Search Volume storage
  SearchVol <- matrix(NA, nrow = param$ngrps, ncol = ngrid) # Search volume

  # Simpson's Rule matrices for growth, diffusion and mortality integrals
  simp_phyto <- array(1, dim = ngridPP)
  simp_phyto[c(seq(2, ngridPP-1,2))] <- 4
  simp_phyto[c(seq(3, ngridPP-1,2))] <- 2
  sm_phyto <- matrix(simp_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE) * (param$dx/3)

  simp_dynam <- array(1, dim = ngrid)
  simp_dynam[c(seq(2, ngrid-1,2))] <- 4
  simp_dynam[c(seq(3, ngrid-1,2))] <- 2
  sm_dynam <- matrix(simp_dynam, nrow = ngrid, ncol = ngrid, byrow = TRUE) * (param$dx/3)

  ## Temperature Effect Matrix
  # Effect of temperature on feeding and predation rate
  #temp_flag <- 2.4^((environ$sst)/10)
  #temp_cil <- 2.8^((environ$sst)/10)
  #temp_gel <- 1.75^((environ$sst)/10)
  #temp_larv <- 2.2^((environ$sst)/10)
  #temp_chaet <- 2.44^((environ$sst)/10)
  #temp_crust <- 2.57^((environ$sst)/10)
  #temp_ocop <- 2.8^((environ$sst)/10)
  #temp_ccop <- 2.8^((environ$sst)/10)
  #temp_fish <- 2.6^((environ$sst)/10)

  #temp_zoo <- c(temp_flag, temp_cil, temp_larv, temp_ocop, temp_ccop, temp_crust, temp_chaet, temp_larv, temp_gel)
  #temp_fish <- rep(temp_fish, num_fish)

  #temp_zoo <- rep(exp(23.93 - 0.59/(8.62e-05*(273+environ$sst))), num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
  #temp_fish <- rep(exp(23.93 - 0.59/(8.62e-05*(273+environ$sst))), num_fish) # exp(25.55 - 0.63/(8.62e-05*(273+environ$sst)))/exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))

  ### Q10 OF 2 FOR ALL ZOO AND FISH
  temp_zoo <- rep(2.^((param$sst - 30)/10), param$num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
  temp_fish <- rep(2.^((param$sst - 30)/10), param$num_fish)
  model$temp_eff <- matrix(c(temp_zoo, temp_fish), nrow = param$ngrps, ncol = ngrid)

  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  for(i in 1:param$ngrps){
    ## Senescence mortality
    if(i < 10){
      model$M_sb[i,] <- param$ZSpre*(w/(10^(param$Groups$Wmat[i])))^param$ZSexp
      model$M_sb[i, 10^(param$Groups$Wmax[i]) < w] <- 0
      model$M_sb[i, 10^(param$Groups$Wmat[i]) > w] <- 0
    }

    if(i > 9){
      model$M_sb[i,] <- 0.1*param$ZSpre*(w/(10^(param$Groups$Wmat[i])))^param$ZSexp
      model$M_sb[i, 10^(param$Groups$Wmax[i]) < w] <- 0
      model$M_sb[i, 10^(param$Groups$Wmat[i]) > w] <- 0
    }

    ### Search volume
    SearchVol[i,] <- (param$Groups$SearchCoef[i])*(w^(param$Groups$SearchExp[i]))
    SearchVol[i, 10^(param$Groups$Wmax[i]) < w] <- 0
    SearchVol[i, 10^(param$Groups$W0[i]) > w] <- 0

    ### Predation Kernels
    if(is.na(param$Groups$PPMRscale[i]) == FALSE){ # If group has an m-value (zooplankton)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
      D.z <- 2*(3*w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas <- (exp(0.02*log(D.z)^2 - param$Groups$PPMRscale[i] + 1.832))^3 # Wirtz's equation
      beta_mat_phyto <- matrix(betas, nrow = ngrid, ncol = ngridPP)
      beta_mat_dynam <- matrix(betas, nrow = ngrid, ncol = ngrid)

      # Calculate feeding kernels
      sp_phyto_predkernel <- exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                            phyto_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)
      sp_dynam_predkernel <- exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                            dynam_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)

      # The feeding kernal of filter feeders is not expected to change much with increasing size so we fix it here

      # if (param$fixed_filterPPMR == TRUE){
        if(i == 3){
          sp_phyto_predkernel <- matrix(sp_phyto_predkernel[44,], nrow = ngrid, ncol = ngridPP, byrow = TRUE)
          sp_dynam_predkernel <- matrix(sp_dynam_predkernel[44,], nrow = ngrid, ncol = ngrid, byrow = TRUE)
        }
        if(i == 8){
          sp_phyto_predkernel <- matrix(sp_phyto_predkernel[61,], nrow = ngrid, ncol = ngridPP, byrow = TRUE)
          sp_dynam_predkernel <- matrix(sp_dynam_predkernel[61,], nrow = ngrid, ncol = ngrid, byrow = TRUE)
        }
      # }

    } else { # If group does not have an m-value (fish)
      beta_mat_phyto <- matrix(param$Groups$PPMR[i], nrow = ngrid, ncol = ngridPP)
      beta_mat_dynam <- matrix(param$Groups$PPMR[i], nrow = ngrid, ncol = ngrid)

      # Calculate feeding kernels
      sp_phyto_predkernel <- exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                            phyto_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)
      sp_dynam_predkernel <- exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                            dynam_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)
    }

    ### GROWTH INTEGRAL CONSTANTS
    # Predators are rows, prey are columns
    model$phyto_growthkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP) *
      sp_phyto_predkernel * gg_log_t_phyto * sm_phyto
    model$dynam_growthkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*gg_log_t_dynam*sm_dynam

    ### DIET INTEGRAL CONSTANTS
    # Predators are rows, prey are columns
    model$phyto_dietkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
      sp_phyto_predkernel*diet_log_t_phyto*sm_phyto
    model$dynam_dietkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*diet_log_t_dynam*sm_dynam

    ### DIFFUSION INTEGRAL CONSTANTS
    # Predators are rows, prey are columns
    model$phyto_diffkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
      sp_phyto_predkernel*diff_log_t_phyto*sm_phyto
    model$dynam_diffkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*diff_log_t_dynam*sm_dynam

    ### MORTALITY INTEGRAL CONSTANTS
    # Prey are rows, predators are columns
    model$dynam_mortkernel[i,,] <- matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid, byrow = TRUE)*
      t(sp_dynam_predkernel)*sm_dynam
  }

  # no_sen = which(param$Groups$species == c("Flagellates", "Ciliates")) # no senescence mortality for flagellates and ciliates
  #model$M_sb[c(ngrps)] = 0
  model$M_sb <- model$temp_eff * model$M_sb # Incorporate temp effect on senscence mortality

  ## Incorporate carnivory (groups that can't eat phyto), temperature effects and gross growth efficiency (assim)
  model$phyto_growthkernel <- sweep(sweep(model$phyto_growthkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto, "*")
  model$phyto_diffkernel <- sweep(sweep(model$phyto_diffkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto^2, "*")
  model$phyto_dietkernel <- sweep(sweep(model$phyto_dietkernel, c(1,2), phyto_theta, "*"), 1, 1, "*")

  # Dim 1 = pred group, dim2 = pred sizes, dim 3 = prey sizes
  model$dynam_growthkernel <- sweep(model$dynam_growthkernel, c(1,2), model$temp_eff, '*')
  model$dynam_diffkernel <- sweep(model$dynam_diffkernel, c(1,2), model$temp_eff^2, '*')

  # We still need four dimensions for diet matrix
  model$dynam_dietkernel <- sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_dietkernel, "*"),
                                       c(1,3), 1, "*"), c(1,2), model$temp_eff, "*")

  # We won't sweep through temp_effect here, but do it in fZooMSS_Run function. This is to make sure
  # it can still work in the future if we have different temperature effects for different groups,
  # so then dynam_mortkernel is dim1 = prey sizes, dim 2 = pred groups, dim 3 = pred size
  model$dynam_mortkernel <- aperm(model$dynam_mortkernel, c(2,1,3))

  #### Because phyto spectrum is constant, we can solve the phyto component of growth, and diffusion before time loop
  model$ingested_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2)) # Ingested phyto
  model$diff_phyto <- model$temp_eff^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2)) # Diffusion from phyto
  model$diet_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP, "*"), dims = 2)) # Diet of total phyto

  model$diet_phyto_all_sizes <- sweep(sweep(model$phyto_dietkernel, 3, model$nPP, "*"),  c(1,2), model$temp_eff, "*") # Diet of total phyto, with all size classes of phyto maintained

  ## Diet of phyto from pico, nano and micro size classes
  model$diet_pico_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) < -11.5), "*"), dims = 2))
  model$diet_nano_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -11.5 & log10(model$w_phyto) < -8.5), "*"), dims = 2))
  model$diet_micro_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -8.5), "*"), dims = 2))

  return(model)
} # End of Setup function