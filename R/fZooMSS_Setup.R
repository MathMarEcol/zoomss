#' Setup ZooMSS Model Structure and Feeding Kernels
#'
#' @title Initialize ZooMSS model components and calculate feeding interactions
#' @description Sets up the ZooMSS model structure by calculating feeding kernels, mortality 
#'   rates, and other model components that remain static during the simulation.
#' @details This function initializes the core ZooMSS model structure by calculating:
#'   
#'   **Static Components (calculated once):**
#'   - Feeding preference kernels based on predator-prey size ratios
#'   - Search volumes and encounter rates between size classes
#'   - Baseline mortality rates (senescence, fishing)
#'   - Initial abundance distributions for all functional groups
#'   
#'   **Dynamic Component Structures (updated during run):**
#'   - Phytoplankton feeding kernels (structure calculated here, values updated with environment)
#'   - Growth and diffusion kernels for zooplankton and fish interactions
#'   - Diet and mortality tracking arrays
#'   
#'   **Model Architecture:**
#'   - Size-structured populations across logarithmic size classes  
#'   - Multiple functional groups with different feeding behaviors
#'   - Environmental coupling through phytoplankton and temperature
#'   
#'   The function separates static calculations (done once for efficiency) from 
#'   dynamic calculations (updated each time step in fZooMSS_Run).
#'
#' @param param Complete parameter list created by fZooMSS_Params containing:
#'   - Groups: Functional group definitions and biological parameters
#'   - Model dimensions (ngrps, ngrid, time parameters)  
#'   - Environmental forcing time series
#'   - Physical and biological constants
#'
#' @return Model object containing:
#'   \itemize{
#'     \item param: Input parameters (passed through)
#'     \item dynam_*: Dynamic feeding kernel arrays for group interactions
#'     \item phyto_*: Phytoplankton feeding kernel arrays
#'     \item nPP: Initial phytoplankton abundance spectrum
#'     \item M_sb_base: Baseline senescence mortality rates
#'     \item fish_mort: Fishing mortality rates
#'     \item assim_eff: Assimilation efficiency matrix
#'     \item temp_eff: Temperature effect matrix (initialized)
#'     \item N: Initial abundance arrays
#'     \item Additional model structure components
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Create parameters for model setup
#' params <- fZooMSS_Params(Groups, input_params)
#' 
#' # Initialize model structure
#' model <- fZooMSS_Setup(params)
#' 
#' # Model is now ready for time integration with fZooMSS_Run
#' results <- fZooMSS_Run(model)
#' }
#'
fZooMSS_Setup <- function(param){

  ## Dynamic prey availability matrix: dim1 is predators, dim2 is predator size classes,
  ## dim3 is prey groups, dim 4 is prey size classes.
  dynam_theta <- array(1, dim = c(param$ngrps, param$ngrid, param$ngrps, param$ngrid))


  ## Makes the model object, full of constant functions for model
  model <- list(
    param = param,
    
    # Dynamic component kernels (these don't change with environment, only with abundances)
    dynam_growthkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngrid)), # predation on zoo and fish
    dynam_diffkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngrid)), # diffusion from zoo and fish consumption
    dynam_dietkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngrid)), # diet from zoo and fish
    dynam_mortkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngrid)), # mortality from predation on dynamic component
    
    # Phytoplankton kernels (these don't change with environment, only with phytoplankton spectrum)
    phyto_growthkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngridPP)), # predation on phytoplankton
    phyto_diffkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngridPP)), # diffusion from phytoplankton consumption
    phyto_dietkernel = array(NA, dim = c(param$ngrps, param$ngrid, param$ngridPP)), # diet from phytoplankton
    
    # Phytoplankton spectrum - will be updated dynamically with current time step parameters
    nPP = 10^(param$phyto_int_ts[1])*(param$w_phyto^(param$phyto_slope_ts[1])), # Phytoplankton abundance spectrum (dynamic)
    
    # Static mortality and group parameters
    M_sb_base = matrix(0, nrow = param$ngrps, ncol = param$ngrid), # base senescence mortality (before temp effect)
    fish_mort = matrix(0, nrow = param$ngrps, ncol = param$ngrid), # fishing mortality

    # Assimilation efficiency (constant)
    assim_eff = matrix(param$Groups$GrossGEscale * param$Groups$Carbon, nrow = param$ngrps, ncol = length(param$w)),

    # Temperature effects matrix - initialize with first timestep values
    temp_eff = matrix(1, nrow = param$ngrps, ncol = param$ngrid), # Will be updated dynamically in run
    M_sb = matrix(0, nrow = param$ngrps, ncol = param$ngrid), # senescence mortality with temp effect

    # Phytoplankton feeding parameters (for maintaining compatibility)
    phyto_theta = matrix(1, nrow = param$ngrps, ncol = param$ngrid, byrow = TRUE),
    
    # Abundance storage for model - needed for Run function
    N = array(NA, dim = c(param$nsave, param$ngrps, param$ngrid)), # dynamic abundance spectrum
    Z = array(NA, dim = c(param$nsave, param$ngrps, param$ngrid)), # Total mortality  
    gg = array(NA, dim = c(param$nsave, param$ngrps, param$ngrid)), # Growth
    diet = array(NA, dim = c(param$nsave, c(param$ngrps), c(param$ngrps+3))) # diet
  )

  # Set phyto_theta for carnivores
  model$phyto_theta[which(param$Groups$FeedType == 'Carnivore'),] <- 0 # Carnivorous groups can't eat phyto

  # GGE for different groups
  assim_phyto <- (param$Groups$GrossGEscale) * param$cc_phyto # Phytoplankton

  #### INITIAL DYNAMIC POPULATION ABUNDANCES
  # Use the first time step for initial conditions
  a_dynam <- 10^(param$phyto_int_ts[1])*(param$w[1]^(param$phyto_slope_ts[1]+1)) # calculate coefficient for initial dynamic spectrum

  # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
  tempN <- matrix(a_dynam*(param$w)^-1, nrow = param$ngrps, ncol = param$ngrid, byrow = TRUE)
  props_z <- param$Groups$Prop[param$zoo_grps] # Zooplankton proportions
  tempN[param$zoo_grps,] <- props_z * tempN[param$zoo_grps,] # Set abundances of diff zoo groups based on smallest size class proportions
  tempN[param$fish_grps,] <- (1/param$num_fish) * tempN[param$fish_grps,] # Set abundandances of fish groups based on smallest size class proportions

  # For each group, set densities at w > Winf and w < Wmin to 0
  tempN[unlist(tapply(param$w_log10, seq_along(param$w), function(wx,Winf) Winf < wx, Winf = param$Groups$Wmax))] <- 0
  tempN[unlist(tapply(param$w_log10, seq_along(param$w), function(wx,Wmin) Wmin > wx, Wmin = param$Groups$W0))] <- 0
  model$N[1,,] <- tempN

  # Fishing mortality (yr^-1)
  for(g in 1:12){
    model$fish_mort[g,match(param$Groups$Fmort_W0[g], param$w_log10):match(param$Groups$Fmort_Wmax[g], param$w_log10)] <- param$Groups$Fmort[g]
  }

  ### MATRICES FOR LOG TRANSFORM OF EQUATION (these don't change with environment)
  # Dynamic component matrices
  gg_log_t_dynam <- ((param$w^-1) %*% t(param$w))/log(10) # Growth
  diff_log_t_dynam <- ((param$w^-2) %*% t(param$w^2))/log(10) # Diffusion
  diet_log_t_dynam <- matrix(param$w, nrow = length(param$w), ncol = length(param$w), byrow = TRUE) # Diet/ingestion

  # Phytoplankton component matrices
  gg_log_t_phyto <- ((param$w^-1) %*% t(param$w_phyto))/log(10) # Growth from phytoplankton
  diff_log_t_phyto <- ((param$w^-2) %*% t(param$w_phyto^2))/log(10) # Diffusion from phytoplankton
  diet_log_t_phyto <- matrix(param$w_phyto, nrow = length(param$w), ncol = length(param$w_phyto), byrow = TRUE) # Diet from phytoplankton

  ### PREDATION KERNELS FOR DYNAMIC SPECTRUM (constant across time)
  dynam_pred_weight_matrix <- matrix(param$w, nrow = param$ngrid, ncol = param$ngrid)
  dynam_prey_weight_matrix <- matrix(param$w, nrow = param$ngrid, ncol = param$ngrid, byrow = TRUE)

  ### PREDATION KERNELS FOR PHYTOPLANKTON (constant across time)  
  phyto_pred_weight_matrix <- matrix(param$w, nrow = param$ngrid, ncol = param$ngridPP)
  phyto_prey_weight_matrix <- matrix(param$w_phyto, nrow = param$ngrid, ncol = param$ngridPP, byrow = TRUE)

  ## Search Volume storage
  SearchVol <- matrix(NA, nrow = param$ngrps, ncol = param$ngrid) # Search volume

  # Simpson's Rule matrices for growth, diffusion and mortality integrals
  simp_dynam <- array(1, dim = param$ngrid)
  simp_dynam[c(seq(2, param$ngrid-1,2))] <- 4
  simp_dynam[c(seq(3, param$ngrid-1,2))] <- 2
  sm_dynam <- matrix(simp_dynam, nrow = param$ngrid, ncol = param$ngrid, byrow = TRUE) * (param$dx/3)

  # Simpson's Rule for phytoplankton 
  simp_phyto <- array(1, dim = param$ngridPP)
  simp_phyto[c(seq(2, param$ngridPP-1,2))] <- 4
  simp_phyto[c(seq(3, param$ngridPP-1,2))] <- 2
  sm_phyto <- matrix(simp_phyto, nrow = param$ngrid, ncol = param$ngridPP, byrow = TRUE) * (param$dx/3)

  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  for(i in 1:param$ngrps){
    ## Base senescence mortality (before temperature effect)
    if(param$Groups$Type[i] == "Zooplankton"){
      model$M_sb_base[i,] <- param$ZSpre*(param$w/(10^(param$Groups$Wmat[i])))^param$ZSexp
      model$M_sb_base[i, 10^(param$Groups$Wmax[i]) < param$w] <- 0
      model$M_sb_base[i, 10^(param$Groups$Wmat[i]) > param$w] <- 0
    }

    if(param$Groups$Type[i] == "Fish"){
      model$M_sb_base[i,] <- 0.1*param$ZSpre*(param$w/(10^(param$Groups$Wmat[i])))^param$ZSexp
      model$M_sb_base[i, 10^(param$Groups$Wmax[i]) < param$w] <- 0
      model$M_sb_base[i, 10^(param$Groups$Wmat[i]) > param$w] <- 0
    }

    ### Search volume
    SearchVol[i,] <- (param$Groups$SearchCoef[i])*(param$w^(param$Groups$SearchExp[i]))
    SearchVol[i, 10^(param$Groups$Wmax[i]) < param$w] <- 0
    SearchVol[i, 10^(param$Groups$W0[i]) > param$w] <- 0

    ### Predation Kernels for dynamic spectrum (these don't change with environment)
    if(!is.na(param$Groups$PPMRscale[i])){ # If group has a PPMR scaling value (m-value)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
      D.z <- 2*(3*param$w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas <- (exp(0.02*log(D.z)^2 - param$Groups$PPMRscale[i] + 1.832))^3 # Wirtz's equation
      beta_mat_dynam <- matrix(betas, nrow = param$ngrid, ncol = param$ngrid)

      # Calculate feeding kernels
      sp_dynam_predkernel <- exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                             dynam_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)

      # The feeding kernel of filter feeders is not expected to change much with increasing size so we fix it here
      if(param$Groups$FeedType[i] == "FilterFeeder"){
        w0idx <- which(param$Groups$W0[i] == param$w_log10)
        sp_dynam_predkernel <- matrix(sp_dynam_predkernel[w0idx,], nrow = param$ngrid, ncol = param$ngrid, byrow = TRUE)
        rm(w0idx)
      }

    } else { # If group is fish
      beta_mat_dynam <- matrix(param$Groups$PPMR[i], nrow = param$ngrid, ncol = param$ngrid)

      # Calculate feeding kernels
      sp_dynam_predkernel <- exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                             dynam_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
        sqrt(2*pi*param$Groups$FeedWidth[i]^2)
    }

    ### PHYTOPLANKTON FEEDING KERNELS  
    # Calculate phytoplankton predation kernel - consistent with model PPMR scaling
    if(!is.na(param$Groups$PPMRscale[i])){ # If group has a PPMR scaling value (m-value)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
      D.z <- 2*(3*param$w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas <- (exp(0.02*log(D.z)^2 - param$Groups$PPMRscale[i] + 1.832))^3 # Wirtz's equation
      beta_mat_phyto <- matrix(betas, nrow = param$ngrid, ncol = param$ngridPP)
    } else { # If group is fish
      beta_mat_phyto <- matrix(param$Groups$PPMR[i], nrow = param$ngrid, ncol = param$ngridPP)
    }
    
    sp_phyto_predkernel <- exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                           phyto_pred_weight_matrix)/param$Groups$FeedWidth[i])^2)/
      sqrt(2*pi*param$Groups$FeedWidth[i]^2)

    # For filter feeders, fix the phytoplankton feeding kernel too
    if(param$Groups$FeedType[i] == "FilterFeeder"){
      w0idx <- which(param$Groups$W0[i] == param$w_log10)
      sp_phyto_predkernel <- matrix(sp_phyto_predkernel[w0idx,], nrow = param$ngrid, ncol = param$ngridPP, byrow = TRUE)
    }

    ### GROWTH INTEGRAL CONSTANTS FOR PHYTOPLANKTON
    model$phyto_growthkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngridPP) *
      sp_phyto_predkernel * gg_log_t_phyto * sm_phyto

    ### DIET INTEGRAL CONSTANTS FOR PHYTOPLANKTON
    model$phyto_dietkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngridPP)*
      sp_phyto_predkernel*diet_log_t_phyto*sm_phyto

    ### DIFFUSION INTEGRAL CONSTANTS FOR PHYTOPLANKTON
    model$phyto_diffkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngridPP)*
      sp_phyto_predkernel*diff_log_t_phyto*sm_phyto

    ### GROWTH INTEGRAL CONSTANTS FOR DYNAMIC SPECTRUM
    model$dynam_growthkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngrid)*
      sp_dynam_predkernel*gg_log_t_dynam*sm_dynam

    ### DIET INTEGRAL CONSTANTS FOR DYNAMIC SPECTRUM  
    model$dynam_dietkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngrid)*
      sp_dynam_predkernel*diet_log_t_dynam*sm_dynam

    ### DIFFUSION INTEGRAL CONSTANTS FOR DYNAMIC SPECTRUM
    model$dynam_diffkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngrid)*
      sp_dynam_predkernel*diff_log_t_dynam*sm_dynam

    ### MORTALITY INTEGRAL CONSTANTS FOR DYNAMIC SPECTRUM
    # Prey are rows, predators are columns
    model$dynam_mortkernel[i,,] <- matrix(SearchVol[i,], nrow = param$ngrid, ncol = param$ngrid, byrow = TRUE)*
      t(sp_dynam_predkernel)*sm_dynam
  }

  # Transform mortkernel dimensions to match original model expectations
  model$dynam_mortkernel <- aperm(model$dynam_mortkernel, c(2,1,3))

  ## Incorporate carnivory (groups that can't eat phyto), temperature effects and gross growth efficiency (assim)
  model$phyto_growthkernel <- sweep(sweep(model$phyto_growthkernel, c(1,2), model$phyto_theta, "*"), 1, assim_phyto, "*")
  model$phyto_diffkernel <- sweep(sweep(model$phyto_diffkernel, c(1,2), model$phyto_theta, "*"), 1, assim_phyto^2, "*")
  model$phyto_dietkernel <- sweep(sweep(model$phyto_dietkernel, c(1,2), model$phyto_theta, "*"), 1, 1, "*")

  ## Initialize M_sb with temperature effects applied (consistent with model structure)
  model$M_sb <- model$temp_eff * model$M_sb_base # Incorporate temp effect on senescence mortality

  ## Convert diet kernel to 4D for proper calculation (consistent with model structure)
  # We need four dimensions for diet matrix: pred groups x pred sizes x prey groups x prey sizes
  model$dynam_dietkernel <- sweep(dynam_theta, c(1,2,4), model$dynam_dietkernel, "*")

  return(model)
} # End of Setup function