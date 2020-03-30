#library(mizer)
# inner_project_loop <- function(no_sp, no_w, n, A, B, S, w_min_idx) {
#   .Call('_mizer_inner_project_loop', PACKAGE = 'mizer', no_sp, no_w, n, A, B, S, w_min_idx)
# }
library(Rcpp)

#multi_zoo_slurm_filter <- function(sst, chlo, a, b, phyto_max, dt, tmaxx){
multi_zoo_slurm_filter <- function(enviro_vector){

  inner_project_loop <- Rcpp::cppFunction('NumericMatrix inner_project_loop(int no_sp, int no_w,
                NumericMatrix n, NumericMatrix A, NumericMatrix B,
                NumericMatrix S, NumericVector w_min_idx, NumericVector w_max_idx) {

    for (int i = 0; i < no_sp; i++) {
        for (int j = w_min_idx[i]+1; j < w_max_idx[i]; j++) {
            n(i,j) = (S(i,j) + A(i,j)*n(i,j-1)) / B(i,j);
        }
        n(i,w_max_idx[i]) = 0;
    }
    return n;
}')
  ## Set up Model Parameter List, this imports "Groups", but also sets parameters that are
  ## fixed across all groups, or required to run the model
  params <- function(fileGroups, enviroo, tmax, f_mort){

    groups = fileGroups  # Read in functional group specific parameters from file
    nutrition = groups$carbon # Extract carbon content (nutr. quality) of each group
    environ = enviroo # Environmental information

    # Set up parameter list
    param = list(groups = groups,
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
                 wMax_phyto = 10^environ$phyto_max		# maximum phytoplankton size class
    )

    param$isave = 1	# how often to save results every 'isave' time steps
    param$fish_grps = which(is.na(groups$prop) == TRUE) # Which rows are fish
    param$zoo_grps = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
    param$num_zoo = sum(is.na(param$groups$prop) == FALSE) # How many zooplankton
    param$num_fish = sum(is.na(param$groups$prop) == TRUE) # How many fish
    return(param)
  }

  ## The Setup function calculates the feeding kernels, dvm (if included), temp effects...
  ## everything that can be solved before we start iterating through time to solve MvF-D
  Setup = function(param){

    # Pull out some useful parameters - just a shortcut
    grp = param$groups # Functional group info table
    ngrps = param$ngrps # number of functional groups
    dt = param$dt # time step
    dx = param$dx # log10 weight step
    environ = param$environ # environmental info
    fish_grps = param$fish_grps # which rows of grp are fish
    zoo_grps = param$zoo_grps # which rows of grp are zoo
    num_zoo = param$num_zoo # how many zooplankton groups
    num_fish = param$num_fish # how many fish groups

    # Set up dynamic weight grid
    w = 10^(seq(from = log10(param$w0), to =  log10(param$wMax), dx))
    ngrid = length(w) # total number of size classes for zoo and fish

    # Set up phytoplankton size classes
    w_phyto = 10^(seq(from = log10(param$w0_phyto), to = log10(param$wMax_phyto), dx))
    ngridPP = length(w_phyto) # total number of size classes for phyto

    # Number of time slots to save
    nsave   = floor(param$tmax/(param$dt*param$isave))


    # Diel Vertical Migration - change availability of phyto to zoo
    # and zoo to fish based on slope of phytoplankton (calculated as
    # proportion of day searching for food and available for predation)
    dvm_max = param$day/24 # maximum proportion of day spent migrating
    ESD_sizes = 2*(3/(4*pi)*w)^(1/3) # convert g wet weight to ESD (cm)
    dvm = dvm_max*(1.02*(ESD_sizes) - 0.02) # size-dependent amount of time spent away from surface
    dvm[which(w < 10^-5.4)] = 0 # Microzoo don't migrate (ESD < 0.02cm)
    dvm[which(w > 10^-0.3)] = dvm_max # Macrozoo migrate max time (ESD > 10cm)

    #plot(log10(w), dvm, ylab = "Proportion of time spent DVMing (away)") # If you want to have a look at how dvm is parameterised across body size

    dvm_mat = matrix(dvm, nrow = ngrid, ncol = ngrid, byrow = TRUE) # Matrix of dvm, nrow = number of pred size classes
    # ncol = number of pre size classes

    # This works out the proportion of time an predator of size w will have access to a prey of size w', for all w and w'
    dvm_mat = 1 -  dvm_mat
    dvm_mat[lower.tri(dvm_mat)] = 0
    dvm_mat = t(dvm_mat) + dvm_mat
    diag(dvm_mat) = diag(dvm_mat)/2

    # Dynamic prey availability matrix: dim1 is predators, dim2 is predator size classes,
    # dim3 is prey groups, dim 4 is prey size classes.
    dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
    dynam_theta = sweep(dynam_theta, c(2,4), dvm_mat,"*")

    # Phyto availability matrix: rows are predators, columns are their size classes,
    # entries are time spent feeding on phytoplankton for the size class
    phyto_theta = matrix(1-dvm, nrow = ngrps, ncol = ngrid, byrow = TRUE)

    ### REMOVE DVM
    dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
    phyto_theta = matrix(1, nrow = ngrps, ncol = ngrid, byrow = TRUE)
    ###

    carn_grps = which(grp$type == 'C')
    phyto_theta[carn_grps,] = 0 # Carnivorous groups can't eat phyto

    cc_phyto = 0.1   # Carbon content of phytoplankton size classes

    ## Makes the model object, full of constant functions for model
    model = list(
      param = param,
      environ = environ,
      ngrid = ngrid,
      ngridPP = ngridPP,

      # Phytoplankton abundance
      nPP = 10^(environ$a)*(w_phyto^(environ$b)), # phyto abundance spectrum

      # Grid parameters
      w = w,
      dx = dx,
      w_phyto = w_phyto,

      # Group parameters storage
      phyto_growthkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # predation on phytoplankton
      phyto_diffkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # diffusion from phytoplankton consumption
      phyto_dietkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # diet from phytoplankton
      dynam_growthkernel =  array(NA, dim = c(ngrps, ngrid, ngrid)), # predation on zoo and fish
      dynam_diffkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # diffusion from zoo and fish consumption
      dynam_dietkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # diet from zoo and fish
      dynam_mortkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # mortality from predation on dynamic component
      M_sb = matrix(0, nrow = ngrps, ncol = ngrid), # senescence mortality
      fish_mort = matrix(0, nrow = ngrps, ncol = ngrid), # fishing mortality

      # Output storage
      N = array(0, dim = c(nsave, ngrps, ngrid)), # dynamic abundance spectrum
      Z = array(0, dim = c(nsave, ngrps, ngrid)), # total mortality
      gg = array(0, dim = c(nsave, ngrps, ngrid)), # growth
      diet = array(0, dim = c(nsave, c(ngrps), c(ngrps+3))), # diet
      Biomass = matrix(0, nrow = nsave, ncol = ngrps), # biomass of each group
      Diff = array(0, dim = c(nsave, ngrps, ngrid)) # save diffusion
    )

    # GGE for different groups
    assim_phyto =  (param$groups$alpha)*cc_phyto # Phytoplankton
    assim_dynam =  matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE) # rows are predators, columns are prey

    #### INITIAL DYNAMIC POPULATION ABUNDANCES
    a_dynam = 10^(environ$a)*(w[1]^(environ$b+1)) # calculate coefficient for initial dynamic spectrum, so that N(w_phyto) equals
    # N(w_dynam) at w[1]

    # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
    tempN = matrix(a_dynam*(w)^-1, nrow = ngrps, ncol = ngrid, byrow = TRUE)
    props_z = grp$prop[zoo_grps] # Zooplankton proportions
    tempN[zoo_grps,] = props_z*tempN[zoo_grps,] # Set abundances of diff zoo groups based on smallest size class proportions
    tempN[fish_grps,] = (1/num_fish)*tempN[fish_grps,] # Set abundandances of fish groups based on smallest size class proportions

    # For each group, set densities at w > Winf and w < Wmin to 0
    tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Winf) Winf < wx, Winf = (grp$Wmax)))] = 0
    tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Wmin) Wmin > wx, Wmin = (grp$W0)))] = 0
    model$N[1,,] = tempN

    # Fishing mortality
    model$fish_mort[fish_grps, c(w >= 1)] = param$f_mort

    ### MATRICES FOR LOG TRANSFORM OF EQUATION
    # Predators are rows, phyto prey weights are columns
    gg_log_t_phyto = ((w^-1) %*% t(w_phyto))/log(10) # Growth
    diff_log_t_phyto = ((w^-2) %*% t(w_phyto^2))/log(10) # Diffusion
    diet_log_t_phyto = matrix(w_phyto, nrow = length(w), ncol = length(w_phyto), byrow = TRUE) # Diet/Ingestion

    # Predators are rows, dynam prey weights are columns
    gg_log_t_dynam = ((w^-1) %*% t(w))/log(10) # Growth
    diff_log_t_dynam = ((w^-2) %*% t(w^2))/log(10) # Diffusion
    diet_log_t_dynam = matrix(w, nrow = length(w), ncol = length(w), byrow = TRUE) # Diet/ingestion

    ### PREDATION KERNELS FOR PHYTOPLANKTON SPECTRUM AND DYNAMIC SPECTRUM
    phyto_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngridPP)
    dynam_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid)
    phyto_prey_weight_matrix = matrix(w_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)
    dynam_prey_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid, byrow = TRUE)

    ## Search Volume storage
    SearchVol = matrix(NA, nrow = ngrps, ncol = ngrid) # Search volume

    # Simpson's Rule matrices for growth, diffusion and mortality integrals
    simp_phyto = array(1, dim = ngridPP)
    simp_phyto[c(seq(2,ngridPP-1,2))] = 4
    simp_phyto[c(seq(3,ngridPP-1,2))] = 2
    sm_phyto = matrix(simp_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)*(dx/3)

    simp_dynam = array(1, dim = ngrid)
    simp_dynam[c(seq(2,ngrid-1,2))] = 4
    simp_dynam[c(seq(3,ngrid-1,2))] = 2
    sm_dynam = matrix(simp_dynam, nrow = ngrid, ncol = ngrid, byrow = TRUE)*(dx/3)

    ## Temperature Effect Matrix
    # Effect of temperature on feeding and predation rate
    #temp_flag <- 2.4^((environ$sst)/10)
    #temp_cil <-  2.8^((environ$sst)/10)
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
    temp_zoo <- rep(2.^((environ$sst - 30)/10), num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
    temp_fish <- rep(2.^((environ$sst - 30)/10), num_fish)
    temp_effect <- matrix(c(temp_zoo, temp_fish), nrow = ngrps, ncol = ngrid)

    #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
    for(i in 1:ngrps){
      ## Senescence mortality
      if(i < 10){
        model$M_sb[i,] = param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
        model$M_sb[i, 10^(grp$Wmax[i]) < w] = 0
        model$M_sb[i, 10^(grp$Wmat[i]) > w] = 0
      }

      if(i > 9){
        model$M_sb[i,] = 0.1*param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
        model$M_sb[i, 10^(grp$Wmax[i]) < w] = 0
        model$M_sb[i, 10^(grp$Wmat[i]) > w] = 0
      }

      ### Search volume
      SearchVol[i,] = (grp$gamma[i])*(w^(grp$q[i]))
      SearchVol[i, 10^(grp$Wmax[i]) < w] = 0
      SearchVol[i, 10^(grp$W0[i]) > w] = 0

      ### Predation Kernels
      if(is.na(grp$m[i]) == FALSE){ # If group has an m-value (zooplankton)
        # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
        D.z = 2*(3*w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
        betas =  (exp(0.02*log(D.z)^2 - grp$m[i] + 1.832))^3 # Wirtz's equation
        beta_mat_phyto = matrix(betas, nrow = ngrid, ncol = ngridPP)
        beta_mat_dynam = matrix(betas, nrow = ngrid, ncol = ngrid)

        # Calculate feeding kernels
        sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                              phyto_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
        sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                              dynam_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)

        if(i == 3){
          sp_phyto_predkernel <- matrix(sp_phyto_predkernel[44,], nrow = ngrid, ncol = ngridPP, byrow = TRUE)
          sp_dynam_predkernel <- matrix(sp_dynam_predkernel[44,], nrow = ngrid, ncol = ngrid, byrow = TRUE)
        }
        if(i == 8){
          sp_phyto_predkernel <- matrix(sp_phyto_predkernel[61,], nrow = ngrid, ncol = ngridPP, byrow = TRUE)
          sp_dynam_predkernel <- matrix(sp_dynam_predkernel[61,], nrow = ngrid, ncol = ngrid, byrow = TRUE)
        }

      } else { # If group does not have an m-value (fish)
        beta_mat_phyto = matrix(grp$beta[i], nrow = ngrid, ncol = ngridPP)
        beta_mat_dynam = matrix(grp$beta[i], nrow = ngrid, ncol = ngrid)

        # Calculate feeding kernels
        sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                              phyto_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
        sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                              dynam_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
      }

      ### GROWTH INTEGRAL CONSTANTS
      # Predators are rows, prey are columns
      model$phyto_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*gg_log_t_phyto*sm_phyto
      model$dynam_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*gg_log_t_dynam*sm_dynam

      ### DIET INTEGRAL CONSTANTS
      # Predators are rows, prey are columns
      model$phyto_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*diet_log_t_phyto*sm_phyto
      model$dynam_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*diet_log_t_dynam*sm_dynam

      ### DIFFUSION INTEGRAL CONSTANTS
      # Predators are rows, prey are columns
      model$phyto_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*diff_log_t_phyto*sm_phyto
      model$dynam_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*diff_log_t_dynam*sm_dynam

      ### MORTALITY INTEGRAL CONSTANTS
      # Prey are rows, predators are columns
      model$dynam_mortkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid, byrow = TRUE)*
        t(sp_dynam_predkernel)*sm_dynam
    }

    no_sen = which(grp$species == c("Flagellates", "Ciliates")) # no senescence mortality for flagellates and ciliates
    #model$M_sb[c(ngrps)] = 0
    model$M_sb = temp_effect*model$M_sb # Incorporate temp effect on senscence mortality

    ## Incorporate dvm, temperature effects and gross growth efficiency (assim)
    model$phyto_growthkernel = sweep(sweep(model$phyto_growthkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto, "*")
    model$phyto_diffkernel = sweep(sweep(model$phyto_diffkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto^2, "*")
    model$phyto_dietkernel =  sweep(sweep(model$phyto_dietkernel, c(1,2), phyto_theta, "*"), 1, 1, "*")

    # Dim 1 = pred group, dim2 = pred sizes, dim 3 = prey group, dim 4 = prey sizes
    model$dynam_growthkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_growthkernel, "*"),
                                           c(1,3), assim_dynam, "*"), c(1,2), temp_effect, "*")
    model$dynam_diffkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_diffkernel, "*"), c(1,3), assim_dynam^2, "*"),
                                   c(1,2), temp_effect^2, "*")
    model$dynam_dietkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_dietkernel, "*"),
                                         c(1,3), 1, "*"), c(1,2), temp_effect, "*")

    # Dim 1 = prey group, dime 2 = prey sizes, dim 3 = pred group, dim 4 = pred sizes
    model$dynam_mortkernel = sweep(aperm(sweep(aperm(dynam_theta, c(3,1,4,2)),
                                               c(2,3,4), model$dynam_mortkernel, "*"), c(1,3,2,4)),
                                   c(3,4), temp_effect, "*")

    #### Because phyto spectrum is constant, we can solve the phyto component of growth, and diffusion before time loop
    model$ingested_phyto = temp_effect*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2)) # Ingested phyto
    model$diff_phyto = temp_effect^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2)) # Diffusion from phyto
    model$diet_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP, "*"), dims = 2)) # Diet of total phyto

    ## Diet of phyto from pico, nano and micro size classes
    model$diet_pico_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) < -11.5), "*"), dims = 2))
    model$diet_nano_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -11.5 & log10(model$w_phyto) < -8.5), "*"), dims = 2))
    model$diet_micro_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -8.5), "*"), dims = 2))

    return(model)
  } # End of Setup function


  ## Run model forward in time
  Project <- function(model, fish_on){

    # Pull out some useful parameters - just a shortcut
    param <- model$param
    grp <- param$group
    ngrid <- model$ngrid
    ngridPP <- model$ngridPP
    ngrps <- param$ngrps
    dt <- param$dt
    fish_grps <- param$fish_grps
    zoo_grps <- param$zoo_grps
    dx <- model$dx
    w <- model$w
    w_phyto <- model$w_phyto
    dw <- model$dw
    assims <- param$nutrition/0.1*(grp$alpha)
    w0idx <- which(grp$W0 > min(grp$W0) & is.na(grp$prop) == FALSE)
    #w0idx <- which(is.na(grp$prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    wMaxidx <- which(grp$Wmax < max(grp$Wmax) & is.na(grp$prop) == FALSE)
    w0maxs <- rep(0, length(w0idx))
    props_z <- grp$prop[w0idx] # Zooplankton proportions

    for(i in 1:length(w0idx)){ # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(w), digits = 2) == grp$W0[w0idx[i]])
      w0maxs[i] <- which(round(log10(w), digits = 2) == grp$Wmax[wMaxidx[i]])
    }

    # Handy stuff
    idx <- 2:ngrid # size class sequence
    itimemax  <- param$tmax / dt  #max index of time array

    # Matrices for MvF and MvF-D numeric solution
    A.iter <- C.iter <- S.iter <- A <- B <- C <- S <- matrix(0,nrow=ngrps,ncol=ngrid)

    # Temporary Matrices that get updated each time step
    # some of these saved for output
    N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes
    nPP <- model$nPP # Abundances of phytoplankton spectrum

    pb = txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3) # Initial progress bar

    # BIG TIME LOOP
    for (itime in 1:itimemax){

      setTxtProgressBar(pb, itime) # Update progress bar

      ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps
      # dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
      N_array <- aperm(replicate(ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(ngrps, N_array), c(4,1,2,3))


      ### GROWTH
      gg <- (model$ingested_phyto +
               rowSums(rowSums(model$dynam_growthkernel*N_array, dims = 3), dims = 2))

      ### MORTALITY
      # Predation mortality
      M2 <- (rowSums(rowSums(model$dynam_mortkernel*N_array, dims = 3), dims = 2))

      # Total dynamic spectrum mortality
      Z = M2 + model$M_sb  + model$fish_mort


      ### DIFFUSION
      #diff <- (model$diff_phyto + rowSums(rowSums(model$dynam_diffkernel*N_array, dims = 3), dims = 2))

      ### MvF WITH DIFFUSION ALGORITHM

      idx.iter <- 2:ngrid
      idx <- 2:(ngrid-1)

      # Numerical implementation matrices (for MvF without diffusion)
      A.iter[,idx.iter] <- dt/dx*gg[,idx.iter-1]
      #A.iter[, idx.iter] <- sweep(gg[, idx.iter - 1, drop = FALSE] * dt, 2, dx, "/")
      C.iter[,idx.iter] <- 1 + dt*Z[,idx.iter] + dt/dx*gg[,idx.iter]
      #C.iter[, idx.iter] <- 1 + sweep(gg[, idx.iter, drop = FALSE] * dt, 2, dx, "/") + Z[, idx.iter, drop = FALSE] * dt
      S.iter[,idx.iter] <- N[,idx.iter]
      N.iter <- N

      curr_min_size <- vector()
      curr_max_size <- vector()
      #C.iter[,w0idx] <- 1 + gg[w0idx,w0mins] * dt / dx[w0idx] + gg[w0idx,w0mins] * dt
      for(i in 1:ngrps){
        ## Set size range index for current group
        curr_min_size[i] = which(round(log10(w), digits = 2) == param$groups$W0[i])
        curr_max_size[i] = which(round(log10(w), digits = 2) == param$groups$Wmax[i])
        # idx_curr = (curr_min_size+1):curr_max_size
      }

         N2 <- inner_project_loop(no_sp = ngrps, no_w = ngrid, n = N.iter,
                                  A = A.iter, B = C.iter, S = S.iter,
                                  w_min_idx = curr_min_size, w_max_idx = curr_max_size)
      # #


      # for(j in w0mins[i]+1:log10(param$wMax)){## Find the abundance at the next size class with standard MvF
      #     N[i,j] <- (S.iter[i,j] + A.iter[i,j]*N[i,j-1])/(C.iter[i,j])
      #
      #  MvF without diffusion for last size class
          # for (i in length(w0maxs)) {
          #  N[i,w0maxs[i]] = 0
          #  N[i,curr_min_size] <- N[i,curr_min_size] # Keep starting sizes constant
          # }
      # # #
      #    }
      # }

      #### Keep smallest fish community size class as equal to equivalent zooplankton size class

      ### Keep smallest zooplankton size class abundnace
      ### for each group locked to others in size spectrum
      for(i in 0:length(w0idx)){
        w_min_curr = w0mins[i]
        exclude_mins = w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] = props_z[i]*sum(N[-exclude_mins, w_min_curr])
      }


      fish_mins = unlist(lapply(param$groups$W0[fish_grps],
                                function(x){which(round(log10(model$w), digits = 2) == x)}))

      if(length(fish_grps) > 1){
        N[fish_grps,fish_mins] = (1/3)*(colSums(N[-fish_grps,fish_mins]))
      }else{
        N[fish_grps, fish_mins] = sum(N[-fish_grps, fish_mins])
      }

      if(fish_on == FALSE){
        N[fish_grps,] = 0 # switch off fish_groups
      }

      # Save results:
      if((itime %% param$isave) == 0){
        isav=itime/param$isave

        ## Phytoplankton diet
        pico_phyto_diet = rowSums(model$diet_pico_phyto*N) # Pico-phytoplankton
        nano_phyto_diet = rowSums(model$diet_nano_phyto*N) # Nano-phytoplankton
        micro_phyto_diet = rowSums(model$diet_micro_phyto*N) # Micro-phytoplankton

        phyto_diet = cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)

        ## Functional group diet
        dynam_diet =  rowSums(aperm(rowSums(sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

        model$diet[isav,,1:3] = phyto_diet
        model$diet[isav,,c(4:15)] = dynam_diet

        model$N[isav,,] <- N # Save abundance

        ## NEED TO CHECK WHAT THIS IF STATEMENT IS FOR, IT DOESN'T APPEAR TO DO ANYTHING
        if(length(zoo_grps) > 1){
          model$N[isav,c(1,2),c(60,61)] <- 0
        }

        ## Save biomass
        #  model$Biomass[isav,] <- rowSums(model$N[isav,,] # Save biomass
        #                                  *matrix(model$w, nrow = ngrps, ncol = ngrid, byrow = TRUE))

        ## Save mortality rates
        #	model$Z[isav,,] <- M2 # Save total predation mortality rates

        ## Save growth
        model$gg[isav,,] <-  (model$diet_phyto +
                                rowSums(rowSums(model$dynam_dietkernel*N_array, dims = 3), dims = 2))
      }

    } # End of time loop

    return(model)

  }


  # enviro_vector <- data.frame("sst" = sst, "chlo" = chlo, "a" = a, "b" = b,
  #                             "phyto_max" = phyto_max, "dt" = dt)

  ##################### RUN THE MODEL ################################################
  param <- params(Groups, enviro_vector, tmax = tmaxx, f_mort = 0) # Set up parameter list
  model <- Setup(param) # Set up model equation stuff
  modelss <- Project(model, fish_on = TRUE) # Run the model

  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances = colMeans(modelss$N[(ceiling(0.5*dim(modelss$N)[1])):(dim(modelss$N)[1]),,], dim = 1)
  ave_diets = colMeans(modelss$diet[(ceiling(0.5*dim(modelss$diet)[1])):(dim(modelss$diet)[1]),,], dim = 1)
  ave_growth = colMeans(modelss$gg[(ceiling(0.5*dim(modelss$gg)[1])):(dim(modelss$gg)[1]),,], dim = 1)
  ave_pred = colMeans(modelss$Z[(ceiling(0.5*dim(modelss$Z)[1])):(dim(modelss$Z)[1]),,], dim = 1)
  results = list("abundances" = ave_abundances, "diets" = ave_diets, "growth" = ave_growth) #, "diets" = ave_diets
  return(results)
}

################ run it ####################


Groups <- read.csv("~/../Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/TestGroups.csv")
enviro_vector <- data.frame("sst" = 21.6957, "chlo" = 0.3763819, "a" = -1.27508, "b" = -0.9491563, "phyto_max" = -7.1, "dt" = .1)
tmaxx <- 50
#res <- multi_zoo_slurm_filter(21.7,0.709,-.234, -.869,-7.7,0.1,200)
#res2 <- multi_zoo_slurm_filter(enviro_data[1,,drop=FALSE])
start_time <- Sys.time()
res <- multi_zoo_slurm_filter(enviro_vector)
end_time <- Sys.time()
two <- end_time - start_time
two
two*(1000/tmaxx*enviro_vector$dt/0.1)
w <- seq(from=-10.7,to=7,0.1)
x11()
plot(w,log10(colSums(abs(res$abundances))))
View(res$abundances)
##
resp <- apply(enviro_data, 1, multi_zoo_slurm_filter)#, mc.cores = 1, SIMPLIFY = FALSE)



res <- list()
enviro_data$dt <- 0.01
tmaxx <- 10
for (i in 1:dim(enviro_data)[1]){ # %do% #dim(enviro_data)[1]) %do%
  .start_time <- Sys.time()
  enviro_v <- as.data.frame(enviro_data[i,,drop=FALSE])
  res[[i]] <- multi_zoo_slurm_filter(enviro_v)
  .end_time <- Sys.time()
  print(.end_time - .start_time)
}
