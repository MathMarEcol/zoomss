#' Run ZooMSS Model Forward in Time
#'
#' @title Execute the main ZooMSS simulation loop with dynamic environmental forcing
#' @description Runs the ZooMSS model forward in time, updating environmental conditions
#'   and population dynamics at each time step using the McKendrick-von Foerster framework.
#' @details This is the core simulation engine of ZooMSS that:
#'
#'   **Environmental Dynamics:**
#'   - Updates phytoplankton abundance spectrum based on chlorophyll time series
#'   - Applies temperature effects on zooplankton and fish metabolism
#'   - Recalculates feeding kernels with current environmental conditions
#'
#'   **Population Dynamics:**
#'   - Solves McKendrick-von Foerster equation for size-structured growth
#'   - Updates feeding interactions between all size classes and groups
#'   - Calculates mortality from predation, senescence, and fishing
#'   - Handles recruitment and boundary conditions for each functional group
#'
#'   **Time Integration:**
#'   - Processes model through all time steps with adaptive environmental forcing
#'   - Saves output at specified intervals for memory efficiency
#'   - Maintains mass balance and numerical stability throughout simulation
#'
#'   Unlike static models, this version dynamically updates phytoplankton spectra
#'   and temperature effects at each time step based on provided environmental data.
#'
#' @param model Model object created by zoomss_setup containing:
#'   - param: Complete parameter list with environmental time series
#'   - Feeding kernels and biological rate parameters
#'   - Initial conditions and model structure
#'
#' @return List containing complete model output:
#'   \itemize{
#'     \item param: Model parameters used in simulation
#'     \item N: Abundance time series (time x groups x size classes)
#'     \item gg: Growth rate time series
#'     \item diet: Diet composition time series
#'     \item Z: Mortality rate time series
#'     \item time: Time values corresponding to saved results (accounting for isave)
#'     \item w: Size class weights (g)
#'     \item Additional time series data and model results
#'   }
#'
#' @examples
#' \dontrun{
#' # Set up model parameters and structure
#' params <- zoomss_params(Groups, input_params)
#' model <- zoomss_setup(params)
#'
#' # Run the simulation
#' results <- zoomss_run(model)
#'
#' # Access final abundances
#' final_abundances <- results$N[dim(results$N)[1],,]
#' }
#'
#' @noRd
#'
zoomss_run <- function(model){

  # Pull out some useful parameters - just a shortcut
  param <- model$param
  dt <- model$param$dt[1]  # Use scalar value instead of vector
  dx <- model$param$dx     # This should already be scalar
  ngrps <- model$param$ngrps
  ngrid <- model$param$ngrid
  w <- model$param$w
  W0 <- param$Groups$W0
  Wmax <- param$Groups$Wmax
  fish_grps <- model$param$fish_grps
  assim_eff <- model$assim_eff
  dynam_dietkernel <- model$dynam_dietkernel
  dynam_growthkernel <- model$dynam_growthkernel
  dynam_mortkernel <- model$dynam_mortkernel
  dynam_diffkernel <- model$dynam_diffkernel

  # Get pre-calculated phytoplankton and temperature time series from param
  phyto_int <- param$phyto_int
  phyto_slope <- param$phyto_slope
  temp_eff_zoo <- param$temp_eff_zoo
  temp_eff_fish <- param$temp_eff_fish

  curr_min_size <- vector()
  curr_max_size <- vector()
  for (i in 1:ngrps){
    curr_min_size[i] <- which(round(log10(w), digits = 2) == W0[i])
    curr_max_size[i] <- which(round(log10(w), digits = 2) == Wmax[i])
  }

  idx_iter <- 2:ngrid
  idx <- 2:(ngrid-1)
  itimemax  <- param$itimemax  # Number of time points to simulate
  diff_off <- identical(param$diffusion, "off") # diffusion off: skip its calculations

  if(length(param$zoo_grps) > 1){ # If there's only one zoo group, then you do not need w0idx. All this stuff gives you info about all zoo groups except the smallest zoo group.
    w0idx <- which(W0 > min(W0) & is.na(param$Groups$Prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    props_z <- param$Groups$Prop[w0idx] # Zooplankton proportions

    for(i in 1:length(w0idx)){
      # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(w), digits = 2) == W0[w0idx[i]])
    }
  }

  # Matrices for MvF and MvF-D numeric solution
  A_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  C_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  S_iter <- matrix(0, nrow = ngrps, ncol = ngrid)

  A <- matrix(0, nrow = ngrps, ncol = ngrid)
  B <- matrix(0, nrow = ngrps, ncol = ngrid)
  C <- matrix(0, nrow = ngrps, ncol = ngrid)
  S <- matrix(0, nrow = ngrps, ncol = ngrid)

  # Temporary Matrices that get updated each time step some of these saved for output
  N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes

  pb <- progress::progress_bar$new(
    format = "ZooMSS Time Loop [:bar] :percent eta: :eta elapsed: :elapsed",
    total = itimemax,
    width = 60,
    show_after = 0
  )

  # BIG TIME LOOP
  for (itime in 1:itimemax){

    pb$tick() # Update progress bar

    # DYNAMIC ENVIRONMENTAL FORCING - Update for current time step
    current_phyto_int <- phyto_int[itime]      # Phytoplankton intercept (varies with chlorophyll)
    current_phyto_slope <- phyto_slope[itime]  # Phytoplankton slope (varies with chlorophyll)

    # Update phytoplankton abundance spectrum with current environmental parameters
    model$nPP <- 10^(current_phyto_int)*(param$w_phyto^(current_phyto_slope))

    # Update temperature effects matrix for all groups based on current SST
    # Temperature affects feeding rates, growth, and mortality
    if(param$num_zoo > 0) {
      for(i in 1:param$num_zoo) {
        zoo_grp_idx <- param$zoo_grps[i]
        model$temp_eff[zoo_grp_idx, ] <- temp_eff_zoo[itime, i]
      }
    }

    if(param$num_fish > 0) {
      for(i in 1:param$num_fish) {
        fish_grp_idx <- param$fish_grps[i]
        model$temp_eff[fish_grp_idx, ] <- temp_eff_fish[itime, i]
      }
    }

    # Calculate phytoplankton feeding using CURRENT environmental conditions
    # These must be recalculated each time step because both nPP and temp_eff change
    current_ingested_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2))
    current_diff_phyto <- if (diff_off) NULL else model$temp_eff^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2))

    # Update senescence mortality with current temperature effects (no compounding - always from base)
    model$M_sb <- model$M_sb_base * model$temp_eff

    # Calculate multipliers for growth, predation, and diffusion
    growth_multiplier <- colSums(N * assim_eff) # 1 x n_sizes
    predation_multiplier <- N * model$temp_eff # n_species x n_sizes (now using updated temp_eff)
    diffusion_multiplier <- if (diff_off) NULL else colSums(N * (assim_eff^2)) # 1 x n_sizes

    ### DO GROWTH
    # Apply temperature effects to growth kernel (consistent with model temperature scaling)
    temp_growth_kernel <- sweep(dynam_growthkernel, c(1,2), model$temp_eff, '*')
    dim(temp_growth_kernel) <- c(ngrps*ngrid, ngrid)
    cs <- .colSums(growth_multiplier * t(temp_growth_kernel), m = ngrid, n = ngrps*ngrid)
    dim(cs) <- c(ngrps, ngrid)

    gg <- current_ingested_phyto + cs

    ### DO MORTALITY

    sw2 <- sweep(dynam_mortkernel, c(2,3), predation_multiplier, '*') # n_sizes x n_species x n_sizes
    ap2 <- aperm(sw2, c(2,3,1))
    M2 <- .colSums(colSums(ap2),ngrid,ngrid) # 1 x n_sizes
    Z <- sweep(model$M_sb + model$fish_mort, 2, M2, '+') # Total dynamic spectrum mortality (n_species x n_sizes)
    rm(sw2, ap2)


    ### DO DIFFUSION
    # Apply temperature effects to diffusion kernel (consistent with model temperature scaling)
    if (diff_off) {
      diff <- matrix(0, nrow = ngrps, ncol = ngrid)
    } else {
      temp_diff_kernel <- sweep(dynam_diffkernel, c(1,2), model$temp_eff^2, '*')
      dim(temp_diff_kernel) <- c(ngrps*ngrid, ngrid)
      cs <- .colSums(diffusion_multiplier * t(temp_diff_kernel), m = ngrid, n = ngrps*ngrid)
      dim(cs) <- c(ngrps, ngrid)
      diff <- current_diff_phyto + cs
    }

    # Express the diffusion term in log10 weight (see zoomss_diffusion_terms)
    dterms <- zoomss_diffusion_terms(diff, model$temp_eff, param$diffusion)
    diff <- dterms$diff

    ### MvF WITH DIFFUSION ALGORITHM
    # Numerical implementation matrices (for MvF without diffusion)
    A_iter[,idx_iter] <- dt/dx * gg[,idx_iter-1] # Growth stuff
    C_iter[,idx_iter] <- 1 + dt * Z[,idx_iter] + dt/dx * gg[,idx_iter] # Mortality
    S_iter[,idx_iter] <- N[,idx_iter] # N at.....
    N_iter <- N # Current Abundance
    N_iter[1,1] <- N[1,1] # This forces R to make a copy of the variable. Otherwise N is linked to N_iter in the Rcpp code and they change together.

    # Numerical implementation matrices (for MvF WITH diffusion)
    A[,idx] <- dt/dx * (gg[,idx-1] + diff[,idx-1] * (dterms$drift+1/(2*dx))) # Growth stuff
    B[,idx] <- diff[,idx+1] * dt/(2*dx^2) # Diffusion term
    C[,idx] <- 1 + dt * Z[,idx] + dt/dx*(gg[,idx] + diff[,idx] * (dterms$drift+1/dx)) # Mortality
    S[,idx] <- N[,idx]

    # The original Base R code for the MvF equation
    N <- zoomss_mvf(ngrps, curr_min_size, curr_max_size,
                           A_iter, C_iter, N_iter, S_iter,
                            A, B, C, N, S)

    #### Keep smallest fish community size class as equal to equivalent zooplankton size class
    ### Keep smallest zooplankton size class abundnace for each group locked to others in size spectrum
    if(length(param$zoo_grps) > 1){ # If you only have one zoo group, it will be locked to phyto spectrum so you do not need to do this
      for(i in 1:length(w0idx)){
        w_min_curr <- w0mins[i]
        exclude_mins <- w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] <- props_z[i] * sum(N[-exclude_mins, w_min_curr])
      }
    }

    fish_mins <- unlist(lapply(W0[fish_grps],
                               function(x){which(round(log10(w), digits = 2) == x)}))

    if(length(fish_grps) > 1 & length(param$zoo_grps) > 1){
      N[fish_grps,fish_mins] <- (1/length(fish_grps))*(colSums(N[-fish_grps,fish_mins]))
    }else{
      N[fish_grps, fish_mins] <- (1/length(fish_grps))*sum(N[-fish_grps, fish_mins])
    }


    # Save results (at regular intervals AND always at the final time step):
    # Save output at regular intervals only
    save_this_step <- (itime %% param$isave) == 0

    if(save_this_step){
      isav <- itime/param$isave

      ## Phytoplankton diet - calculate using CURRENT dynamic conditions
      # Use current phytoplankton spectrum and temperature effects (both change each time step)
      current_pico_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) < -11.5), "*"), dims = 2))
      current_nano_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) >= -11.5 & log10(model$param$w_phyto) < -8.5), "*"), dims = 2))
      current_micro_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) >= -8.5), "*"), dims = 2))

      # Calculate total consumption by multiplying diet potential by current abundances
      pico_phyto_diet <- rowSums(current_pico_diet*N)   # Pico-phytoplankton consumption
      nano_phyto_diet <- rowSums(current_nano_diet*N)   # Nano-phytoplankton consumption
      micro_phyto_diet <- rowSums(current_micro_diet*N) # Micro-phytoplankton consumption

      ## Functional group diet - use proper kernel-based calculation (consistent with model structure)
      ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps: dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
      N_array_temp <- aperm(replicate(ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(ngrps, N_array_temp), c(4,1,2,3))
      # Apply temperature effects to diet kernel like phytoplankton diets
      temp_dynam_dietkernel <- sweep(dynam_dietkernel, c(1,2), model$temp_eff, "*")
      dynam_diet <- rowSums(aperm(rowSums(sweep(temp_dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

      # Save current time accounting for isave interval
      model$time[isav] <- param$time[itime]

      model$N[isav,,] <- N # Save N by taxa and size
      model$Z[isav,,] <-  Z ## Save mortality
      model$gg[isav,,] <-  gg ## Save growth

      model$diet[isav,,1:3] <- cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)
      model$diet[isav,,c(4:(dim(param$Groups)[1]+3))] <- dynam_diet

    }
  } # End of time loop
  return(model)

}

#' Diffusion term in log10 weight
#'
#' @title Convert the diffusion term to the log10 weight grid
#' @description The model equation is written in weight $w$, with $N$ the abundance per unit
#'   weight, $G$ the growth rate, $Z$ the mortality rate and $D(w)$ the diffusion coefficient
#'   (which carries the temperature effect $\tau$ as in the published model, squared):
#'   \deqn{\partial N/\partial t = -\partial(GN)/\partial w + \tfrac12 \partial^2(DN)/\partial w^2 - ZN.}
#'   The model solves it on $x = \log_{10} w$ for the abundance per unit $x$, $n = N w \log(10)$.
#'   With $g = G/(\log(10) w)$, the solver's equation is
#'   \deqn{\partial n/\partial t = -\partial(gn)/\partial x + \tfrac12\partial^2(\mathtt{diff}\, n)/\partial x^2
#'     - \mathtt{drift}\,\partial(\mathtt{diff}\, n)/\partial x - Zn,}
#'   where `diff` and `drift` are returned by this function. The exact transform has
#'   $\mathtt{diff} = D/(\log(10)^2 w^2)$ and $\mathtt{drift} = -\log(10)/2$, a drift towards
#'   smaller sizes, because equal steps in $x$ are unequal steps in $w$.
#' @details Write $D_0$ for the diffusion coefficient without the temperature effect, so that
#'   $D = \tau^2 D_0$ in the published model. `diffusion` is one of
#'   \describe{
#'     \item{"published"}{As in the published code: `diff` is used as built in zoomss_setup and
#'       zoomss_run, $\mathtt{diff} = \tau^2 D_0/(\log(10) w^2)$ (one division by $\log(10)$ where the
#'       transform needs two), and $\mathtt{drift} = +\log(10)/2$ (the wrong sign).}
#'     \item{"on"}{The exact transform of the published equation: `diff` is divided by
#'       $\log(10)$, so that $\mathtt{diff} = \tau^2 D_0/(\log(10)^2 w^2)$, and
#'       $\mathtt{drift} = -\log(10)/2$. Temperature is squared as in the published model.}
#'     \item{"tau_once"}{As "on", but with the temperature effect entering diffusion once:
#'       `diff` is divided by $\log(10)\tau$, so that $\mathtt{diff} = \tau D_0/(\log(10)^2 w^2)$.
#'       Diffusion is the feeding rate times the square of the weight gained per meal;
#'       temperature scales the feeding rate, and each meal still adds the same weight.}
#'     \item{"off"}{Diffusion off.}
#'   }
#'   With $\mathtt{drift} = -\log(10)/2$ the coefficient `A` in zoomss_run, which multiplies
#'   the abundance in the size class below, is non-negative only while
#'   $dx < 1/\log(10) \approx 0.43$. This holds at the default $dx = 0.1$.
#' @param diff Matrix (groups x size classes) of diffusion coefficients as built in zoomss_run,
#'   which include the temperature effect squared.
#' @param temp_eff Matrix (groups x size classes) of temperature effects.
#' @param diffusion One of "on", "tau_once", "published" or "off".
#' @return List with `diff`, the diffusion coefficients to use in the solver, and `drift`, the
#'   coefficient of the first-derivative term.
#' @noRd
zoomss_diffusion_terms <- function(diff, temp_eff, diffusion = "on"){
  switch(diffusion,
         published = list(diff = diff, drift = log(10)/2),
         on = list(diff = diff / log(10), drift = -log(10)/2),
         tau_once = list(diff = diff / (log(10) * temp_eff), drift = -log(10)/2),
         off = list(diff = diff * 0, drift = 0),
         stop("diffusion must be one of \"on\", \"tau_once\", \"published\" or \"off\""))
}
