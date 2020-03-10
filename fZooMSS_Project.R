## Run model forward in time
Project <- function(model, fish_on){

  # Pull out some useful parameters - just a shortcut
  param <- model$param
  grp <- param$group
  ngrid <- model$ngrid
  ngrps <- param$ngrps
  dt <- param$dt
  fish_grps <- param$fish_grps
  zoo_grps <- param$zoo_grps
  dx <- model$dx
  w <- model$w

  idx.iter <- 2:ngrid
  idx <- 2:(ngrid-1)
  itimemax  <- param$tmax / dt  #max index of time array

  if(length(zoo_grps) > 1){ # If there's only one zoo group, then you do not need w0idx. All this stuff gives you info about all zoo groups except the smallest zoo group.
    w0idx <- which(grp$W0 > min(grp$W0) & is.na(grp$prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    props_z <- grp$prop[w0idx] # Zooplankton proportions

    for(i in 1:length(w0idx)){ # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(w), digits = 2) == grp$W0[w0idx[i]])
    }
  }


  # Matrices for MvF and MvF-D numeric solution
  A.iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  C.iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  S.iter <- matrix(0, nrow = ngrps, ncol = ngrid)

  A <- matrix(0, nrow = ngrps, ncol = ngrid)
  B <- matrix(0, nrow = ngrps, ncol = ngrid)
  C <- matrix(0, nrow = ngrps, ncol = ngrid)
  S <- matrix(0, nrow = ngrps, ncol = ngrid) # Previous abbundance

  # Temporary Matrices that get updated each time step some of these saved for output
  N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes

  pb <- txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3) # Initial progress bar

  # N_array <- array(0,c(12,178,12,178)) # preallocate
  idx_array <- sort(array(matrix(1:2136),4562496))

  # BIG TIME LOOP
  for (itime in 1:itimemax){

    setTxtProgressBar(pb, itime) # Update progress bar
    # browser()

    ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps
    # dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes

    # temp_array <- array(N,4562496) # Convert N to vector and replicate at the same time
    N_array <- array(array(N,4562496)[idx_array], c(12, 178, 12, 178)) # Then rearrange the vector to the order it should be

    gg <- rowSums(rowSums(model$dynam_growthkernel*N_array, dims = 3), dims = 2) ### GROWTH
    M2 <- rowSums(rowSums(model$dynam_mortkernel*N_array, dims = 3), dims = 2) ### MORTALITY: Predation mortality
    diff <- rowSums(rowSums(model$dynam_diffkernel*N_array, dims = 3), dims = 2) ### DIFFUSION

    gg <- gg + model$ingested_phyto
    diff <- diff + model$diff_phyto
    Z <- M2 + model$M_sb  + model$fish_mort ### MORTALITY: Total dynamic spectrum mortality

    ### MvF WITH DIFFUSION ALGORITHM

    # Numerical implementation matrices (for MvF without diffusion)
    A.iter[,idx.iter] <- dt/dx * gg[,idx.iter-1] # Growth stuff
    C.iter[,idx.iter] <- 1 + dt * Z[,idx.iter] + dt/dx * gg[,idx.iter] # Mortality
    S.iter[,idx.iter] <- N[,idx.iter] # N at
    N.iter <- N # Current Abundance

    # Numerical implementation matrices (for MvF WITH diffusion)
    A[,idx] <- dt/dx * (gg[,idx-1] + diff[,idx-1] * (log(10)/2+1/(2*dx)))
    B[,idx] <- diff[,idx+1] * dt/(2*dx^2) # Diffusion term
    C[,idx] <- 1 + dt * Z[,idx] + dt/dx*(gg[,idx] + diff[,idx] * (log(10)/2+1/dx))
    S[,idx] <- N[,idx]

    for(i in 1:ngrps){
      ## Set size range index for current group
      curr_min_size = which(round(log10(w), digits = 2) == param$groups$W0[i])
      curr_max_size = which(round(log10(w), digits = 2) == param$groups$Wmax[i])
      idx_curr = (curr_min_size+1):curr_max_size

      for(j in idx_curr){## Find the abundance at the next size class with standard MvF
        N.iter[i,j] <- (S.iter[i,j] + A.iter[i,j]*N[i,j-1])/(C.iter[i,j])
        if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
          k <- j - 1
          N[i,k] <- (S[i,k] + A[i,k] * N[i,k-1] + B[i,k] * N.iter[i,k+1]) / C[i,k]
        }

        # MvF without diffusion for last size class
        if(j == idx_curr[length(idx_curr)]){
          N[i,j] <- 0
          # N[i,curr_min_size] <- N.iter[i,curr_min_size] # Keep starting sizes constant
        }
      }
    }

    #### Keep smallest fish community size class as equal to equivalent zooplankton size class
    ### Keep smallest zooplankton size class abundnace
    ### for each group locked to others in size spectrum

    if(length(zoo_grps) > 1){ # If you only have one zoo group, it will be locked to phyto spectrum so you do not need to do this
      for(i in 1:length(w0idx)){
        w_min_curr = w0mins[i]
        exclude_mins = w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] = props_z[i] * sum(N[-exclude_mins, w_min_curr])
      }
    }

    fish_mins <- unlist(lapply(param$groups$W0[fish_grps],
                              function(x){which(round(log10(model$w), digits = 2) == x)}))

    if(length(fish_grps) > 1 & length(zoo_grps) > 1){
      N[fish_grps,fish_mins] = (1/length(fish_grps))*(colSums(N[-fish_grps,fish_mins]))
    }else{
      N[fish_grps, fish_mins] = (1/length(fish_grps))*sum(N[-fish_grps, fish_mins])
    }

    if(fish_on == FALSE){
      N[fish_grps,] <- 0 # switch off fish_groups
    }

    # Save results:
    if((itime %% param$isave) == 0){
      isav <- itime/param$isave

      ## Phytoplankton diet
      pico_phyto_diet <- rowSums(model$diet_pico_phyto*N) # Pico-phytoplankton
      nano_phyto_diet <- rowSums(model$diet_nano_phyto*N) # Nano-phytoplankton
      micro_phyto_diet <- rowSums(model$diet_micro_phyto*N) # Micro-phytoplankton

      phyto_diet <- cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)

      ## Functional group diet
      dynam_diet <- 0 # dynam_diet <- rowSums(aperm(rowSums(sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

      model$diet[isav,,1:3] <- phyto_diet
      model$diet[isav,,c(4:(dim(grp)[1]+3))] <- dynam_diet

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
      model$gg[isav,,] <- (model$diet_phyto +
                             rowSums(rowSums(model$dynam_dietkernel*N_array, dims = 3), dims = 2))
    }

    # rm(N_array)
  } # End of time loop

  return(model)

}