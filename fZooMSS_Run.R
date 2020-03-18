## Run model forward in time using previosly defined setup and model lists

## Last Updated Tuesday 17th March 2020
##

fZooMSS_Project <- function(model){

  # Pull out some useful parameters - just a shortcut
  param <- model$param
  dt <- param$dt
  dx <- param$dx

  idx.iter <- 2:model$ngrid
  idx <- 2:(model$ngrid-1)
  itimemax  <- param$tmax / param$dt  #max index of time array

  if(length(param$zoo_grps) > 1){ # If there's only one zoo group, then you do not need w0idx. All this stuff gives you info about all zoo groups except the smallest zoo group.
    w0idx <- which(param$Groups$W0 > min(param$Groups$W0) & is.na(param$Groups$Prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    props_z <- param$Groups$Prop[w0idx] # Zooplankton proportions

    for(i in 1:length(w0idx)){ # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(model$w), digits = 2) == param$Groups$W0[w0idx[i]])
    }
  }

  # Matrices for MvF and MvF-D numeric solution
  A.iter <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)
  C.iter <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)
  S.iter <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)

  A <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)
  B <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)
  C <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid)
  S <- matrix(0, nrow = model$param$ngrps, ncol = model$ngrid) # Previous abundance

  # Temporary Matrices that get updated each time step some of these saved for output
  N <- matrix(model$N[1,,], nrow = model$param$ngrps, ncol = model$ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes

  pb <- txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3) # Initial progress bar

  # BIG TIME LOOP
  for (itime in 1:itimemax){

    setTxtProgressBar(pb, itime) # Update progress bar

    growth_multiplier <- colSums(N*model$assim_eff)
    predation_multiplier <- N*model$temp_eff
    diffusion_multiplier <- colSums(N*(model$assim_eff^2))

    ### RFH - Apply is slow as it implements a loop. Turns out colSums and aperm is 50 % faster in these cases
    sw <- sweep(model$dynam_growthkernel, 3, growth_multiplier, '*')
    ap <- colSums(aperm(sw, c(3,1,2)))
    gg <- model$ingested_phyto + ap

    sw2 <- sweep(model$dynam_mortkernel, c(2,3), predation_multiplier, '*')
    M2 <- colSums(colSums(aperm(sw2, c(2,3,1))))
    rm(sw, sw2, ap)

    # Total dynamic spectrum mortality
    Z <- sweep(model$M_sb + model$fish_mort, 2, M2, '+')

    sw <- sweep(model$dynam_diffkernel, 3, diffusion_multiplier, '*')
    ap <- colSums(aperm(sw, c(3,1,2)))
    diff <- model$diff_phyto + ap

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

    ### RFH - I have a C++ implementation for this but it won't help with speed at this point I don't think.
    for(i in 1:model$param$ngrps){

      ### RFH - Can't these next few rows be done in fZooMSS_Setup and stored as a vector?
      ## Set size range index for current group
      curr_min_size <- which(round(log10(model$w), digits = 2) == param$Groups$W0[i])
      curr_max_size <- which(round(log10(model$w), digits = 2) == param$Groups$Wmax[i])
      idx_curr <- (curr_min_size+1):curr_max_size

      for(j in idx_curr){ ## Find the abundance at the next size class with standard MvF
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

    if(length(param$zoo_grps) > 1){ # If you only have one zoo group, it will be locked to phyto spectrum so you do not need to do this
      for(i in 1:length(w0idx)){
        w_min_curr <- w0mins[i]
        exclude_mins <- w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] <- props_z[i] * sum(N[-exclude_mins, w_min_curr])
      }
    }

    fish_mins <- unlist(lapply(param$Groups$W0[param$fish_grps],
                              function(x){which(round(log10(model$w), digits = 2) == x)}))

    if(length(param$fish_grps) > 1 & length(param$zoo_grps) > 1){
      N[param$fish_grps,fish_mins] = (1/length(param$fish_grps))*(colSums(N[-param$fish_grps,fish_mins]))
    }else{
      N[param$fish_grps, fish_mins] = (1/length(param$fish_grps))*sum(N[-param$fish_grps, fish_mins])
    }

    if(param$fish_on == FALSE){
      N[param$fish_grps,] <- 0 # switch off fish_groups
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
      ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps
      # dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
      N_array <- aperm(replicate(model$ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(model$param$ngrps, N_array), c(4,1,2,3))

      dynam_diet =  rowSums(aperm(rowSums(sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

      model$diet[isav,,1:3] = phyto_diet
      model$diet[isav,,c(4:(dim(param$Groups)[1]+3))] = dynam_diet

      model$N[isav,,] <- N # Save abundance

      ## Save Abbundance
      model$Abundance[isav,] <- rowSums(model$N[isav,,])

      ## Save biomass
      model$Biomass[isav,] <- rowSums(model$N[isav,,]*matrix(model$w, nrow = model$param$ngrps, ncol = model$ngrid, byrow = TRUE))

      ## Save mortality rates
      model$M2[isav,,] <- M2 # Save predation mortality rates

      ## Save growth
      model$gg[isav,,] <-  model$ingested_phyto + apply(sweep(model$dynam_growthkernel, 3, growth_multiplier, '*'), c(1,2), sum)
    }

  } # End of time loop

  return(model)

}