## Calculates average trophic level of each group in each grid square, returns a matrix where
## rows correspond to grid squares, columns to each group

Extract_ZooMSS_Species <- function(N, TestGroups, w){
  species <- N %>%
    reduce(`+`) # Elementwise summing of the model abundance
  species <- species/length(N)

  rownames(species) <- TestGroups$species
  out <- as_tibble(t(species))
  out <- out %>%
    add_column("Weight" = w)
}

# Return Abundance/Biomass by weight class

Extract_ZooMSS <- function(N, w){

  Abundance <- matrix(NA, nrow = length(N), ncol = length(w)) # Preallocate
  Biomass <- matrix(NA, nrow = length(N), ncol = length(w)) # Preallocate

  # Loop through and extract the model output.
  for (i in 1:length(N)) {
    Biomass[i,] <- apply(sweep(N[[i]], 2, w, '*'), 2, sum) # Use sweep to multiply by the weights
    Abundance[i,] <- apply(N[[i]], 2, sum)
  }
  out <- tibble("Biomass" = Biomass, "Abundance" = Abundance)

}


NBSS_data <- function(LOPC, Bins, chl_lims){
  LOPC2 <- LOPC %>%
    dplyr::filter(Chl >= chl_lims[1] & Chl <= chl_lims[2]) %>%
    dplyr::select(starts_with("NB_"), Chl) %>%
    pivot_longer(cols = starts_with("NB_")) %>%
    mutate(name = str_remove(name,"NB_"),
           name = as.numeric(name))

  Bins2 <- Bins[LOPC2$name,]
  colnames(Bins2) <- "Bins"
  LOPC2 <- bind_cols(LOPC2, Bins2)
  LOPC2 <- LOPC2 %>%
    drop_na(value) %>%
    dplyr::filter(value > 0) %>%
    dplyr::rename(NB = value)

  NBSS <- LOPC2 %>%
    group_by(name) %>%
    summarise(NB = mean(NB, na.rm = T),
              Bins = Bins[1],
              n = n())

  rm(Bins, Bins2)
  out <- list()
  out[[1]] <- NBSS
  out[[2]] <- LOPC2
  return(out)

  # return(NBSS)
}



NBSS_ZooMSS <- function(res, chl_lims, enviro_data, w, plot_params){

  w_mg <- w *1e3 # Convert weight to mg (as is normal for zooplankton NBSS)
  subs_mdl <- enviro_data$chlo >= chl_lims[1] & enviro_data$chlo <= chl_lims[2] # What cells are within the chl range?
  limits = 10^ c(log10(w_mg[1])-plot_params$dx/2, log10(w_mg) + plot_params$dx/2) # Set limits of each bin

  model_out <- Extract_ZooMSS(res[subs_mdl], w_mg)

  model_NBSS <- tibble("Bins" = w*1000)

  model_NBSS <- model_NBSS %>%
    mutate(
      Biomass = colSums(model_out$Biomass) / length(res), # Get average of all biomass rows
      Abundance = colSums(model_out$Abundance) / length(res), # Get average of all abundance rows
      BinLower = limits[1:length(limits)-1], # Lower limit of the bins
      BinUpper = limits[2:length(limits)], # Upper limit of the bins
      BinWidth = BinUpper - BinLower,
      NB = Biomass/BinWidth) %>%
    filter(log10(BinUpper) <= plot_params$MaxSize & log10(BinLower) > plot_params$MinSize) %>%
    filter(Biomass > 0)

  rm(subs_mdl, model_out)
  # return(model_NBSS)

  out <- list()
  out[[1]] <- model_NBSS
  out[[2]] <- model_out

  return(out)
}




TrophicLevel_Calc <- function(diets, TestGroups){

  ################### TROPHIC LEVEL
  stor = matrix(0, nrow = length(diets), ncol = 12) # Output storage, I have 12 groups, so 12 columns and 1 row for each grid square
  colnames(stor) <- TestGroups$species

  starter_tl = c(1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)### Starting trophic levels, phyto is 1, for all other 12 groups is 2

  for(i in 1:length(diets)){ # Loop over number of grid squares
    start_tl = starter_tl # Initialise trophic level vector
    diet = diets[[i]] # Pull out diet matrix for current grid square

    ### You only need these two lines if diet matrix has a column for small, medium and large phytoplankton
    diet[,3] = rowSums(diet[,c(1:3)])
    diet = diet[,-c(1,2)]

    diet_prop = round((diet/rowSums(diet)), 3) # Calculate diet composition, round to 3 decimal places

    for(k in 1:100){ # Gauss-Siedel iterative loop to calculate trophic levels
      for(j in 1:dim(diet_prop)[1]){ # Loop over each group
        curr_tl =  sum(start_tl*diet_prop[j,], na.rm = TRUE) # Calculate trophic level of group j
        start_tl[j+1] = curr_tl + 1 # Update trophic level
      }
    }
    stor[i,] = start_tl[2:13] # Save non phyto trophic levels

  } # Close enviro grid square loop

  return(stor)
}


biom_ratios <- function(res, TestGroups, enviro_data){

  biom_store = matrix(0, nrow = length(res), ncol = 3)
  colnames(biom_store) = c('Phyto', 'Zoo', 'Fish')
  ratio_store = matrix(0, nrow = length(res), ncol = 3)
  colnames(ratio_store) = c('ZooPhyto', 'FishZoo', 'FishPhyto')

  for(i in 1:length(res)){

    ### CALCULATE PHYTO BIOMASS
    # w0_phyto = -14.5		# minimum phytoplankton size class (1um)
    wMax_phyto = min(enviro_data$phyto_max[i])
    wMax_phyto = min(-8, wMax_phyto)
    w_phyto = 10^(seq(-14.5, wMax_phyto,0.1))
    a = enviro_data$a[i]
    b = enviro_data$b[i]
    biom_store[i,1] = sum((10^a)*(w_phyto^(b+1)))

    ### CALCULATE ZOO AND FISH BIOMASS
    w0 = min(TestGroups$W0)
    wMax = max(TestGroups$Wmax)
    w = 10^seq(w0, wMax, 0.1) # Vector of body sizes
    w[w > 1e5] = 0 # Cut off large fish values (heavy senescence mortality on these size classes - model artifact)

    curr_abund = res[[i]]
    curr_biom = rowSums(sweep(curr_abund, 2, w, '*')) # Calculate biomass from abundance and sum across each group

    zoo_groups = which(is.na(TestGroups$prop) == FALSE) # Which rows are zooplankton
    fish_groups = which(is.na(TestGroups$prop) == TRUE) # Which rows are fish

    biom_store[i,2] = sum(curr_biom[zoo_groups]) # Sum biomass over all zoo groups
    biom_store[i,3] = sum(curr_biom[fish_groups]) # Sum biomass over all fish groups

    ratio_store[i,1] = biom_store[i,2]/biom_store[i,1] # Calculate zoo to phyto biom ratio
    ratio_store[i,2] = biom_store[i,3]/biom_store[i,2] # Calculate fish to zoo biom ratio
    ratio_store[i,3] = biom_store[i,3]/biom_store[i,1] # Calculate fish to phyto biom ratio
  }

  return(ratio_store)

}

Biomass_Change = function(model_abunds, control_abunds, w) {
  model_fish_bioms <- 0 # vector(nrow = length(model_abunds), dim(enviro_data)[1]) # Preallocate Fish Biomass
  control_fish_bioms <- 0 # vector(nrow = length(control_abunds), dim(enviro_data)[1]) # Preallocate Fish Proportions

  # fish_groups <- 10:12 # Which rows are fish
  control_fish_groups <- ((dim(control_abunds[[1]])[1])-2):dim(control_abunds[[1]])[1]
  model_fish_groups <- ((dim(model_abunds[[1]])[1])-2):dim(model_abunds[[1]])[1]

  for(j in 1:length(model_abunds)){
    control_fish_bioms[j] = sum(apply(sweep(control_abunds[[j]], 2, w, '*'), 1, sum)[control_fish_groups]) # Get the reference/control data
    model_fish_bioms[j] = sum(apply(sweep(model_abunds[[j]], 2, w, '*'), 1, sum)[model_fish_groups]) # The run we are comparing to
  }

  # fish_props = (model_fish_bioms/control_fish_bioms)*100-100
  fish_props = ((model_fish_bioms-control_fish_bioms)/control_fish_bioms)*100
  return(fish_props)
}

Pie_Preparation = function(bdf){
  bdf <- bdf %>%
    dplyr::select(-c(Fish_Small, Fish_Med, Fish_Large, SST, Chl, x, y))

  ZS <- bdf %>%
    mutate(ZoopSum = dplyr::select(., Larvaceans:Jellyfish) %>%
             apply(1, sum, na.rm=TRUE))

  bdf <- bdf/ZS$ZoopSum
  bdf <- bdf %>%
    pivot_longer(names_to = "Taxa", values_to = "Biomass", cols  = everything()) %>%
    group_by(Taxa) %>%
    summarise(Mean = round(mean(Biomass,na.rm = TRUE), digits = 2)) %>%
    ungroup %>%
    mutate(Taxa = factor(Taxa))
  rm(ZS)

  bdf <- left_join(x = bdf, y = TestGroups[ , c("species", "Colour", "Feeding")], by = c("Taxa"="species"), all.x=TRUE)
  bdf <- bdf %>%
    arrange(Feeding) %>% # Group feeding types together
    mutate(Taxa = factor(Taxa, levels = Taxa), # Fix factor levels to this order too
    )
}

Biomass_Proportion = function(){
  biomass_brick <- readRDS("../GlobalZoopAbundModels/ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")

  bdf <- as.data.frame(biomass_brick, xy = TRUE)

  bdf <- bdf %>%
    dplyr::select(-c(Fish_Small, Fish_Med, Fish_Large)) %>%
    mutate(Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

  bdf2 <- pivot_longer(bdf, names_to = "Taxa", values_to = "Biomass",
                       cols = c("Larvaceans", "OmniCopepods","CarnCopepods", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish"))

  bdf3 <- bdf2 %>%
    filter(!is.na(Chl) & !is.na(SST)) %>%
    mutate(Chl = round(Chl,1), # Round to 2 decimal places
           Chl = as.factor(Chl)) %>%
    dplyr::group_by(Chl)

  bdf4 <- bdf3 %>%
    summarise(TotalBio = sum(Biomass,na.rm = T),
              Larvaceans = sum(Biomass[Taxa == "Larvaceans"],na.rm = T)/TotalBio,
              OmniCopepods = sum(Biomass[Taxa == "OmniCopepods"],na.rm = T)/TotalBio,
              CarnCopepods = sum(Biomass[Taxa == "CarnCopepods"],na.rm = T)/TotalBio,
              Euphausiids = sum(Biomass[Taxa == "Euphausiids"],na.rm = T)/TotalBio,
              Chaetognaths = sum(Biomass[Taxa == "Chaetognaths"],na.rm = T)/TotalBio,
              Salps = sum(Biomass[Taxa == "Salps"],na.rm = T)/TotalBio,
              Jellyfish = sum(Biomass[Taxa == "Jellyfish"],na.rm = T)/TotalBio) %>%
    ungroup()

  bdf5 <- gather(bdf4, key = "Taxa", value = "Proportion", -c(TotalBio, Chl)) %>%
    mutate(Chl = as.numeric(as.character(Chl)),
           Taxa = as.factor(Taxa),
           Taxa = factor(Taxa, levels = c("OmniCopepods", "Euphausiids", "Jellyfish", "Chaetognaths", "CarnCopepods", "Salps", "Larvaceans")))
}


PhytoBio_Proportion = function(){
  biomass_brick <- readRDS("../GlobalZoopAbundModels/ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")

  source("../ZooMSS/Ancillary/Supplementary/EnvironmentalData/fZooMSS_CalculatePhytoParam.R")

  bdf <- as.data.frame(biomass_brick, xy = TRUE)

  phyto <- fZooMSS_CalculatePhytoParam(10^bdf$Chl)
  bdf <- bind_cols(bdf, phyto)

  bdf <- bdf %>%
    dplyr::select(c(Flagellates, Ciliates, Chl, SST, pico_biom, nano_biom, micro_biom)) %>%
    rename("PicoPhytoplankton" = pico_biom, "NanoPhytoplankton" = nano_biom, "MicroPhytoplankton" = micro_biom) %>%
    mutate(Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

  bdf2 <- pivot_longer(bdf, names_to = "Taxa", values_to = "Biomass",
                       cols = c(Flagellates, Ciliates, PicoPhytoplankton, NanoPhytoplankton, MicroPhytoplankton))

  bdf3 <- bdf2 %>%
    filter(!is.na(Chl) & !is.na(SST)) %>%
    mutate(Chl = round(Chl,1), # Round to 2 decimal places
           Chl = as.factor(Chl)) %>%
    dplyr::group_by(Chl)

  bdf4 <- bdf3 %>%
    summarise(TotalBio = sum(Biomass,na.rm = T),
              Flagellates = sum(Biomass[Taxa == "Flagellates"],na.rm = T)/TotalBio,
              Ciliates = sum(Biomass[Taxa == "Ciliates"],na.rm = T)/TotalBio,
              PicoPhytoplankton = sum(Biomass[Taxa == "PicoPhytoplankton"],na.rm = T)/TotalBio,
              NanoPhytoplankton = sum(Biomass[Taxa == "NanoPhytoplankton"],na.rm = T)/TotalBio,
              MicroPhytoplankton = sum(Biomass[Taxa == "MicroPhytoplankton"],na.rm = T)/TotalBio) %>%
    ungroup()

  bdf5 <- gather(bdf4, key = "Taxa", value = "Proportion", -c(TotalBio, Chl)) %>%
    mutate(Chl = as.numeric(as.character(Chl)),
           Taxa = as.factor(Taxa),
           Taxa = factor(Taxa, levels = c("PicoPhytoplankton", "NanoPhytoplankton", "MicroPhytoplankton", "Flagellates", "Ciliates")))
}




Diet_Prepare = function(diets, enviros){

  diets <- diets[enviros]

  ### These lines create diet matrices for subset (oligo, eutro etc). The matrices rows are predators, columns are prey.
  ### The first 3 columns are phytoplankton broken into pico, nano and microphytoplankton
  diets_mat <- apply(array(unlist(diets), dim = c(nrow(diets[[1]]), ncol(diets[[1]]), length(diets))), c(1,2), mean)

  ### Add together pico, nano, microphytoplankton, append to diet matrices and remove pico, nano and microphytoplankton columns
  diets_mat <- cbind(rowSums(diets_mat[,c(1:3)]), diets_mat[,-c(1:3)])

  ### Add dummy row for phytoplankton predation (all zero)
  diets_mat <- rbind(rep(0, 13), diets_mat)

  ### Flip matrices, so rows are prey and columns are predators
  diets_mat = t(diets_mat)

  namess <- c("Phytoplankton", "Flagellates", "Ciliates", "Larvaceans", "OmniCopepods",
              "CarnCopepods", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish",
              "PlanktivorousFish", "Fish_Med", "Fish_Large")
  colnames(diets_mat) <- namess
  rownames(diets_mat) <- namess

  return(diets_mat)
}


Chord_Prepare = function(diets, enviros, TestGroups) {

  diets <- diets[enviros]

  ### These lines create diet matrices for subset (oligo, eutro etc). The matrices rows are predators, columns are prey.
  ### The first 3 columns are phytoplankton broken into pico, nano and microphytoplankton
  diets_mat <- apply(array(unlist(diets), dim = c(nrow(diets[[1]]), ncol(diets[[1]]), length(diets))), c(1,2), mean)

  ### Add together pico, nano, microphytoplankton, append to diet matrices and remove pico, nano and microphytoplankton columns
  diets_mat <- cbind(rowSums(diets_mat[,c(1:3)]), diets_mat[,-c(1:3)])

  ### Add dummy row for phytoplankton predation (all zero)
  diets_mat <- rbind(rep(0, 13), diets_mat)

  ### Flip matrices, so rows are prey and columns are predators, then add together fish groups to get total fish diet
  diets_mat = t(diets_mat)
  diets_mat[,11] = rowSums(diets_mat[,c(11:13)])
  diets_mat[11,] = colSums(diets_mat[c(11:13),])
  diets_mat = diets_mat[1:11, 1:11] ### Remove seperated fish groups

  namess <- c("Phytoplankton", TestGroups$species[1:9], "Fish")
  colnames(diets_mat) <- namess
  rownames(diets_mat) <- namess

  return(diets_mat)
}

PPMR_plot = function(TestGroups, res, enviros){

  min_size = min(TestGroups$W0) # smallest size class
  max_size = max(TestGroups$Wmax) # largest size class
  w = 10^(seq(from = min_size, to = max_size, 0.1)) # all size classes

  # Calculate PPMR (beta) table, where dim1 = group, dim2 = body size with
  # value being PPMR for that body size (this is not realised PPMR - not
  # emergent from diet but calculated from m-values and Wirtz, 2012 equation)
  D.z = 2*(3*(w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
  zoo_m = TestGroups$m # pull out m-values from parameter table
  betas =  log10(t(sapply(zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3}))) # Convert m to betas, using Wirtz 2012 equation
  betas = betas[-which(is.na(TestGroups$m)),] # remove fish rows

  ## Modify beta matrix for larvaceans and salps - all size classes for these groups feed on same prey, so log10PPMR increases by 0.1 for each 0.1 log10 size interval
  betas[3,45:75] <- betas[3,44] + seq(0.1,3.1,0.1) # Larvaceans (index 44 in w vector is smallest size class, 75 is maximum size class)
  betas[8,61:121] <- betas[8,61] + seq(0.1,6.1,0.1) # Salps (index 61 in w vector is smallest size class, 121 is maximum size class

  # Calculate ave abundances across oligo/eutro grid squares, then calculate ave
  # biomass and proportion of total zoo biomass that is from each group size class
  abunds = res[enviros]
  ave = matrix(0, nrow = dim(TestGroups)[1], ncol = length(w))
  for(i in 1:length(abunds)){
    ave = ave + abunds[[i]]/length(abunds)
  }
  ave_biom = sweep(ave, 2, w, "*") # Calculate oligo biomass for zoo groups
  ave_biom = ave_biom[-which(is.na(TestGroups$m)),] # remove rows for fish
  beta_props = ave_biom/sum(ave_biom) # Calculate fraction of zoo biomass in each group, in each size class

  out <- list()
  out[[1]] <- betas
  out[[2]] <- beta_props
  names(out) <- c("betas", "beta_props")

  temp <- density(betas, weights = beta_props)

  out <- tibble("x" = temp$x, "y" = temp$y, "mn_beta" = sum(beta_props*betas))

  spbeta_props = ave_biom/rowSums(ave_biom) # Species specific proportions
  spPPMR <- tibble("Species" = TestGroups$species[-which(is.na(TestGroups$m))], "Betas" = rowSums(spbeta_props*betas), "y" = NA) # Get species-specific PPMR

  for (s in 1:length(spPPMR$Species)){
    spPPMR$y[s] <- out$y[which.min(abs(out$x - spPPMR$Betas[s]))]
  }
  out2 <- list()
  out2[[1]] <- out
  out2[[2]] <- spPPMR

  return(out2)
}


### BIOM DATAFRAME
biom_mat_zoo <- function(N, groups, cut_point1, cut_point2){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton (excluding flag and cils)
  num_zoo = length(zoo_groups)

  biom_mat <- matrix(0, nrow = length(N), ncol = num_zoo)


  for(i in 1:length(N)){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(w),2) >= cut_point1 & round(log10(w),2) <= cut_point2)
    zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props
  }

  colnames(biom_mat) = as.character(groups$species)[3:9]
  biom_mat = as.data.frame(biom_mat)
  biom_mat
}

# ### BIOM DATAFRAME - Flagellates and Ciliates
# biom_mat_hetero <- function(res, model, TestGroups, cut_point1, cut_point2){
#   zoo_groups = which(is.na(TestGroups$prop) == FALSE)[c(1,2)] # Which rows are flagellates and ciliates
#   num_zoo = length(zoo_groups)+1 # Add one column for the phytoplankton
#
#   biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
#
#   for(i in 1:dim(enviro_data)[1]){
#     # Calculate flag and cil biomass
#     N_ave = res[[i]]
#     B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
#     weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
#     zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
#
#     # Calculate phytoplankton biomass
#     w0_phyto = -14.5    # minimum phytoplankton size class (1um)
#     wMax_phyto = min(enviro_data$phyto_max[i])
#     w_phyto = 10^(seq(-14.5, wMax_phyto,0.1))
#     a = enviro_data$a[i]
#     b = enviro_data$b[i]
#     phyto_biom = sum((10^a)*(w_phyto^(b+1)))
#
#     # Concatenate and save
#     all_bioms = c(phyto_biom,zoo_bioms)
#     biom_mat[i,] <- all_bioms
#   }
#
#   colnames(biom_mat) = c('Phyto',as.character(TestGroups$species)[1:2])
#   biom_mat = as.data.frame(biom_mat)
# }

### RATIO OF ZOO TO PHYTO, FISH TO ZOO AND FISH TO PHYTO
ratio_frame <- function(N, groups, enviro, cut_point1, cut_point2){
  enviro_data = enviro
  chlo <- enviro_data$chlo
  c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio (0.02 Chl:C)
  phyto_carb <- c_chl*chlo/1000 # (convert to grams carbon)
  # phyto_ww <- phyto_carb*0.1


  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1:2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)

  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)

  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props
  }

  colnames(biom_mat) = as.character(groups$species)[3:9]
  biom_mat = as.data.frame(biom_mat)
  zoo_ww = rowSums(biom_mat)
  zoo_carb = rowSums(sweep(biom_mat, 2, groups$carbon[3:9], "*"))

  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  num_fish = length(fish_groups)

  biom_matf <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_fish)

  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    fish_bioms = rowSums(B_ave[fish_groups,])
    fish_props = fish_bioms
    biom_matf[i,] <- fish_props
  }

  colnames(biom_matf) = as.character(groups$species)[10:12]
  biom_matf = as.data.frame(biom_matf)
  fish_ww = rowSums(biom_matf)
  fish_carb = rowSums(sweep(biom_matf, 2, groups$carbon[10:12], "*"))

  ratio_frame = data.frame("ZP" = zoo_carb/phyto_carb, "FZ" = fish_carb/zoo_carb, "FP" = fish_carb/phyto_carb)
  ratio_frame

}

## CALCULATE SLOPES OF PHYTO, ZOO AND FISH
slope_frame <- function(N, groups, start, finish){
  ########## RESULTS TABLES
  slopes <- matrix(0, dim(enviro_data)[1],3)
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  num_fish = length(fish_groups)

  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    enviro <- enviro_data[i,]
    tot_zoo = colSums(N_ave[zoo_groups,])
    tot_fish = N_ave[fish_groups,]

    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }
    tot_anim = tot_zoo + tot_fish
    # Average Slope
    fish_start = which(round(log10(model$w), digits = 2) == (param$groups$W0[dim(groups)[1]]))
    fish_finish = which(round(log10(model$w), digits = 2) == (param$groups$Wmat[dim(groups)[1]]))
    max_phyto = round(log10(param$wMax_phyto), digits = 2)
    zoo_start =  which(round(log10(model$w), digits = 2) == start)
    zoo_finish = which(round(log10(model$w), digits = 2) == finish)
    zoo_slope2 = lm(log10(tot_anim[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$
      coefficients[2]
    #zoo_slope2 =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
    #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
    fish_slope2 = round(lm(log10(tot_fish[fish_start:fish_finish])~log10(model$w[fish_start:fish_finish]))$
                          coefficients[2], digits = 2)
    #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
    #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
    phyto_slope = enviro$b

    slopes[i,1] <- phyto_slope
    slopes[i,2] <- zoo_slope2
    slopes[i,3] <- fish_slope2
  }

  colnames(slopes) <- c('Phyto_Slope', 'Zoo_Slope', 'Fish_Slope')
  slopes
}

### CALCULATE TOTAL FISH BIOMASS IN EACH GRID SQUARE
biom_mat_fish = function(res){
  l_res = length(res)
  fish_biom = rep(0, l_res)

  for(i in 1:l_res){
    N = res[[i]]
    fish = N[c(10:12),c(1:168)]
    w_fish = w[1:168]
    fish_biom[i] = sum(sweep(fish, 2, w_fish, "*"))
  }
  fish_biom
}

#first, take out the by-species total abundances for each grid square
summary_by_species <- function(N, Groups, enviro_data){
  #setup matrix with a row per square, and enough columns for the enviro data and species abundance/biomass
  res_summary <- matrix(nrow=dim(enviro_data)[1],ncol=(dim(enviro_data)[2]+2*dim(Groups)[1]))
  colnames(res_summary) <- c(colnames(enviro_data),paste0(as.character(Groups$species),"_abund"),paste0(as.character(Groups$species),"_biom"))
  #paste in the environment data
  for (i in 1:dim(enviro_data)[1]) {
    for (j in 1:dim(enviro_data)[2]) {
      res_summary[i,j]=enviro_data[i,j]
    }
  }
  #abundance
  for (j in 1:dim(Groups)[1]) {
    for (i in 1:dim(enviro_data)[1]) {
      res_summary[i,dim(enviro_data)[2]+j]<-sum(N[[i]][j,])
    }
  }
  #biomass
  dx = 0.1         # log10 weight step
  w0 = 10^(min(Groups$W0))    # minimum dynamic size class
  wMax = 10^(max(Groups$Wmax))# maximum dynamic size class
  w = 10^(seq(from = log10(w0), to =  log10(wMax), dx))
  for (i in 1:dim(enviro_data)[1]) {
    res_summary[i,dim(enviro_data)[2]+dim(Groups)[1]+1:dim(Groups)[1]]<-apply(sweep(N[[i]], 2, w, '*'),1,sum)
  }


  #Biomass <-  apply(sweep(N[[i]], 2, model$w, '*'),c(1,2),sum)

  res_summary <- as_tibble(res_summary)
  return(res_summary)
}
