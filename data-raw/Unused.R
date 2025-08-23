#' Calculate Species-Level Biomass
#'
#' @title Sum biomass across size classes for each functional group
#' @description Converts abundance to biomass and sums across all size classes
#'   to provide total biomass per functional group per spatial cell.
#' @details This function combines abundance-to-biomass conversion with size-class
#'   aggregation in one step, providing species-level biomass summaries useful
#'   for spatial analyses and ecological comparisons.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector (param$w)
#'
#' @return List of vectors with total biomass per functional group (grams wet weight)
#' @export
#'
zExtractBiomassSpecies = function(res, vmdl) {
  # if (dim(res[[1]])[2] != length(mdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) apply(sweep(x, 2, vmdl$param$w, '*'), 1, sum))
  return(Biomass)
}




#' Convert Abundance to Carbon Biomass
#'
#' @title Convert ZooMSS abundances to carbon biomass across all size classes
#' @description Converts abundance data to carbon biomass by multiplying by body weights
#'   and then by carbon content factors for each functional group.
#' @details This function performs a two-step conversion:
#'   1. Abundance to wet weight biomass (using body weights)
#'   2. Wet weight to carbon biomass (using group-specific carbon content)
#'
#'   Carbon biomass is essential for biogeochemical analyses and comparisons
#'   with field data that are often reported in carbon units.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector and carbon content factors
#'
#' @return List of carbon biomass matrices (grams carbon)
#' @export
#'
zCarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- purrr::map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  return(Biomass)
}

#' Convert Abundance to Species-Level Carbon Biomass
#'
#' @title Convert abundances to carbon biomass and sum across size classes
#' @description Converts abundance data to carbon biomass and then sums across all
#'   size classes to provide total carbon biomass per functional group.
#' @details This function combines carbon biomass conversion with size-class aggregation:
#'   1. Converts abundance to wet weight biomass
#'   2. Converts to carbon biomass using group-specific factors
#'   3. Sums across all size classes for each functional group
#'
#'   Provides species-level carbon biomass useful for ecological stoichiometry
#'   and biogeochemical cycle analyses.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector and carbon content factors
#'
#' @return List of vectors with total carbon biomass per functional group (grams carbon)
#' @export
#'
zSpeciesCarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- purrr::map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  Biomass <- reduceSize(Biomass)

  return(Biomass)
}

#' Calculate Size-Class Biomass
#'
#' @title Sum biomass across functional groups for each size class
#' @description Converts abundance to biomass and sums across all functional groups
#'   to provide total biomass per size class per spatial cell.
#' @details This function provides size-class-level biomass by summing across
#'   functional groups. Useful for analyzing community size structure and
#'   comparing size spectrum patterns between locations.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param w Vector of body weights for each size class (grams)
#'
#' @return List of vectors with total biomass per size class (grams wet weight)
#' @export
#'
zExtractBiomassSize = function(res, w) {
  if (dim(res[[1]])[2] != length(w)){print("error")}
  Biomass <- purrr::map(res, function(x) apply(sweep(x, 2, w, '*'), 2, sum))
  return(Biomass)
}




#' Convert List to Tibble Format
#'
#' @title Convert ZooMSS list output to tibble with species names
#' @description Converts ZooMSS list output to a tibble format with proper column
#'   names based on functional group species names.
#' @details This function converts aggregated ZooMSS output (typically from
#'   reduceSize or similar functions) into a tibble format suitable for
#'   analysis and visualization. Currently designed for 2D data (species x cells).
#'
#' @param li List of vectors/matrices from ZooMSS aggregation functions
#' @param vmdl ZooMSS model object containing species names (param$Groups$Species)
#'
#' @return Tibble with columns named by species and rows representing spatial cells
#' @export
#'
zConvert2Tibble <- function(li, vmdl){
  df <- tibble::as_tibble(matrix(unlist(li), nrow=length(li), byrow=T), .name_repair = "unique") %>%
    dplyr::rename_with(~vmdl$param$Groups$Species)
  return(df)
}





#' Create Diet Matrix in Long Tibble Format
#'
#' @title Convert diet matrix to long format for analysis and visualization
#' @description Converts ZooMSS diet matrix from wide format to long (tidy) format
#'   with predator-prey relationships clearly defined.
#' @details This function transforms diet matrices into a long format suitable for
#'   analysis and visualization of feeding relationships. The resulting tibble
#'   contains predator-prey pairs with diet fraction values, making it easy to
#'   analyze trophic interactions and create food web visualizations.
#'
#' @param mat Diet matrix from ZooMSS output (predators x prey)
#' @param mdl ZooMSS model object containing species names for labelling
#'
#' @return Long tibble with columns: Predator, Prey, Diet
#' @export
#'
makeDietTibble <- function(mdl){

  mdl <- out


  tibble::as_tibble(mdl$diet, .name_repair = "unique")

  suppressMessages(
    out <- tibble::as_tibble(mat, .name_repair = "unique") %>%
      dplyr::rename_with(~c("Phyto_Small", "Phyto_Med", "Phyto_Large", mdl$param$Groups$Species)) %>%
      dplyr::mutate(Predator = mdl$param$Groups$Species) %>%
      tidyr::pivot_longer(cols = .data$Phyto_Small:.data$Fish_Large, names_to = "Prey", values_to = "Diet")
  )
  return(out)
}










