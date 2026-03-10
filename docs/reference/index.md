# Package index

## Data

Data related functions and datasets

- [`GroupInputs`](GroupInputs.md) : Default functional groups for the
  ZooMSS model
- [`getGroups()`](getGroups.md) : Load default or custom functional
  groups for ZooMSS model
- [`calculatePhytoParam()`](calculatePhytoParam.md) : Calculate
  phytoplankton abundance spectrum from chlorophyll data
- [`createEnviroData()`](createEnviroData.md) : Generate synthetic
  environmental data for ZooMSS testing
- [`createInputParams()`](createInputParams.md) : Create input
  parameters data frame for ZooMSS model runs
- [`validateGroups()`](validateGroups.md) : Validate ZooMSS functional
  groups data structure and values

## Model Runs

Functions for running the model

- [`zoomss_model()`](zoomss_model.md) : Main ZooMSS model function for
  complete simulations

## Plotting

Functions for plotting.

- [`plotEnvironment()`](plotEnvironment.md) : Plot environmental forcing
  data
- [`plotPPMR()`](plotPPMR.md) : Visualize predator-prey mass ratio
  patterns in ZooMSS results
- [`plotSizeSpectra()`](plotSizeSpectra.md) : Visualize abundance size
  spectra across functional groups
- [`plotTimeSeries()`](plotTimeSeries.md) : Unified function to
  visualize time series changes for different metrics

## Data wrangling

Helper functions to convert units and data format

- [`averageTimeSeries()`](averageTimeSeries.md) : Calculate mean of
  final portion of ZooMSS time series
- [`extractPPMR()`](extractPPMR.md) : Calculate predator-prey mass ratio
  data for visualization
- [`extractSizeRange()`](extractSizeRange.md) : Extract specific size
  class range from model variable
- [`extractTrophicLevels()`](extractTrophicLevels.md) : Compute trophic
  levels for functional groups using diet composition
- [`getBiomass()`](getBiomass.md) : Convert ZooMSS abundance matrices to
  biomass by multiplying by body weights
- [`getGroups()`](getGroups.md) : Load default or custom functional
  groups for ZooMSS model
- [`reduceAll()`](reduceAll.md) : Aggregate abundances across all groups
  and size classes
- [`reduceSize()`](reduceSize.md) : Aggregate ZooMSS abundances across
  all size classes
- [`reduceSpecies()`](reduceSpecies.md) : Aggregate ZooMSS abundances
  across all species
