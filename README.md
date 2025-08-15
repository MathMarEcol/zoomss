
<!-- README.md is generated from README.Rmd. Please edit that file -->

# zoomss

<!-- badges: start -->

<!-- badges: end -->

## Installation

You can install the development version of zoomss from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("MathMarEcol/zoomss")
```

# Zooplankton Model of Size Spectrum (ZooMSS)

## Overview of ZooMSS

The Zooplankton Model of Size Spectra (ZooMSS) is a functional
size-spectrum model of the marine ecosystem (following Heneghan et
al. 2016) to resolve phytoplankton, nine zooplankton functional groups
(heterotrophic flagellates and ciliates, omnivorous and carnivorous
copepods, larvaceans, euphausiids, salps, chaetognaths and jellyfish)
and three size-based fish groups. Zooplankton functional groups are
resolved using their body-size ranges, size-based feeding
characteristics and carbon content, and the zooplankton community
emerges from the model across global environmental gradients, depending
on the functional traits of the different groups.

We developed the Zooplankton Model of Size Spectra (ZooMSSv2) based on
the prototype of Heneghan et al. (2016). ZooMSS uses the functional
size-spectrum framework (Blanchard et al., 2017) to resolve the body
size ranges, size-based feeding characteristics and carbon content of
nine zooplankton groups and three fish groups. The model supports
time-varying environmental conditions enabling studies of seasonal
cycles, climate change scenarios, and ecosystem responses to
environmental variability.

ZooMSS represents the marine ecosystem as three communities:
phytoplankton, zooplankton and fish. The zooplankton community consists
of nine of the most abundant zooplankton groups, and the fish community
was made up of a small, medium and large group. Dynamics of the
phytoplankton are not explicitly resolved in the model, rather the mean
size structure of the phytoplankton community is estimated directly from
satellite chlorophyll a observations (Brewin et al., 2010; Barnes et
al., 2011; Hirata et al., 2011). Abundances of the zooplankton and fish
communities are driven by size-dependent processes of growth and
mortality, with the temporal dynamics of each functional group governed
by separate second-order McKendrick-von Foerster equations.

## How to run ZooMSS using `R`

### Input Files

ZooMSS requires two input files:

1.  **GroupInputs.csv** - Contains all taxa-specific parameter values,
    including size ranges and functional group properties.

2.  **Environmental data** - A data frame with time series of
    environmental conditions:

    - For time-varying conditions: Time series with `time_step`, `sst`,
      and `chlo` columns
    - For constant conditions: Single values replicated across all time
      steps

### Running the Model

#### Constant Environmental Conditions

For studies with unchanging environmental conditions:

``` r
# Create constant environmental conditions
enviro_data <- data.frame(
  time_step = 1:5000,
  sst = rep(15, 5000),  # Constant temperature
  chlo = rep(0.5, 5000) # Constant chlorophyll
)

# Include in your setup parameters
input_params <- cbind(base_params, enviro_data)
source("setup_RUN_NAME.R")
```

#### Time-varying Environmental Conditions

For dynamic environmental studies:

``` r
# Create time-varying environmental conditions
enviro_data <- data.frame(
  time_step = 1:5000,
  sst = 15 + 4 * sin(2 * pi * (1:5000) / 365),  # Seasonal cycle
  chlo = 0.5 + 0.3 * cos(2 * pi * (1:5000) / 365)
)

# Include in your setup parameters
input_params <- cbind(base_params, enviro_data)
```

### Environmental Data Structure

Provide environmental time series as a data frame with: - `time_step`:
Sequential time step number (1, 2, 3, …) - `sst`: Sea surface
temperature (°C) - `chlo`: Chlorophyll concentration (mg/m³)

The model processes all environmental conditions dynamically, whether
they vary through time or remain constant.

### Plotting

The model includes several built-in plotting functions for analysis and
visualization:

#### Basic Plots

- `zPlot_BiomassTimeSeries(dat)` - Plot biomass time series with
  optional species selection and proportional stacking
- `zPlot_SizeSpectra(dat)` - Plot overall species-resolved size
  spectra  
- `zPlot_PPMR(dat)` - Plot predator-prey mass ratios

#### Time Series Analysis (when SaveTimeSteps = TRUE)

- `zPlot_AbundTimeSeries(dat)` - Plot abundance time series
- `zPlot_GrowthTimeSeries(dat)` - Plot growth rate time series
- `zPlot_PredTimeSeries(dat)` - Plot predation time series

where `dat` is the model output list.

## Performance and Computational Considerations

## Publications

1.  Heneghan, R.F., Everett, J.D., Blanchard, J.L., Richardson,
    A.J., 2016. Zooplankton Are Not Fish: Improving Zooplankton Realism
    in Size-Spectrum Models Mediates Energy Transfer in Food Webs.
    Front. Mar. Sci. 3, 1–15. <https://doi.org/10.3389/fmars.2016.00201>

2.  Heneghan, R.F., Everett, J.D., Sykes, P., Batten, S.D., Edwards, M.,
    Takahashi, K., Suthers, I.M., Blanchard, J.L., Richardson, A.J., in
    review, A global size-spectrum model of the marine ecosystem that
    resolves zooplankton composition. Ecological Modelling

## Getting Help

If you encounter problems running the model, raise an issue on GitHub:
<https://github.com/MathMarEcol/ZoopSizeSpectraModel/issues>

If you find errors or want to improve the model, we’d love you to make
the changes and submit a pull request for us to review and approve.
