# Execute the main ZooMSS simulation loop with dynamic environmental forcing

Runs the ZooMSS model forward in time, updating environmental conditions
and population dynamics at each time step using the McKendrick-von
Foerster framework.

## Usage

``` r
zoomss_run(model)
```

## Arguments

- model:

  Model object created by zoomss_setup containing:

  - param: Complete parameter list with environmental time series

  - Feeding kernels and biological rate parameters

  - Initial conditions and model structure

## Value

List containing complete model output:

- param: Model parameters used in simulation

- N: Abundance time series (time x groups x size classes)

- gg: Growth rate time series

- diet: Diet composition time series

- Z: Mortality rate time series

- time: Time values corresponding to saved results (accounting for
  isave)

- w: Size class weights (g)

- Additional time series data and model results

## Details

Run ZooMSS Model Forward in Time

This is the core simulation engine of ZooMSS that:

**Environmental Dynamics:**

- Updates phytoplankton abundance spectrum based on chlorophyll time
  series

- Applies temperature effects on zooplankton and fish metabolism

- Recalculates feeding kernels with current environmental conditions

**Population Dynamics:**

- Solves McKendrick-von Foerster equation for size-structured growth

- Updates feeding interactions between all size classes and groups

- Calculates mortality from predation, senescence, and fishing

- Handles recruitment and boundary conditions for each functional group

**Time Integration:**

- Processes model through all time steps with adaptive environmental
  forcing

- Saves output at specified intervals for memory efficiency

- Maintains mass balance and numerical stability throughout simulation

Unlike static models, this version dynamically updates phytoplankton
spectra and temperature effects at each time step based on provided
environmental data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set up model parameters and structure
params <- zoomss_params(Groups, input_params)
model <- zoomss_setup(params)

# Run the simulation
results <- zoomss_run(model)

# Access final abundances
final_abundances <- results$N[dim(results$N)[1],,]
} # }
```
