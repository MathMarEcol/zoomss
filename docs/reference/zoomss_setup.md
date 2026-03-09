# Initialize ZooMSS model components and calculate feeding interactions

Sets up the ZooMSS model structure by calculating feeding kernels,
mortality rates, and other model components that remain static during
the simulation.

## Usage

``` r
zoomss_setup(param)
```

## Arguments

- param:

  Complete parameter list created by zoomss_params containing:

  - Groups: Functional group definitions and biological parameters

  - Model dimensions (ngrps, ngrid, time parameters)

  - Environmental forcing time series

  - Physical and biological constants

## Value

Model object containing:

- param: Input parameters (passed through)

- dynam_xxx: Dynamic feeding kernel arrays for group interactions (where
  xxx = growthkernel, diffkernel, dietkernel, mortkernel)

- phyto_xxx: Phytoplankton feeding kernel arrays (where xxx =
  growthkernel, diffkernel, dietkernel)

- nPP: Initial phytoplankton abundance spectrum

- M_sb_base: Baseline senescence mortality rates

- fish_mort: Fishing mortality rates

- assim_eff: Assimilation efficiency matrix

- temp_eff: Temperature effect matrix (initialized)

- N: Initial abundance arrays

- time: Time array for storing time values (initialized as NA)

- Additional model structure components

## Details

Setup ZooMSS Model Structure and Feeding Kernels

This function initializes the core ZooMSS model structure by
calculating:

**Static Components (calculated once):**

- Feeding preference kernels based on predator-prey size ratios

- Search volumes and encounter rates between size classes

- Baseline mortality rates (senescence, fishing)

- Initial abundance distributions for all functional groups

**Dynamic Component Structures (updated during run):**

- Phytoplankton feeding kernels (structure calculated here, values
  updated with environment)

- Growth and diffusion kernels for zooplankton and fish interactions

- Diet and mortality tracking arrays

**Model Architecture:**

- Size-structured populations across logarithmic size classes

- Multiple functional groups with different feeding behaviors

- Environmental coupling through phytoplankton and temperature

The function separates static calculations (done once for efficiency)
from dynamic calculations (updated each time step in zoomss_run).

## Examples

``` r
if (FALSE) { # \dontrun{
# Create parameters for model setup
params <- zoomss_params(Groups, input_params)

# Initialize model structure
model <- zoomss_setup(params)

# Model is now ready for time integration with zoomss_run
results <- zoomss_run(model)
} # }
```
