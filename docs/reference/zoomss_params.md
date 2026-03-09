# Initialize and validate ZooMSS model parameters

Sets up the complete parameter list for ZooMSS model runs, including
functional group parameters, model dimensions, and environmental forcing
data.

## Usage

``` r
zoomss_params(Groups, input_params, isave)
```

## Arguments

- Groups:

  Data frame containing functional group definitions with columns:
  Species, Type, W0 (log min size), Wmax (log max size), and various
  biological parameters

- input_params:

  Data frame with model parameters including: time (time vector in
  years), sst (sea surface temperature), and chl (chlorophyll). The time
  vector can start at any value and the model automatically calculates
  dt (time step) and tmax (maximum time).

- isave:

  Save frequency in time steps (default: 50)

## Value

List containing comprehensive model parameters:

- Groups: Functional group definitions

- ngrps: Number of functional groups

- ngrid: Number of size classes

- w: Size class weights (g)

- tmax, dt, isave: Temporal parameters

- zoo_grps, fish_grps: Indices for different organism types

- phyto_int, phyto_slope: Time series of phytoplankton parameters

- temp_eff_zoo, temp_eff_fish: Time series of temperature effects

- Additional biological and physical parameters

## Details

Set Up ZooMSS Model Parameters

This function creates a comprehensive parameter object that contains:

**Static Parameters (fixed across time steps):**

- Model dimensions (number of groups, size classes, time steps)

- Biological parameters (growth efficiency, mortality rates)

- Size class definitions and ranges for each functional group

- Phytoplankton size spectrum parameters

**Dynamic Parameters (calculated from environmental data):**

- Phytoplankton abundance time series based on chlorophyll

- Temperature effects on metabolism for zooplankton and fish

- Environmental forcing validation and interpolation

The function validates that environmental time series data covers the
full simulation period and pre-calculates time-varying parameters to
optimize model performance during the main simulation loop.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load functional groups
data(Groups)

# Create environmental time series
env_data <- createEnviroData(10, 0.01)
input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

# Generate parameter list
params <- zoomss_params(Groups, input_params, isave = 50)
} # }
```
