# Calculate predator-prey mass ratio data for visualization

Calculates predator-prey mass ratio (PPMR) values and biomass weightings
for creating PPMR distribution plots in ZooMSS analysis.

## Usage

``` r
extractPPMR(mdl)
```

## Arguments

- mdl:

  ZooMSS results object containing abundance data (N) and model
  parameters

## Value

For 2D input: List containing PPMR density data and species-specific
values for plotting For 3D input: Array where first dimension is time,
containing PPMR results for each timestep

## Details

Calculate PPMR Data for Plotting

This function computes theoretical and realized PPMR patterns by:

- Calculating size-dependent PPMR values using Wirtz 2012 equations

- Weighting by biomass to show community-level patterns

- Computing species-specific PPMR values

- Handling special cases for filter feeders (larvaceans, salps)

- Processing either time-averaged abundances (2D) or full time series
  (3D)

The function dynamically determines size class ranges for larvaceans and
salps based on their W0 and Wmax values. For 3D abundance arrays, it
calculates PPMR for each time step separately.

This is a helper function primarily used by plotPPMR for visualization.
PPMR analysis provides insights into food web structure and predation
patterns.
