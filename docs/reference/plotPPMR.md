# Visualize predator-prey mass ratio patterns in ZooMSS results

Creates a plot showing the distribution of predator-prey mass ratios
(PPMR) across functional groups, providing insights into the trophic
structure of the ecosystem.

## Usage

``` r
plotPPMR(mdl, idx)
```

## Arguments

- mdl:

  ZooMSS results object containing model outputs and parameters

- idx:

  The time index to plot

## Value

ggplot object showing PPMR distribution with species-specific overlays

## Details

Plot Predator-Prey Mass Ratio (PPMR)

This function calculates and visualizes PPMR patterns by:

- Computing theoretical PPMR values for each functional group and size
  class

- Weighting by biomass to show realized community patterns

- Creating a density plot of PPMR distribution across the community

- Overlaying species-specific PPMR values as points

PPMR is a key ecological metric that describes the size relationship
between predators and their prey, providing insight into food web
structure and energy transfer efficiency in marine ecosystems.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running ZooMSS model
results <- zoomss_model(input_params, Groups)
ppmr_plot <- plotPPMR(results)
print(ppmr_plot)
} # }
```
