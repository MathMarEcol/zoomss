# Calculate phytoplankton abundance spectrum from chlorophyll data

Converts chlorophyll concentration data to phytoplankton size spectrum
parameters (slope, intercept, maximum size) using established
oceanographic relationships.

## Usage

``` r
calculatePhytoParam(dat)
```

## Arguments

- dat:

  Data frame containing chlorophyll data (chl column in mg/m^3) and
  optionally phytoplankton biomass (phy column in g/m^3)

## Value

Data frame with added columns:

- phyto_slope: Power law slope for phytoplankton size spectrum

- phyto_int: Log10 intercept for phytoplankton abundance

- phyto_max: Maximum phytoplankton size (log10 grams)

- pico_biom, nano_biom, micro_biom: Biomass in each size class

## Details

Calculate Phytoplankton Size Spectrum Parameters

This function implements the Brewin et al. (2015) algorithm to partition
chlorophyll among picophytoplankton, nanophytoplankton, and
microphytoplankton size classes, then calculates:

- Size spectrum slope and intercept parameters

- Maximum phytoplankton size based on micro proportion

- Biomass estimates for each size class

These parameters drive the dynamic phytoplankton spectrum in ZooMSS that
serves as the base of the food web. The function can work with either
chlorophyll-only data (using empirical relationships) or direct
phytoplankton biomass measurements.

## References

Brewin, R.J.W., et al. (2015). A three-component model of phytoplankton
size class for the Atlantic Ocean. Ecological Modelling, 306, 90-101.

Maranon, E., et al. (2014). Resource supply overrides temperature as a
controlling factor of marine phytoplankton growth. PLoS ONE, 9(6),
e99312.
