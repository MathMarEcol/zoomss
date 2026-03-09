# Convert ZooMSS abundance matrices to biomass by multiplying by body weights

Converts abundance data to wet weight biomass by multiplying abundances
by the corresponding body weights for each size class. Optionally
converts to carbon biomass.

## Usage

``` r
getBiomass(mdl, units = "ww")
```

## Arguments

- mdl:

  ZooMSS model object containing abundance array (N) and weight vector
  (param\$w)

- units:

  Character string specifying biomass units: "ww" for wet weight
  (default) or "carbon" for carbon biomass

## Value

3D array of biomass values with same dimensions as N. Units depend on
the units parameter:

- "ww": grams wet weight

- "carbon": grams carbon

## Details

Convert Abundance to Biomass

This function transforms abundance matrices to biomass by applying the
weight vector across size classes. Essential for analyses requiring
biomass units rather than abundance counts. Works with 3D arrays (time,
groups, size_classes). Can convert to either wet weight or carbon
biomass units.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run ZooMSS model
results <- zoomss_model(input_params, Groups)

# Convert abundances to wet weight biomass
biomass_ww <- getBiomass(results, units = "ww")

# Convert abundances to carbon biomass
biomass_carbon <- getBiomass(results, units = "carbon")
} # }
```
