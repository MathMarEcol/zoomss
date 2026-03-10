# Visualize abundance size spectra across functional groups

Creates a log-log plot of abundance versus body size for all functional
groups, showing the classic size spectrum pattern in marine ecosystems.

## Usage

``` r
plotSizeSpectra(mdl, by = "abundance", n_years = 10)
```

## Arguments

- mdl:

  ZooMSS results object containing model outputs and parameters

- by:

  Character string specifying the metric to plot. Options: "abundance",
  "biomass", "mortality", "growth" (default: "abundance")

- n_years:

  The number of years (from the end) over which to average the size
  spectra

## Value

ggplot object showing log abundance vs log body weight by species

## Details

Plot Size Spectra for ZooMSS Results

This function visualizes the abundance size spectrum by:

- Converting abundance data to long format with body weights

- Filtering out zero abundances to focus on active size classes

- Creating log-log plots colored by functional group

- Using species-specific colors defined in the Groups parameter table

Size spectra are fundamental patterns in marine ecology, typically
showing declining abundance with increasing body size. This
visualization helps assess model realism and identify dominant size
classes within each functional group.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running ZooMSS model
results <- zoomss_model(input_params, Groups)
size_plot <- plotSizeSpectra(results)
print(size_plot)
} # }
```
