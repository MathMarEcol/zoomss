# Unified function to visualize time series changes for different metrics

Creates time series plots showing how abundance, biomass, mortality, or
growth rates of functional groups change throughout the ZooMSS
simulation period.

## Usage

``` r
plotTimeSeries(
  mdl,
  by = "abundance",
  type = "line",
  transform = "identity",
  species = NULL
)
```

## Arguments

- mdl:

  ZooMSS results object containing model outputs with time series data

- by:

  Character string specifying the metric to plot. Options: "abundance",
  "biomass", "mortality", "growth" (default: "abundance")

- type:

  Character vector of plot type. Use `line` for the default line plot,
  `stack` or `fill` (as per geom_area) for stacked or proportional
  plots. (default: "line")

- transform:

  Character vector of the required y-axis transformation. Options from
  `scale_*_continuous` (Default: "identity).

- species:

  Character vector of species names to include. If NULL, all species
  included (default: NULL, applies to all metrics)

## Value

ggplot object showing the requested time series by species

## Details

Plot Time Series Data for ZooMSS Results

This function creates time series visualizations by:

- **Abundance**: Summing abundances across size classes, log-transformed
  y-axis

- **Biomass**: Calculating biomass (abundance × weight), with optional
  stacking and proportional scaling

- **Mortality**: Averaging predation mortality rates across size classes

- **Growth**: Averaging growth rates across size classes,
  log-transformed y-axis

All plots use species-specific colors and filter out zero values. Time
series plots help identify:

- Equilibration time for model runs

- Seasonal or cyclical patterns in ecological metrics

- Relative patterns between functional groups

- Model stability and convergence behavior

## Examples

``` r
if (FALSE) { # \dontrun{
# After running ZooMSS model
results <- zoomss_model(input_params, Groups)

# Plot different metrics
abundance_plot <- plotTimeSeries(results, by = "abundance", transform = "log10")
biomass_plot <- plotTimeSeries(results, by = "biomass", transform = "log10")
mortality_plot <- plotTimeSeries(results, by = "mortality")
growth_plot <- plotTimeSeries(results, by = "growth")

stacked_plot <- plotTimeSeries(results, by = "biomass", type = "stack")
prop_plot <- plotTimeSeries(results, by = "biomass", type = "fill")

# Focus on specific species (works for all metrics)
copepod_plot <- plotTimeSeries(results, by = "biomass",
                              species = c("OmniCopepods", "CarnCopepods"))
abundance_copepods <- plotTimeSeries(results, by = "abundance",
                                    species = c("OmniCopepods", "CarnCopepods"))
mortality_copepods <- plotTimeSeries(results, by = "mortality",
                                    species = c("OmniCopepods", "CarnCopepods"))
growth_copepods <- plotTimeSeries(results, by = "growth",
                                 species = c("OmniCopepods", "CarnCopepods"))
} # }
```
