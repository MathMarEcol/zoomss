# Plot environmental forcing data

Creates plots of sea surface temperature and chlorophyll time series for
visualizing environmental forcing data used in ZooMSS model runs.

## Usage

``` r
plotEnvironment(env_data)
```

## Arguments

- env_data:

  Environmental data frame with time, sst, chl columns

## Value

ggplot object (if patchwork available) or list of two ggplot objects

## Details

Plot Environmental Time Series

This function creates two separate plots with different y-axes scales:

- SST plot (red line) with temperature in deg C

- Chlorophyll plot (green line) with concentration in mg/m^3

The plots can be combined using the patchwork package if available,
otherwise separate plots are returned as a list. This helps users
visualize the environmental forcing that drives ZooMSS model dynamics.

## Examples

``` r
# Create sample data and plot
env_data <- data.frame(
  time = 1:100,
  dt = 0.01,
  sst = 15 + 3*sin(2*pi*(1:100)/50),
  chl = 0.5 + 0.2*cos(2*pi*(1:100)/50)
)
plots <- plotEnvironment(env_data)
```
