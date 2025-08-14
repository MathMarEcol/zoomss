# Zooplankton Model of Size Spectrum (ZooMSS)

## Overview of ZooMSS

The Zooplankton Model of Size Spectra (ZooMSS) is a functional size-spectrum model of the marine ecosystem (following Heneghan et al. 2016, ZooMSSv1) to resolve phytoplankton, nine zooplankton functional groups (heterotrophic flagellates and ciliates, omnivorous and carnivorous copepods, larvaceans, euphausiids, salps, chaetognaths and jellyfish) and three size-based fish groups. Zooplankton functional groups are resolved using their body-size ranges, size-based feeding characteristics and carbon content, and the zooplankton community emerges from the model across global environmental gradients, depending on the functional traits of the different groups.

We developed the Zooplankton Model of Size Spectra version 2 (ZooMSSv2) based on the prototype of Heneghan et al. (2016). ZooMSSv2 uses the functional size-spectrum framework (Blanchard et al., 2017) to resolve the body size ranges, size-based feeding characteristics and carbon content of nine zooplankton groups and three fish groups. The model supports time-varying environmental conditions enabling studies of seasonal cycles, climate change scenarios, and ecosystem responses to environmental variability.

## Model Parameters: Static vs Dynamic

This implementation distinguishes between static (fixed) and dynamic (time-varying) model components:

### Static Parameters (Fixed)
- Group-specific biological parameters (feeding preferences, growth rates)
- Feeding kernels for predator-prey interactions
- Model grid structure and basic physical constants
- Suitable for components that don't change with environmental conditions

### Dynamic Parameters (Time-varying)
- Sea surface temperature (SST) and chlorophyll concentrations
- Phytoplankton size structure and abundance
- Temperature-dependent biological rates
- Environmental forcing that drives ecosystem dynamics

ZooMSSv2 represents the marine ecosystem as three communities: phytoplankton, zooplankton and fish. The zooplankton community consists of nine of the most abundant zooplankton groups, and the fish community was made up of a small, medium and large group. Dynamics of the phytoplankton are not explicitly resolved in the model, rather the mean size structure of the phytoplankton community is estimated directly from satellite chlorophyll a observations (Brewin et al., 2010; Barnes et al., 2011; Hirata et al., 2011). Abundances of the zooplankton and fish communities are driven by size-dependent processes of growth and mortality, with the temporal dynamics of each functional group governed by separate second-order McKendrick-von Foerster equations.

## How to run ZooMSS using `R`

### Model Structure

The model consists of 5 main files that work together:

1. **setup_RUN_NAME.R** - Set up your specific simulation parameters. This is typically the only file you need to edit.

2. **fZooMSS_Model.R** - The top-level model function which sequentially runs the following 3 functions and saves the output.

3. **fZooMSS_Params.R** - Set all parameter values and process environmental time series for dynamic forcing.

4. **fZooMSS_Setup.R** - Calculate all static components that can be done in advance of the time loop (feeding kernels, initial conditions).

5. **fZooMSS_Run.R** - Contains the main time-dependent simulation loop with dynamic environmental updates.

### Input Files

ZooMSS requires two input files:

1. **TestGroups.csv** - Contains all taxa-specific parameter values, including size ranges and functional group properties.

2. **Environmental data** - A data frame with time series of environmental conditions:
   - For time-varying conditions: Time series with `time_step`, `sst`, and `chlo` columns
   - For constant conditions: Single values replicated across all time steps

### Running the Model

#### Constant Environmental Conditions
For studies with unchanging environmental conditions:
```r
# Create constant environmental conditions
enviro_data <- data.frame(
  time_step = 1:5000,
  sst = rep(15, 5000),  # Constant temperature
  chlo = rep(0.5, 5000) # Constant chlorophyll
)

# Include in your setup parameters
input_params <- cbind(base_params, enviro_data)
source("setup_RUN_NAME.R")
```

#### Time-varying Environmental Conditions
For dynamic environmental studies:
```r
# Create time-varying environmental conditions
enviro_data <- data.frame(
  time_step = 1:5000,
  sst = 15 + 4 * sin(2 * pi * (1:5000) / 365),  # Seasonal cycle
  chlo = 0.5 + 0.3 * cos(2 * pi * (1:5000) / 365)
)

# Include in your setup parameters
input_params <- cbind(base_params, enviro_data)
```

### Environmental Data Structure

Provide environmental time series as a data frame with:
- `time_step`: Sequential time step number (1, 2, 3, ...)
- `sst`: Sea surface temperature (°C)
- `chlo`: Chlorophyll concentration (mg/m³)

The model processes all environmental conditions dynamically, whether they vary through time or remain constant.

### Plotting

The model includes several built-in plotting functions for analysis and visualization:

#### Basic Plots
* `fZooMSS_Plot_BiomassTimeSeries(dat)` - Plot biomass time series with optional species selection and proportional stacking
* `gg <- fZooMSS_Plot_SizeSpectra(dat)` - Plot overall species-resolved size spectra  
* `gg <- fZooMSS_Plot_PPMR(dat)` - Plot predator-prey mass ratios

#### Time Series Analysis (when SaveTimeSteps = TRUE)
* `gg <- fZooMSS_Plot_AbundTimeSeries(dat)` - Plot abundance time series
* `gg <- fZooMSS_Plot_GrowthTimeSeries(dat)` - Plot growth rate time series
* `gg <- fZooMSS_Plot_PredTimeSeries(dat)` - Plot predation time series

where `dat` is the model output list, and `gg` is a ggplot output from the plotting routines.

## Performance and Computational Considerations

### Memory Usage
- Base model: ~10-50 MB depending on grid size
- Environmental time series storage: +5-20 MB 
- Long simulations (>10,000 time steps): Consider chunking environmental data

### Runtime
- Constant environmental conditions: Baseline performance
- Time-varying environmental conditions: ~20-50% longer depending on environmental variability
- Most overhead from phytoplankton parameter recalculation when conditions change

### Recommended Settings
- `dt = 0.01` years (3.65 days) for seasonal resolution
- `dt = 0.001` years (~9 hours) for sub-weekly dynamics  
- Maximum recommended simulation: 1000 years with `dt = 0.01`

### ZooMSS Dashboard
Assuming you have saved all the time-steps (using `SaveTimeSteps <- TRUE`) you can explore the output using our newly created **ZooMSS Dashboard**. To access, go here: https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/. Please drop us a line with any comments or ideas.

## Research Applications

This implementation enables research on:

1. **Ecosystem equilibrium studies** - Constant environmental forcing for baseline ecosystem characterization
2. **Seasonal ecosystem dynamics** - How seasonal environmental cycles affect food web structure  
3. **Climate change impacts** - Long-term ecosystem responses to warming and changing productivity
4. **Environmental variability effects** - Ecosystem responses to different levels of environmental variation
5. **Regional comparisons** - Ecosystem responses across different environmental regimes
6. **Extreme events** - Effects of heat waves, productivity crashes, etc.
7. **Ecosystem indicators** - How environmental variability affects ecosystem metrics

## Validation and Testing

The implementation has been validated by:

1. **Constant conditions comparison**: When environmental time series is constant, results match original model exactly
2. **Conservation tests**: Energy and mass conservation maintained under all environmental conditions  
3. **Numerical stability**: Model remains stable under extreme environmental variations
4. **Performance testing**: Computational overhead is reasonable for typical scenarios

## Publications

1. Heneghan, R.F., Everett, J.D., Blanchard, J.L., Richardson, A.J., 2016. Zooplankton Are Not Fish: Improving Zooplankton Realism in Size-Spectrum Models Mediates Energy Transfer in Food Webs. Front. Mar. Sci. 3, 1–15. https://doi.org/10.3389/fmars.2016.00201

2. Heneghan, R.F., Everett, J.D., Sykes, P., Batten, S.D., Edwards, M., Takahashi, K., Suthers, I.M., Blanchard, J.L., Richardson, A.J., in review, A global size-spectrum model of the marine ecosystem that resolves zooplankton composition. Ecological Modelling

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce simulation length or increase time step
2. **Slow performance**: Check if environmental time series are properly formatted
3. **Unrealistic results**: Validate environmental time series for extreme values
4. **Model differences**: Ensure temperature effects use the same formula (2^((sst-30)/10)) across all conditions

### Performance Tips

1. Use constant environmental conditions for initial testing and parameter exploration
2. Start with shorter simulations to test time-varying scenarios
3. Pre-validate environmental time series for reasonable ranges
4. Monitor memory usage for long simulations
5. Ensure environmental data covers the full simulation period


## Getting Help

If you encounter problems running the model, raise an issue on GitHub: https://github.com/MathMarEcol/ZoopSizeSpectraModel/issues

If you find errors or want to improve the model, we'd love you to make the changes and submit a pull request for us to review and approve.

For questions about the implementation or specific technical issues, please include:
- Your environmental data structure and size
- Simulation parameters (dt, tmax, etc.)
- Error messages or unexpected behavior
- Whether the issue occurs with both constant and time-varying conditions


