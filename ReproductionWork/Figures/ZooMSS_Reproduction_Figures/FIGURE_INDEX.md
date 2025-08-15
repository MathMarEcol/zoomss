# ZooMSS Reproduction Comparison - Figure Index

This directory contains comprehensive visualizations comparing the old ZooMSS, new ZooMSS (reproduction disabled), and new ZooMSS (reproduction enabled) across low, medium, and high chlorophyll scenarios.

## Figure Categories

### 1. Size Spectrum Analysis
**Purpose**: Compare abundance and biomass distributions across body size ranges

- `size_spectrum_low_abundance.png` - Abundance spectra for low chlorophyll scenario
- `size_spectrum_low_biomass.png` - Biomass spectra for low chlorophyll scenario
- `size_spectrum_medium_abundance.png` - Abundance spectra for medium chlorophyll scenario
- `size_spectrum_medium_biomass.png` - Biomass spectra for medium chlorophyll scenario
- `size_spectrum_high_abundance.png` - Abundance spectra for high chlorophyll scenario
- `size_spectrum_high_biomass.png` - Biomass spectra for high chlorophyll scenario

**Key Insights**: Shows how reproduction affects the shape and magnitude of size distributions across different environmental conditions.

### 2. Growth Rate Comparisons
**Purpose**: Analyze growth rates across size classes and model versions

- `growth_rates_low.png` - Growth rates under low chlorophyll conditions
- `growth_rates_medium.png` - Growth rates under medium chlorophyll conditions
- `growth_rates_high.png` - Growth rates under high chlorophyll conditions

**Key Insights**: Reveals how energy allocation to reproduction affects growth patterns across size classes.

### 3. Reproduction Allocation Analysis
**Purpose**: Examine energy allocation to reproduction across size ranges (reproduction-enabled model only)

- `reproduction_allocation_low.png` - Reproduction allocation under low chlorophyll
- `reproduction_allocation_medium.png` - Reproduction allocation under medium chlorophyll
- `reproduction_allocation_high.png` - Reproduction allocation under high chlorophyll

**Key Insights**: Shows maturity-dependent reproduction patterns and environmental effects on reproductive investment.

### 4. Diet Composition Analysis
**Purpose**: Compare feeding patterns and diet composition between model versions

- `diet_composition_low.png` - Diet composition under low chlorophyll conditions
- `diet_composition_medium.png` - Diet composition under medium chlorophyll conditions
- `diet_composition_high.png` - Diet composition under high chlorophyll conditions

**Key Insights**: Reveals how reproduction affects feeding behavior and trophic interactions.

### 5. Time Series Dynamics
**Purpose**: Show population dynamics over 120-year simulations

- `time_series_low_fish.png` - Fish biomass dynamics under low chlorophyll
- `time_series_low_total.png` - Total biomass dynamics under low chlorophyll
- `time_series_medium_fish.png` - Fish biomass dynamics under medium chlorophyll
- `time_series_medium_total.png` - Total biomass dynamics under medium chlorophyll
- `time_series_high_fish.png` - Fish biomass dynamics under high chlorophyll
- `time_series_high_total.png` - Total biomass dynamics under high chlorophyll

**Key Insights**: Demonstrates stability differences between pinned and reproduction-driven populations.

### 6. Summary Comparisons
**Purpose**: Overall comparison of final model outcomes

- `summary_fish_biomass.png` - Final fish biomass comparison (log scale)
- `summary_total_biomass.png` - Final total biomass comparison (log scale)
- `summary_fish_zoo_ratio.png` - Fish/zooplankton ratio comparison
- `combined_summary.png` - All summary plots combined

**Key Insights**: Quantifies the overall impact of reproduction mechanism on ecosystem structure.

## Key Findings from Visualizations

### Major Differences Revealed:

1. **Environmental Sensitivity**: 
   - Reproduction model shows ~9x higher sensitivity to chlorophyll changes
   - More dramatic population responses to environmental gradients

2. **Size Structure Effects**:
   - Reproduction alters size spectrum shape and slope
   - Different abundance patterns across size classes

3. **Population Dynamics**:
   - Reproduction-driven populations show more variability
   - Different equilibration patterns over 120-year simulations

4. **Energy Allocation Trade-offs**:
   - Clear reproduction allocation patterns based on maturity size
   - Growth rate modifications in mature size classes

5. **Trophic Interactions**:
   - Diet composition changes between model versions
   - Altered predator-prey dynamics

## Technical Details

- **Resolution**: All figures saved at 300 DPI for publication quality
- **Format**: PNG files for easy viewing and sharing
- **Color Scheme**: 
  - Blue: Old ZooMSS
  - Purple: New ZooMSS (no reproduction)
  - Orange: New ZooMSS (with reproduction)
- **Data Source**: 120-year simulations with original parameter values preserved

## Usage Recommendations

1. **Start with**: `combined_summary.png` for overall comparison
2. **Size Analysis**: Review size spectrum plots to understand structural changes
3. **Dynamics**: Examine time series plots for stability assessment
4. **Mechanisms**: Study reproduction allocation and growth rate plots for mechanistic insights
5. **Ecology**: Analyze diet composition plots for trophic implications

## Additional Resources

- `all_plots.RData` - R objects for further customization and analysis
- `../comprehensive_reproduction_comparison_results.RData` - Raw simulation results
- `../model_version_summary.md` - Detailed model documentation

All figures demonstrate the significant ecological implications of implementing reproduction mechanisms in marine ecosystem models, showing both benefits (increased realism) and challenges (increased variability and environmental sensitivity).