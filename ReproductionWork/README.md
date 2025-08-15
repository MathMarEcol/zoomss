# ZooMSS Reproduction Mechanism Implementation

This directory contains all files related to the implementation and testing of the reproduction mechanism in ZooMSS, which replaces the "pinning" of fish abundances with biologically realistic reproduction dynamics.

## Directory Structure

### `/Analysis/`
Contains R scripts for analyzing reproduction mechanism performance and comparing model versions:
- `analyze_reproduction_results.R` - Analysis of reproduction results
- `comprehensive_reproduction_comparison.R` - Comprehensive comparison across environmental gradients
- `corrected_reproduction_visualization.R` - Visualization of corrected reproduction dynamics
- `comprehensive_visualization_analysis.R` - Comprehensive visualization analysis
- `comprehensive_biomass_timeseries.R` - Biomass time series analysis
- `diagnostic_and_timeseries_plots.R` - Diagnostic plotting functions
- `final_results_summary.R` - Final summary of all analyses
- `fixed_visualization_analysis.R` - Fixed visualization analysis
- `improved_visualization_analysis.R` - Improved visualization analysis
- `monitor_simulation_progress.R` - Progress monitoring utilities

### `/TestScripts/`
Contains R scripts for testing the reproduction mechanism:
- Scripts for validating reproduction implementation
- Diagnostic tests for energy allocation
- Performance comparison tests

### `/Documentation/`
Contains documentation files:
- `REPRODUCTION_IMPLEMENTATION.md` - Main implementation documentation
- `ZooMSS_Reproduction_Mechanism_Report.qmd` - Quarto report
- `ZooMSS_Reproduction_Performance_Analysis.md` - Performance analysis documentation

### `/Figures/`
Contains all figure output directories:
- `ZooMSS_Reproduction_Figures/` - Initial reproduction figures
- `ZooMSS_Reproduction_Figures_Corrected/` - Corrected figures
- `ZooMSS_Reproduction_Figures_Fixed/` - Fixed figures
- `ZooMSS_Reproduction_Figures_Improved/` - Improved figures
- `ZooMSS_Biomass_TimeSeries/` - Biomass time series figures

### `/Data/`
Contains data files and backups:
- `reproduction_diagnostics.RData` - Diagnostic data
- `comprehensive_reproduction_comparison_results.RData` - Comparison results
- `TestGroups_Original_Backup.csv` - Backup of original parameters

## Key Changes to Core Model

The reproduction mechanism implementation involved modifications to:
1. `TestGroups.csv` - Added reproduction toggle column
2. `fZooMSS_Params.R` - Added reproduction parameters
3. `fZooMSS_Setup.R` - Energy allocation and reproduction arrays
4. `fZooMSS_Run.R` - Reproduction calculations and conditional pinning

## Energy Allocation

- **Immature fish**: 100% energy to growth
- **Mature fish**: 70% energy to growth, 30% to reproduction
- **Reproduction toggle**: Can be turned on/off per functional group

## Testing Framework

The implementation includes comprehensive testing across:
- Low, medium, and high chlorophyll environments
- 120+ year simulations
- Comparison of old vs new model versions
- Performance diagnostics and validation
