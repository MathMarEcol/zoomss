# ZooMSS Biomass Time Series Analysis - Figure Index

## Overview
This directory contains comprehensive biomass time series plots from the ZooMSS reproduction mechanism analysis. All plots compare three model versions across three environmental scenarios over 120 years.

**Model Versions:**
- **Old ZooMSS**: Traditional pinning mechanism
- **New ZooMSS (No Repro)**: New framework with reproduction disabled
- **New ZooMSS (With Repro)**: New framework with reproduction enabled

**Environmental Scenarios:**
- **Low**: Oligotrophic conditions (low chlorophyll-a)
- **Medium**: Mesotrophic conditions (medium chlorophyll-a)  
- **High**: Eutrophic conditions (high chlorophyll-a)

---

## Community-Level Time Series

### Total Community Biomass
- **[`low_community.png`](low_community.png)** - Total community biomass over time in low chlorophyll scenario
- **[`medium_community.png`](medium_community.png)** - Total community biomass over time in medium chlorophyll scenario
- **[`high_community.png`](high_community.png)** - Total community biomass over time in high chlorophyll scenario

**Key Insights:**
- Shows dramatic differences between reproduction on/off
- Medium scenario shows strongest reproduction benefits
- Low scenario shows reproduction challenges

### Fish vs Zooplankton Biomass
- **[`low_fish_zoo.png`](low_fish_zoo.png)** - Fish and zooplankton biomass time series in low chlorophyll scenario
- **[`medium_fish_zoo.png`](medium_fish_zoo.png)** - Fish and zooplankton biomass time series in medium chlorophyll scenario
- **[`high_fish_zoo.png`](high_fish_zoo.png)** - Fish and zooplankton biomass time series in high chlorophyll scenario

**Key Insights:**
- Reveals trophic coupling between fish and zooplankton
- Shows environmental dependency of reproduction effects
- Demonstrates ecosystem-level responses

### Fish:Zooplankton Ratios
- **[`low_ratio.png`](low_ratio.png)** - Fish:zooplankton biomass ratio over time in low chlorophyll scenario
- **[`medium_ratio.png`](medium_ratio.png)** - Fish:zooplankton biomass ratio over time in medium chlorophyll scenario
- **[`high_ratio.png`](high_ratio.png)** - Fish:zooplankton biomass ratio over time in high chlorophyll scenario

**Key Insights:**
- Quantifies ecosystem balance changes
- Shows reproduction mechanism creates more balanced ecosystems
- Reveals productivity-dependent trophic structure

---

## Functional Group Time Series

### Fish Groups
- **[`low_fish_groups.png`](low_fish_groups.png)** - Individual fish group biomass time series in low chlorophyll scenario
- **[`medium_fish_groups.png`](medium_fish_groups.png)** - Individual fish group biomass time series in medium chlorophyll scenario
- **[`high_fish_groups.png`](high_fish_groups.png)** - Individual fish group biomass time series in high chlorophyll scenario

**Key Insights:**
- Shows size-specific responses to reproduction mechanism
- Reveals differential impacts across fish functional groups
- Demonstrates model version differences for each fish group

### Zooplankton Groups
- **[`low_zoo_groups.png`](low_zoo_groups.png)** - Individual zooplankton group biomass time series in low chlorophyll scenario
- **[`medium_zoo_groups.png`](medium_zoo_groups.png)** - Individual zooplankton group biomass time series in medium chlorophyll scenario
- **[`high_zoo_groups.png`](high_zoo_groups.png)** - Individual zooplankton group biomass time series in high chlorophyll scenario

**Key Insights:**
- Shows indirect effects of fish reproduction on zooplankton
- Reveals functional group-specific responses
- Demonstrates trophic cascade effects

### All Functional Groups Overview
- **[`low_all_groups.png`](low_all_groups.png)** - All functional groups biomass time series in low chlorophyll scenario
- **[`medium_all_groups.png`](medium_all_groups.png)** - All functional groups biomass time series in medium chlorophyll scenario
- **[`high_all_groups.png`](high_all_groups.png)** - All functional groups biomass time series in high chlorophyll scenario

**Key Insights:**
- Comprehensive view of entire ecosystem response
- Shows relative importance of different functional groups
- Reveals ecosystem-wide impacts of reproduction mechanism

---

## Equilibrium Analysis

### Community Equilibrium Comparison
- **[`equilibrium_community.png`](equilibrium_community.png)** - Final equilibrium total biomass comparison across all scenarios and model versions

**Key Insights:**
- Quantifies final state differences between model versions
- Shows environmental gradient effects on reproduction benefits
- Demonstrates productivity thresholds for reproduction success

### Fish vs Zooplankton Equilibrium
- **[`equilibrium_fish_zoo.png`](equilibrium_fish_zoo.png)** - Final equilibrium fish and zooplankton biomass comparison across scenarios

**Key Insights:**
- Shows differential impacts on fish vs zooplankton
- Reveals trophic structure changes at equilibrium
- Quantifies ecosystem balance shifts

---

## Data Files

### Analysis Results
- **[`biomass_timeseries_analysis.RData`](biomass_timeseries_analysis.RData)** - Complete R workspace with all plots and processed time series data

**Contents:**
- `all_plots`: List of all ggplot objects
- `all_timeseries_data`: Complete time series dataset
- Summary statistics and analysis results

---

## Key Findings Summary

### Environmental Dependency
1. **High Productivity**: +14% fish biomass, +102% total biomass
2. **Medium Productivity**: +209% fish biomass, +281% total biomass (strongest response)
3. **Low Productivity**: -86% fish biomass, -35% total biomass (reproduction struggles)

### Ecosystem Structure Changes
- More balanced fish:zooplankton ratios with reproduction
- Enhanced trophic coupling in productive environments
- Realistic population dynamics replace artificial pinning

### Model Performance
- Full backward compatibility maintained
- 120-year simulations stable across all scenarios
- Biologically realistic reproduction dynamics implemented

---

## Technical Details

**Simulation Parameters:**
- Duration: 120 years
- Time step: 0.01 years
- Model versions: 3
- Environmental scenarios: 3
- Total simulations: 9

**Plot Specifications:**
- Resolution: 300 DPI
- Format: PNG
- Community plots: 12" × 8"
- Functional group plots: 14" × 10"
- Overview plots: 16" × 12"

**Analysis Date:** August 15, 2025  
**Generated by:** `comprehensive_biomass_timeseries.R`  
**Source Data:** `comprehensive_reproduction_comparison_results.RData`

---

## Usage Notes

These plots provide comprehensive documentation of the ZooMSS reproduction mechanism performance across environmental gradients. They demonstrate:

1. **Biological Realism**: Reproduction mechanism creates realistic population dynamics
2. **Environmental Responsiveness**: Strong dependency on ecosystem productivity
3. **Ecosystem Effects**: Community-wide impacts beyond just fish populations
4. **Model Validation**: Successful implementation with maintained backward compatibility

For detailed analysis and interpretation, see [`ZooMSS_Reproduction_Performance_Analysis.md`](../ZooMSS_Reproduction_Performance_Analysis.md).