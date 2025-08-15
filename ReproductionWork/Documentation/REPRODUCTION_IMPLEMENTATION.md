# ZooMSS Reproduction Mechanism Implementation

## Overview

This implementation adds a reproduction mechanism to the ZooMSS marine ecosystem model, allowing fish groups to reproduce instead of relying on the traditional "pinning" mechanism where fish abundance is set equal to zooplankton abundance at the smallest size classes.

## Key Changes Made

### 1. Group Parameters (`TestGroups.csv`, `TestGroups_WithFishing.csv`)
- Added `Reproduction` column (0 = traditional pinning, 1 = reproduction mechanism)
- Fish groups now have `Reproduction = 1`
- Zooplankton groups maintain `Reproduction = 0`

### 2. Parameters (`fZooMSS_Params.R`)
- Added `reproduction_on` parameter that checks if any group has reproduction enabled
- This parameter controls the overall reproduction system state

### 3. Setup Function (`fZooMSS_Setup.R`)
- **Energy Allocation**: Split net growth efficiency between growth and reproduction
  - Immature individuals: 100% energy to growth
  - Mature individuals: 70% energy to growth, 30% to reproduction
  - Smooth transition using sigmoid function based on maturity size (`Wmat`)
  
- **Initial Conditions**: 
  - Without reproduction: Fish abundance pinned to zooplankton levels
  - With reproduction: Fish start at 10% of pinning levels to allow population persistence
  
- **Storage**: Added arrays for reproduction rates and recruitment tracking

### 4. Run Function (`fZooMSS_Run.R`)
- **Reproduction Calculations**:
  - Calculate energy available for reproduction from feeding
  - Apply reproduction efficiency to determine reproduction rates
  - Calculate recruitment to smallest size classes following DBPM approach
  
- **Conditional Fish Pinning**: 
  - Fish pinning only occurs when `reproduction_on = FALSE`
  - When reproduction is enabled, fish abundance is determined by growth, mortality, and recruitment
  
- **Output Storage**: Save reproduction rates and recruitment for analysis

## Usage

### Basic Usage
1. Set `Reproduction = 1` for fish groups in the Groups CSV file
2. Set `Reproduction = 0` for zooplankton groups  
3. Run model normally - reproduction will be automatically enabled

### Advanced Configuration
- Adjust maturity sizes (`Wmat`) to control when reproduction begins
- Modify reproduction energy fraction (currently 30%) in Setup function
- Adjust initial fish abundance scaling factor (currently 10%) in Setup function

## Model Behavior

### Without Reproduction (Traditional)
- Fish abundance at smallest size classes is pinned to zooplankton abundance
- Energy allocation: 100% to growth
- Fish populations maintained artificially

### With Reproduction  
- Fish abundance determined by population dynamics
- Energy split between growth (70%) and reproduction (30%) for mature individuals
- Recruitment from reproductive output maintains populations
- More realistic population dynamics

## Biological Realism

The implementation follows the DBPM (Dynamic Bioclimate Envelope Model) approach:

1. **Energy Budget**: Divides net growth efficiency between somatic growth and reproduction
2. **Maturity**: Reproduction only occurs after reaching maturity size
3. **Recruitment**: New individuals added to smallest size class based on reproductive output
4. **Size-Based**: Reproduction scales with body size and abundance

## Technical Notes

### Energy Allocation Formula
```r
maturity_function <- 1 / (1 + exp(-5 * (log10(w) - maturity_size)))
reproduction_fraction <- 0.3 * maturity_function
growth_fraction <- 1 - reproduction_fraction
```

### Recruitment Calculation
```r
reproduction_biomass <- sum(reproduction_rates * body_mass * abundance * dx)
recruitment_rate <- reproduction_biomass / (dx * min_body_mass)
```

## Validation Results

Testing shows:
- ✓ Model runs successfully with both reproduction on/off
- ✓ Energy allocation works correctly across size classes
- ✓ Fish pinning is disabled when reproduction is enabled
- ✓ Backward compatibility maintained
- ✓ Reproduction produces self-sustaining fish populations

## Future Enhancements

1. **Dynamic Feeding**: More sophisticated calculation of feeding rates for reproduction
2. **Species-Specific Parameters**: Different reproduction strategies per species
3. **Environmental Effects**: Temperature or food-dependent reproduction rates
4. **Spawning Seasonality**: Temporal variation in reproduction
5. **Density Dependence**: Reproduction rates affected by population density

## Files Modified

1. `TestGroups.csv` - Added Reproduction column
2. `TestGroups_WithFishing.csv` - Added Reproduction column  
3. `fZooMSS_Params.R` - Added reproduction_on parameter
4. `fZooMSS_Setup.R` - Energy allocation and initial conditions
5. `fZooMSS_Run.R` - Reproduction calculations and recruitment

## Test Scripts

- `test_reproduction.R` - Basic functionality test
- `diagnostic_reproduction.R` - Detailed analysis of reproduction dynamics
- `test_improved_reproduction.R` - Parameter optimization test
- `final_reproduction_test.R` - Comprehensive validation

The reproduction mechanism is now fully integrated and ready for use in ZooMSS simulations.
