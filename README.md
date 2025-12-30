# GEMS Data Processing Pipeline

## Overview

This repository processes GEMS lander data using the `targets` package in R. The pipeline calibrates RGA (residual gas analyzer) data, processes ADV (acoustic Doppler velocimeter) data, calculates eddy covariance fluxes using MATLAB, and produces O₂ and CO₂ flux estimates.

## Pipeline Structure

### Input Data
- Raw GEMS lander files (RGA, ADV, status)
- SeapHOx O₂ calibration data
- ProOceanus CO₂ calibration data
- Lander movement times and bad data periods

### Main Processing Steps

1. **Data Ingestion**: Read raw GEMS lander files and process to clean adv, rga, and status datasets
1. **RGA Processing**: Load, clean, normalize to argon, bin to 15-min intervals, pair high/low measurements
2. **Calibration**: Fit linear models for O₂ (using SeapHOx) and CO₂ (using ProOceanus)
3. **ADV Processing**: Bin velocity data to 15-min intervals, rotate coordinates
4. **Data Integration**: Join RGA, ADV, and temperature data
5. **Flux Calculation**: 
   - Export ADV data to MATLAB format
   - Run eddy covariance analysis in MATLAB
   - Calculate gradient-based O₂ and CO₂ fluxes using ustar and estimated length scale
6. **Output**: Write calibrated flux dataset and hourly statistics

## Running the Pipeline

### Prerequisites
- R with required packages (see `_targets.R`)
- MATLAB

### Execute Pipeline

```r
# In R console
library(targets)
tar_make()
```

### Output Files
- `data/processed/lander/rga_calibrated.parquet` - Calibrated RGA data
- `data/processed/lander/gems_flux.parquet` - Final flux dataset with O₂ and CO₂ fluxes
