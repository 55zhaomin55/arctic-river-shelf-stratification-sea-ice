# Arctic River–Shelf Stratification and Summer Sea-Ice Reproducibility Package

This repository contains the MATLAB code and processed analysis inputs used to reproduce the figures and principal statistical analyses for the manuscript:

**From River Freshwater to Near-Surface Stratification: Shelf-Dependent Links to Arctic Summer Sea Ice**

Authors: **Min Zhao, Xuezhi Bai, Diandian Bai, and Shunqi Pan**

Version: **1.0.0**

## Overview

The analysis examines six Arctic river–shelf systems during 1979–2024:

- Mackenzie–Beaufort
- Yukon–Chukchi
- Kolyma–East Siberian
- Lena–Laptev
- Yenisei–Kara
- Ob–Kara

The repository supports analyses of river discharge, upper-ocean salinity, near-surface haline stratification, summer sea ice, atmospheric variability, spatial-domain sensitivity, and decadal variability.

The code is organized so that figure-generation scripts use standardized processed input files and relative repository paths. Local workstation paths and intermediate working files are not required to reproduce the archived figures.

## Repository Structure

```text
Arctic_River_Shelf_HSI_Reproducibility/
├── README.md
├── LICENSE
├── CITATION.cff
├── VERSION.txt
├── OPEN_RESEARCH_TEMPLATE.md
├── code/
│   ├── project_config.m
│   ├── run_all_figures.m
│   ├── figures/
│   └── functions/
├── data/
│   ├── DATA_MANIFEST.csv
│   ├── README_DATA.md
│   └── geospatial/
├── output/
└── utilities/
    ├── check_required_inputs.m
    └── prepare_repository_inputs_TEMPLATE.m
```

## Main Processed Input Files

The `data/` directory contains standardized analysis inputs used by the public figure scripts. The principal files include:

```text
annual_river_shelf_metrics_1979_2024.mat
arctic_hsi_climatology_1979_2024.mat
marginal_sea_salinity_profiles_1979_2024.mat
runoff_salinity_distance_bins_1979_2024.mat
domain_sensitivity_hsi_sic_1979_2024.mat
sea_ice_monthly_1979_2024.mat
sea_ice_grid_coordinates.mat
wind_speed_JAS_1979_2024.mat
surface_air_temperature_JAS_1979_2024.mat
arctic_oscillation_monthly_1979_2024.mat
discharge_mackenzie_1972_2024.mat
discharge_yukon_1972_2024.mat
discharge_kolyma_1972_2024.mat
discharge_lena_1972_2024.mat
discharge_yenisei_1972_2024.mat
discharge_ob_1972_2024.mat
```

Geospatial inputs are stored under `data/geospatial/`.

See `data/DATA_MANIFEST.csv` and `data/README_DATA.md` for file-level descriptions, variables, source datasets, and provenance information.

## Software Requirements

The scripts are written in MATLAB. The workflow also uses:

- MATLAB Statistics and Machine Learning Toolbox for several statistical functions
- MATLAB Mapping Toolbox for shapefile operations
- M_Map for polar map plotting

The exact MATLAB release used for the archived analysis should be recorded in the Zenodo/GitHub release notes before final publication.

## Reproducing the Figures

From the repository root, first check that all required standardized inputs are present:

```matlab
run('utilities/check_required_inputs.m')
```

Then add the code directory and run all figure scripts:

```matlab
addpath('code')
run('code/run_all_figures.m')
```

Generated figures are written to:

```text
output/
```

Individual figures can also be reproduced by running the corresponding scripts under:

```text
code/figures/
```

## Data Preparation

The public repository contains only a template for converting local working files into the standardized repository inputs:

```text
utilities/prepare_repository_inputs_TEMPLATE.m
```

The authors' workstation-specific preparation script is not distributed because it contains local absolute file paths. Users who wish to construct equivalent inputs from the original public data products should adapt the template using the source information documented in `data/README_DATA.md` and `data/DATA_MANIFEST.csv`.

## Data Sources

The analyses use data derived from the following public data products:

- NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice Concentration, Version 5
- ORAS5 global ocean reanalysis
- ERA5 atmospheric reanalysis
- Arctic Great Rivers Observatory discharge records
- NOAA Climate Prediction Center Arctic Oscillation index

The manuscript provides the formal dataset citations and access information.

Original third-party datasets are not relicensed by this repository. Processed files supplied for reproducibility retain the attribution and usage requirements of their original data providers.

## Reproducibility Notes

The archived code is intended to reproduce the analysis and figure-generation workflow used for the submitted manuscript. The main statistical procedures include:

- calculation of the Haline Stratification Index in the upper ocean
- Pearson correlation analysis using detrended annual series
- higher-order partial correlation controlling simultaneously for atmospheric variables
- spatial-domain sensitivity analysis
- lagged runoff–salinity analysis
- 10-year moving correlations
- binned HSI–sea-ice concentration analysis

Numerical results should be reproduced from the archived processed inputs rather than transcribed manually from the manuscript.

## Citation

If you use this code or processed dataset, please cite the archived release. The final Zenodo DOI should be inserted below after deposition:

> Zhao, M., Bai, X., Bai, D., & Pan, S. (2026). *Reproducibility code and processed data for “From River Freshwater to Near-Surface Stratification: Shelf-Dependent Links to Arctic Summer Sea Ice”* (Version 1.0.0) [Software and data set]. Zenodo. **DOI TO BE ADDED**

Citation metadata are also provided in `CITATION.cff`.

## License

The MATLAB source code in this repository is released under the **MIT License**. See the `LICENSE` file for details.

The MIT License applies to the source code only. Data files derived from third-party products remain subject to the attribution, licensing, and usage terms of their original data providers.

## Contact

**Min Zhao**  
College of Oceanography, Hohai University, Nanjing, China  
School of Engineering, Cardiff University, Cardiff, United Kingdom

**Corresponding author:**  
Xuezhi Bai  
College of Oceanography, Hohai University  
Email: xuezhi.bai@hhu.edu.cn

## Version History

### Version 1.0.0
Initial archived reproducibility release accompanying the manuscript submission.
