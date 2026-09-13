# Data files required by the figure code

Place the processed data files listed in `DATA_MANIFEST.csv` in this directory. Geospatial files belong in `data/geospatial/`.

The public repository should contain the **processed/refined inputs actually used to generate the figures**. The large third-party raw products (ORAS5, ERA5, NSIDC CDR, ArcticGRO) do not need to be duplicated here if they are already preserved by their official repositories; cite those source datasets in the manuscript and repository metadata.

## Recommended file renaming from the working directory

- `HSI.mat` -> `annual_river_shelf_metrics_1979_2024.mat`
- `ORAS5_HSI_JAS_climatology_1979_2024.mat` -> `arctic_hsi_climatology_1979_2024.mat`
- `ORAS5_marginal_sea_monthly_profiles_1979_2024.mat` -> `marginal_sea_salinity_profiles_1979_2024.mat`
- `ORAS5_S15_distance_bins_1979_2024.mat` -> `runoff_salinity_distance_bins_1979_2024.mat`
- `Figure6_HSI_SIC_domain_sensitivity_final.mat` -> `domain_sensitivity_hsi_sic_1979_2024.mat`
- `sic.mat` -> `sea_ice_monthly_1979_2024.mat`
- `lon_lat.mat` -> `sea_ice_grid_coordinates.mat`
- `wind.mat` -> `wind_speed_JAS_1979_2024.mat`
- `SAT.mat` -> `surface_air_temperature_JAS_1979_2024.mat`
- `AO_index.mat` -> `arctic_oscillation_monthly_1979_2024.mat`

The six monthly discharge files should be standardized to the river-specific filenames listed in `DATA_MANIFEST.csv`; the recommended public variable name is `discharge_m3s` (the loader also accepts the legacy `discharge8`). The six drainage-basin boundary MAT files should likewise be renamed and placed under `data/geospatial/`.

## Important reproducibility checks

1. `450 km / 0°` values in Figure 6 should match the final HSI-SIC correlations used in Figure 8.
2. The Mackenzie-Beaufort partial correlation in Figure S2 should reproduce the Figure 10/Table S2 value (approximately 0.75 in the current analysis).
3. Figure 9, Figure 10, Figure 11, and Figure S3 must all use the same final HSI time series contained in `annual_river_shelf_metrics_1979_2024.mat`.
4. Keep the sampling-center coordinates centralized in `code/functions/river_shelf_systems.m`.
