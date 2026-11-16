# ST-CoFRK: A Multivariate Spatio-Temporal Extension of Fixed Rank Kriging

This repository contains the R implementation of **ST-CoFRK**, developed as part of a Master's thesis in Mathematical Engineering at Politecnico di Milano.

ST-CoFRK extends Fixed Rank Kriging (FRK) to jointly model multiple correlated spatio-temporal processes. The model combines tensor-product spatial and temporal basis functions, a sparse precision matrix accounting for spatial and temporal dependence, and a cross-process covariance matrix. Parameter estimation is performed through an EM-type algorithm.

The repository includes the model implementation, simulation experiments comparing ST-CoFRK with univariate FRK, a parameter recovery study, and an application to air pollution data in Northern Italy.

## File description

### Model implementation

These scripts contain the main functions shared by the simulation and real-data applications.

| File | Description |
|---|---|
| `00_load_all.R` | Checks the required R packages and loads the modules used in the simulation pipeline. |
| `02_generate_multires_st_basis.R` | Constructs spatial bisquare and temporal Gaussian basis functions and their spatio-temporal representation. |
| `03_build_precision_matrices.R` | Constructs the sparse spatial, temporal, and spatio-temporal precision matrices. |
| `04_build_model_object.R` | Builds the model object, including the fixed-effect and random-effect design matrices. |
| `05_em_fit.R` | Implements the EM-type estimation algorithm for model parameters and latent coefficients. |
| `06_update_precision_parameters.R` | Estimates the spatial, temporal, and interaction precision parameters using numerical optimization and stochastic approximations. |
| `07_predict.R` | Computes ST-CoFRK predictions using the fitted model. |

### Simulation experiments

These scripts define the simulation settings and evaluate the predictive performance of ST-CoFRK.

| File | Description |
|---|---|
| `01s_simulation_scenarios.R` | Defines simulation scenarios, generates synthetic spatio-temporal fields, and constructs observation patterns with different sampling conditions and missing-data regions. |
| `08s_fit_FRK_benchmark.R` | Fits univariate spatio-temporal FRK models for comparison with ST-CoFRK. |
| `09s_metrics_and_plots.R` | Computes prediction errors and generates graphical summaries of the simulation results. |
| `10s_compare_models.R` | Summarizes model estimates and predictive performance across simulation runs. |
| `11s_run_example.R` | Main script for running Monte Carlo simulations, comparing ST-CoFRK with FRK, and saving the results. |

### Parameter recovery

| File | Description |
|---|---|
| `1rec_parameter_recovery.R` | Generates data directly from the ST-CoFRK model with known parameters and evaluates parameter recovery across Monte Carlo repetitions. |

### Real-data application

These scripts implement the application to EEA air quality observations in Northern Italy, jointly modeling NO2, O3, PM10, and PM2.5.

| File | Description |
|---|---|
| `01_prepare_real_data.R` | Prepares environmental observations, constructs spatial and temporal BAUs, aggregates measurements, and defines model covariates. |
| `07_predict_northern_italy_grid.R` | Generates predictions over a spatial grid covering Northern Italy using the fitted ST-CoFRK model. |
| `09_real_data_plots.R` | Produces maps, temporal visualizations, observed-versus-fitted comparisons, and cross-process correlation plots. |
| `10_predict_municipality_timeseries.R` | Computes municipality-level predictions over time using the fitted model. |
| `11_run_real_data.R` | Main script for fitting ST-CoFRK to the environmental dataset and producing spatial and temporal predictions. |

## Running the code

The repository is designed to be used through **three main scripts**, which allow users to run the complete experiments without modifying the underlying model implementation.

All scripts should be executed from the project directory in R or RStudio.

### Simulation experiments

The file `11s_run_example.R` runs the Monte Carlo experiments, including ST-CoFRK estimation, comparison with univariate FRK, and generation of predictive metrics and plots.

Users can select the simulation scenario and number of Monte Carlo repetitions directly in this script. Scenario configurations and parameters can be customized in `01s_simulation_scenarios.R`.

```r
source("11s_run_example.R")
```

Results are saved in `results_scenarios/`.

### Parameter recovery

The file `1rec_parameter_recovery.R` evaluates parameter estimation using synthetic data generated directly from the ST-CoFRK model.

Users can modify the true parameter values, simulation dimensions, and number of Monte Carlo repetitions at the beginning of the script.

```r
source("1rec_parameter_recovery.R")
```

Results are saved in `results_parameter_recovery/`.

### Real-data application

The file `11_run_real_data.R` runs the complete analysis of air pollution data in Northern Italy, including model fitting, spatial predictions, and visualizations.

Users can select the pollutants, basis dimensions, and estimation settings directly in this script. Additional preprocessing parameters are defined in `01_prepare_real_data.R`.

```r
source("11_run_real_data.R")
```

The application requires the environmental dataset and geographical files used for spatial prediction and municipality-level analyses.

All remaining scripts contain supporting functions and are loaded by the corresponding main scripts.

## Data Requirements

The simulation experiments (`11s_run_example.R`) and parameter recovery study (`1rec_parameter_recovery.R`) are fully synthetic and **do not require external datasets**.

The real-data application (`11_run_real_data.R`) requires environmental observations and additional geographical datasets. These files are not distributed with the repository and must be obtained and placed in the appropriate directories before running the analysis.

### 1. Air quality observations — EEA

**Source:** [European Environment Agency — Air Quality Download Service](https://eeadmz1-downloads-webapp.azurewebsites.net/)

The application uses air quality observations from Northern Italy for four pollutants: NO2, O3, PM10, and PM2.5.

The script `01_prepare_real_data.R` expects a preprocessed CSV file located at:

`data/dataset_clean.csv`

The CSV must contain the following columns:

| Column | Description |
|---|---|
| `lon` | Monitoring station longitude in decimal degrees (WGS84) |
| `lat` | Monitoring station latitude in decimal degrees (WGS84) |
| `date` | Observation date, preferably in `YYYY-MM-DD` format |
| `pollutant` | Pollutant identifier: `NO2`, `O3`, `PM10`, or `PM25` |
| `value` | Observed pollutant concentration, in consistent units |
| `altitude` | Altitude associated with the monitoring location, in metres |
| `population_density` | Population covariate associated with the monitoring location |

The input file must already contain the spatial covariates `altitude` and `population_density`. These are not automatically retrieved for monitoring stations by the preprocessing script.

The preparation routine:

1. Filters the selected pollutants and invalid observations.
2. Constructs spatial Basic Areal Units (BAUs) of 1.6 km × 1.6 km using the EPSG:3035 coordinate reference system.
3. Constructs temporal BAUs of seven days, aligned to Mondays.
4. Aggregates observations by pollutant, spatial BAU, and week using their arithmetic mean.
5. Aggregates station-level spatial covariates to the BAU level and standardizes them.
6. Constructs annual seasonal covariates and the observation and prediction structures required by ST-CoFRK.

The input must use consistent concentration units and pollutant identifiers. The repository does not include a complete preprocessing pipeline from raw EEA downloads to `dataset_clean.csv`.

### 2. Population data — Eurostat Census Grid 2021

**Source:** [Eurostat — Population Grids](https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/population-grids)

The spatial prediction routine uses the 2021 European population grid at 1 km resolution.

The code expects the following GeoPackage:

`eurostat_tmp/Eurostat_Census-GRID_2021_V3/ESTAT_Census_2021_V3.gpkg`

The GeoPackage must contain:

- A layer named `ESTAT_Census_2021_V3`.
- A numeric attribute `T`, representing total population in the grid cell.
- Valid spatial geometries compatible with the prediction grid.

The prediction routine identifies the Eurostat cell intersecting each prediction location, extracts the population value, and standardizes it using the scaling parameters computed during model fitting.

The expected GeoPackage path, layer name, and field name are hard-coded in `07_predict_northern_italy_grid.R`. Users downloading a different Eurostat release may need to adjust them.

### 3. Altitude data — Copernicus DEM

**Source:** [Copernicus Digital Elevation Model](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)

The application uses Copernicus DEM GLO-90 to obtain altitude values across the prediction domain.

The script expects a GeoTIFF file at:

`data/copernicus_dem_90_northern_italy.tiff`

This must be a georeferenced raster covering the Northern Italy prediction domain. If the DEM is downloaded as multiple tiles, they must first be combined or otherwise prepared as a single compatible raster at the expected location.

The routine extracts elevation values at prediction locations, applies the transformation `log1p(pmax(altitude, 0))`, and standardizes them using the scaling parameters from model fitting.

### 4. Administrative boundaries — ISTAT

**Source:** [ISTAT — Administrative Boundaries](https://www.istat.it/notizia/confini-delle-unita-amministrative-a-fini-statistici-al-1-gennaio-2018-2/)

Municipality-level predictions and some spatial visualizations require the administrative boundaries of Italian municipalities.

The code expects the shapefile:

`Com01012026_g/Com01012026_g_WGS84.shp`

The dataset must contain the following attributes:

| Attribute | Description |
|---|---|
| `PRO_COM_T` | Municipality identification code |
| `COMUNE` | Municipality name |
| `COD_REG` | Italian regional code |

The municipality prediction routine uses these boundaries to identify the spatial BAUs intersecting the selected municipality and compute area-weighted predictions.

The implementation uses municipality boundaries corresponding to the 2026 ISTAT release and selects Northern Italian regions using regional codes 1–8.

**Important:** The expected shapefile name reflects the locally prepared version used in this project. The original ISTAT download may require extraction, coordinate reference system conversion, or filename adjustments. All associated shapefile components (`.shp`, `.shx`, `.dbf`, `.prj`) must be kept together.

Northern Italy's outer boundary is obtained separately through `giscoR` using the Eurostat NUTS-1 regions `ITC` and `ITH`. This boundary does not require a manually supplied shapefile, but it may require an internet connection when not cached.

### 5. Expected directory structure

All paths are relative to the repository's root directory.

```text
ST-CoFRK/
│
├── 00_load_all.R
├── 01_prepare_real_data.R
├── 01s_simulation_scenarios.R
├── ...
├── 11_run_real_data.R
├── 11s_run_example.R
├── 1rec_parameter_recovery.R
│
├── data/
│   ├── dataset_clean.csv
│   └── copernicus_dem_90_northern_italy.tiff
│
├── eurostat_tmp/
│   └── Eurostat_Census-GRID_2021_V3/
│       └── ESTAT_Census_2021_V3.gpkg
│
└── Com01012026_g/
    ├── Com01012026_g_WGS84.shp
    ├── Com01012026_g_WGS84.shx
    ├── Com01012026_g_WGS84.dbf
    └── Com01012026_g_WGS84.prj
```

### 6. R package dependencies

The main model and simulation scripts use:

```r
install.packages(c(
  "MASS", "Matrix", "FRK", "ggplot2",
  "sp", "spacetime", "zoo"
))
```

The real-data application additionally requires:

```r
install.packages(c(
  "sf", "terra", "giscoR", "dplyr", "patchwork"
))
```

The `splines` and `stats` packages are included with standard R installations.

After installing the required packages and preparing the external datasets, open the project in RStudio and execute `11_run_real_data.R` from the project directory.

**Reproducibility note:** The repository provides the model implementation and the routines for processing the expected input data. Exact reproduction of the original real-data application requires the corresponding preprocessed EEA observations, geographic covariates, and compatible source-data versions.

## Research context

This code was developed for the Master's thesis *A Multivariate Spatio-Temporal Extension of Fixed Rank Kriging*, Politecnico di Milano, 2026.