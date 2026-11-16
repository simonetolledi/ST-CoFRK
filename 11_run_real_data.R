
############################################################
# 11_run_real_data.R
#
# ST-CoFRK application to EEA air quality data
# for Northern Italy.
# Fits the multivariate model and generates spatial
# predictions, diagnostic plots, and municipal time series.
############################################################

source("01_prepare_real_data.R")
source("02_generate_multires_st_basis.R")
source("03_build_precision_matrices.R")
source("04_build_model_object.R")
source("05_em_fit.R")
source("06_update_precision_parameters.R")


# Model settings and spatial-temporal basis dimensions.

settings <- real_data_settings

settings$n_basis_space_side <- 10L
settings$n_basis_time <- 20L

# Additional settings can be modified in real_data_settings.


pollutants_to_use = c("NO2", "O3", "PM10", "PM25")

real_data <- prepare_real_data_from_csv(
  settings = settings,
  pollutants = pollutants_to_use
)


basis <- generate_multires_st_basis(
  simulation = real_data,
  spatial_range_multiplier =
    settings$spatial_range_multiplier,
  temporal_scale =
    settings$temporal_scale,
  verbose = TRUE
)


precision <- build_precision_matrices(
  basis = basis,
  kappa_space =
    real_data$precision_parameters$kappa_space,
  rho_time =
    real_data$precision_parameters$rho_time,
  gamma_interaction =
    real_data$precision_parameters$gamma_interaction,
  verbose = TRUE
)


# Process-specific fixed effects: seasonal harmonics,
# population density, and altitude.

mean_formula <- ~ sin_year + cos_year + sin2_year + cos2_year + population_density_scaled + log_altitude_scaled

model <- build_model_object(
  simulation = real_data,
  basis = basis,
  mean_formula = mean_formula,
  verbose = TRUE
)


fit <- em_fit(
  model = model,
  precision = precision,
  basis = basis,
  n_em = 20,
  tolerance = 1e-3,
  initial_sigma_fine = c(1, 1, 1, 1),
  update_beta = TRUE,
  update_sigma_fine = TRUE,
  update_precision = TRUE,
  verbose = TRUE
)


# Predict on selected weeks over the Northern Italy BAU grid.

source("07_predict_northern_italy_grid.R")
predictions <- predict_cofrk_northern_italy(
  real_data = real_data,
  fit = fit,
  basis = basis,
  time_ids = c(55, 72, 81),#c(3, 29, 103, 133, 160, 238),
  verbose = TRUE
)


source("09_real_data_plots.R")

plots <- make_real_data_plots(
  real_data = real_data,
  model = model,
  fit = fit,
  predictions = predictions,
  time_ids = c(55, 72, 81)#c(3, 20, 55, 72, 107, 124)
)

plots$predicted_maps
plots$observed_maps
plots$observed_vs_fitted
plots$temporal_response
plots$process_correlation


# Municipality-level predictions across the study period.

source("10_predict_municipality_timeseries.R")

comune_ts <- predict_cofrk_municipality_timeseries(
  real_data = real_data,
  fit = fit,
  basis = basis,
  comune = "Milano",
  comuni_shp = "Com01012026_g/Com01012026_g_WGS84.shp",
  mean_formula = mean_formula
)

plot_municipality_timeseries(comune_ts, real_data)



