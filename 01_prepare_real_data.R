# ST-CoFRK: preparation of EEA air quality observations.
# Converts station measurements to the spatio-temporal BAU
# structure used by the model. Only station-occupied
# spatial BAUs are included during model fitting.

real_data_settings <- list(
  spatial_bau_size_m = 1600,
  temporal_bau_days = 7L,
  
  n_basis_space_side = 10L,
  n_basis_time = 15L,
  
  spatial_range_multiplier = 1.25,
  temporal_scale = NULL,
  
  # Initial parameters for the spatio-temporal precision.
  kappa_space = 0.2,
  rho_time = 0.2,
  gamma_interaction = 1,
  
  # Baseline measurement-error standard deviations (not provided by EEA).
  sigma_noise = c(
    NO2 = 0.00001,
    O3 = 0.00001,
    PM10 = 0.00001,
    PM25 = 0.00001
  ),
  
  # ETRS89 / LAEA Europe: metric coordinate reference system.
  projected_crs = 3035
)


validate_real_dataset <- function(dataset_clean) {
  
  required_columns <- c(
    "lon", "lat", "date", "pollutant", "value",
    "altitude", "population_density"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dataset_clean)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "dataset_clean non contiene le colonne richieste: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


build_week_id <- function(date, temporal_bau_days = 7L) {
  
  date <- as.Date(date)
  
  # Anchor weekly intervals to the preceding Monday.
  min_date <- min(date, na.rm = TRUE)
  weekday <- as.POSIXlt(min_date)$wday
  days_since_monday <- (weekday + 6L) %% 7L
  origin_date <- min_date - days_since_monday
  
  time_id <- as.integer(
    floor(as.numeric(date - origin_date) / temporal_bau_days)
  ) + 1L
  
  list(
    time_id = time_id,
    origin_date = origin_date
  )
}


build_real_spatial_baus <- function(
    dataset_clean,
    spatial_bau_size_m = 1600,
    projected_crs = 3035
) {
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop(
      "Il pacchetto 'sf' e' necessario per costruire le BAU spaziali reali."
    )
  }
  
  locations <- unique(
    dataset_clean[, c("lon", "lat"), drop = FALSE]
  )
  
  locations <- locations[
    is.finite(locations$lon) & is.finite(locations$lat),
    ,
    drop = FALSE
  ]
  
  sf_locations <- sf::st_as_sf(
    locations,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  sf_locations <- sf::st_transform(
    sf_locations,
    crs = projected_crs
  )
  
  projected_coordinates <- sf::st_coordinates(sf_locations)
  
  locations$easting <- projected_coordinates[, 1L]
  locations$northing <- projected_coordinates[, 2L]
  
  # Align the grid origin to multiples of the BAU width.
  origin_easting <- floor(
    min(locations$easting) / spatial_bau_size_m
  ) * spatial_bau_size_m
  
  origin_northing <- floor(
    min(locations$northing) / spatial_bau_size_m
  ) * spatial_bau_size_m
  
  locations$bau_col <- floor(
    (locations$easting - origin_easting) /
      spatial_bau_size_m
  )
  
  locations$bau_row <- floor(
    (locations$northing - origin_northing) /
      spatial_bau_size_m
  )
  
  bau_key <- paste(
    locations$bau_col,
    locations$bau_row,
    sep = "_"
  )
  
  unique_keys <- unique(bau_key)
  
  locations$space_id <- match(
    bau_key,
    unique_keys
  )
  
  locations$bau_easting <- origin_easting +
    (locations$bau_col + 0.5) * spatial_bau_size_m
  
  locations$bau_northing <- origin_northing +
    (locations$bau_row + 0.5) * spatial_bau_size_m
  
  space_grid <- unique(
    locations[
      ,
      c(
        "space_id",
        "bau_col",
        "bau_row",
        "bau_easting",
        "bau_northing"
      ),
      drop = FALSE
    ]
  )
  
  space_grid <- space_grid[
    order(space_grid$space_id),
    ,
    drop = FALSE
  ]
  
  scale_to_unit <- function(x) {
    x_range <- range(x, na.rm = TRUE)
    
    if (diff(x_range) == 0) {
      return(rep(0.5, length(x)))
    }
    
    (x - x_range[1L]) / diff(x_range)
  }
  
  # Scale projected BAU centers to [0, 1] for basis evaluation.
  space_grid$x <- scale_to_unit(space_grid$bau_easting)
  space_grid$y <- scale_to_unit(space_grid$bau_northing)
  
  list(
    locations = locations,
    space_grid = space_grid,
    origin_easting = origin_easting,
    origin_northing = origin_northing
  )
}


prepare_real_data <- function(
    dataset_clean,
    settings = real_data_settings,
    pollutants = c("NO2", "O3", "PM10", "PM25"),
    verbose = TRUE
) {
  
  validate_real_dataset(dataset_clean)
  dataset_clean <- dataset_clean[
    dataset_clean$pollutant %in% pollutants,
    ,
    drop = FALSE
  ]
  
  dataset_clean$date <- as.Date(dataset_clean$date)
  dataset_clean$pollutant <- as.character(dataset_clean$pollutant)
  dataset_clean$value <- as.numeric(dataset_clean$value)
  
  process_lookup <- setNames(
    seq_along(pollutants),
    pollutants
  )
  
  unknown_pollutants <- setdiff(
    unique(dataset_clean$pollutant),
    names(process_lookup)
  )
  
  if (length(unknown_pollutants) > 0L) {
    stop(
      "Pollutant non riconosciuti: ",
      paste(unknown_pollutants, collapse = ", ")
    )
  }
  
  dataset_clean <- dataset_clean[
    is.finite(dataset_clean$lon) &
      is.finite(dataset_clean$lat) &
      is.finite(dataset_clean$value) &
      !is.na(dataset_clean$date),
    ,
    drop = FALSE
  ]
  
  # Assign stations to spatial BAUs.
  
  spatial <- build_real_spatial_baus(
    dataset_clean = dataset_clean,
    spatial_bau_size_m = settings$spatial_bau_size_m,
    projected_crs = settings$projected_crs
  )
  
  location_key <- paste(
    spatial$locations$lon,
    spatial$locations$lat,
    sep = "_"
  )
  
  data_key <- paste(
    dataset_clean$lon,
    dataset_clean$lat,
    sep = "_"
  )
  
  location_match <- match(
    data_key,
    location_key
  )
  
  if (anyNA(location_match)) {
    stop("Errore nell'assegnazione delle osservazioni alle BAU spaziali.")
  }
  
  dataset_clean$space_id <-
    spatial$locations$space_id[location_match]
  
  # Derive BAU-level spatial covariates from station locations.
  
  dataset_clean$altitude <- as.numeric(
    dataset_clean$altitude
  )
  
  dataset_clean$population_density <- as.numeric(
    dataset_clean$population_density
  )
  
  # Clamp negative altitudes before applying log1p.
  dataset_clean$log_altitude <- log1p(
    pmax(dataset_clean$altitude, 0)
  )
  
  # Use unique station locations to avoid weighting repeated observations.
  location_covariates <- unique(
    dataset_clean[
      ,
      c(
        "lon",
        "lat",
        "space_id",
        "log_altitude",
        "population_density"
      )
    ]
  )
  
  bau_covariates <- stats::aggregate(
    cbind(
      log_altitude,
      population_density
    ) ~ space_id,
    data = location_covariates,
    FUN = mean,
    na.rm = TRUE
  )
  
  # Standardize covariates using the station-occupied BAUs.
  logalt_mean <- mean(
    bau_covariates$log_altitude,
    na.rm = TRUE
  )
  
  logalt_sd <- stats::sd(
    bau_covariates$log_altitude,
    na.rm = TRUE
  )
  
  bau_covariates$log_altitude_scaled <-
    (bau_covariates$log_altitude - logalt_mean) /
    logalt_sd
  
  
  pop_mean <- mean(
    bau_covariates$population_density,
    na.rm = TRUE
  )
  
  pop_sd <- stats::sd(
    bau_covariates$population_density,
    na.rm = TRUE
  )
  
  bau_covariates$population_density_scaled <-
    (bau_covariates$population_density - pop_mean) /
    pop_sd
  
  
  spatial$space_grid <- merge(
    spatial$space_grid,
    bau_covariates,
    by = "space_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  spatial$space_grid <- spatial$space_grid[
    order(spatial$space_grid$space_id),
    ,
    drop = FALSE
  ]
  
  
  
  
  
  # Construct weekly temporal BAUs.
  
  week_info <- build_week_id(
    date = dataset_clean$date,
    temporal_bau_days = settings$temporal_bau_days
  )
  
  dataset_clean$time_id <- week_info$time_id
  
  n_time <- max(dataset_clean$time_id)
  
  time_grid <- data.frame(
    time_id = seq_len(n_time)
  )
  
  time_grid$week_start <- week_info$origin_date +
    (time_grid$time_id - 1L) * settings$temporal_bau_days
  
  time_grid$week_end <- time_grid$week_start +
    settings$temporal_bau_days - 1L
  
  time_grid$year <- factor(
    format(time_grid$week_end, "%Y")
  )
  
  time_grid$t_scaled <- if (n_time == 1L) {
    0
  } else {
    seq(0, 1, length.out = n_time)
  }
  
  # Annual and semiannual seasonal harmonics for the mean model.
  time_grid$week_mid <- time_grid$week_start +
    floor((settings$temporal_bau_days - 1L) / 2L)
  
  time_grid$day_of_year <- as.integer(
    format(time_grid$week_mid, "%j")
  )
  
  time_grid$sin_year <- sin(
    2 * pi * time_grid$day_of_year / 365.25
  )
  
  time_grid$cos_year <- cos(
    2 * pi * time_grid$day_of_year / 365.25
  )
  
  time_grid$sin2_year <- sin(
    4 * pi * time_grid$day_of_year / 365.25
  )
  
  time_grid$cos2_year <- cos(
    4 * pi * time_grid$day_of_year / 365.25
  )
  
  # Aggregate observations by pollutant, spatial BAU, and week.
  
  dataset_clean$process_id <- unname(
    process_lookup[dataset_clean$pollutant]
  )
  
  aggregate_formula <- value ~ process_id + pollutant +
    space_id + time_id
  
  weekly_observed <- stats::aggregate(
    aggregate_formula,
    data = dataset_clean,
    FUN = mean,
    na.rm = TRUE
  )
  
  names(weekly_observed)[
    names(weekly_observed) == "value"
  ] <- "Z_obs"
  
  
  weekly_counts <- stats::aggregate(
    value ~ process_id + pollutant + space_id + time_id,
    data = dataset_clean,
    FUN = length
  )
  
  names(weekly_counts)[
    names(weekly_counts) == "value"
  ] <- "n_daily_values"
  
  weekly_observed <- merge(
    weekly_observed,
    weekly_counts,
    by = c(
      "process_id",
      "pollutant",
      "space_id",
      "time_id"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Build the full process-by-space-by-time index on occupied BAUs.
  
  n_process <- length(process_lookup)
  space_grid <- spatial$space_grid
  
  full_index <- expand.grid(
    process_id = seq_len(n_process),
    space_id = space_grid$space_id,
    time_id = time_grid$time_id,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  process_names <- names(process_lookup)
  
  full_index$process <- process_names[
    full_index$process_id
  ]
  
  full_index$pollutant <- full_index$process
  
  full_index$x <- space_grid$x[
    match(full_index$space_id, space_grid$space_id)
  ]
  
  full_index$log_altitude_scaled <-
    space_grid$log_altitude_scaled[
      match(
        full_index$space_id,
        space_grid$space_id
      )
    ]
  
  full_index$population_density_scaled <-
    space_grid$population_density_scaled[
      match(
        full_index$space_id,
        space_grid$space_id
      )
    ]
  
  
  full_index$y <- space_grid$y[
    match(full_index$space_id, space_grid$space_id)
  ]
  
  full_index$t_scaled <- time_grid$t_scaled[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$sin_year <- time_grid$sin_year[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$cos_year <- time_grid$cos_year[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$sin2_year <- time_grid$sin2_year[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$cos2_year <- time_grid$cos2_year[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$year <- time_grid$year[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$week_start <- time_grid$week_start[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$week_end <- time_grid$week_end[
    match(full_index$time_id, time_grid$time_id)
  ]
  
  full_index$st_bau_id <- match(
    paste(
      full_index$space_id,
      full_index$time_id,
      sep = "_"
    ),
    unique(
      paste(
        full_index$space_id,
        full_index$time_id,
        sep = "_"
      )
    )
  )
  
  observed_key <- paste(
    weekly_observed$process_id,
    weekly_observed$space_id,
    weekly_observed$time_id,
    sep = "_"
  )
  
  full_key <- paste(
    full_index$process_id,
    full_index$space_id,
    full_index$time_id,
    sep = "_"
  )
  
  observed_match <- match(
    full_key,
    observed_key
  )
  
  full_index$Z_obs <- weekly_observed$Z_obs[
    observed_match
  ]
  
  full_index$Z_obs_original <- weekly_observed$Z_obs_original[
    observed_match
  ]
  
  full_index$n_daily_values <-
    weekly_observed$n_daily_values[observed_match]
  
  full_index$observed <- !is.na(full_index$Z_obs)
  
  # Compatibility alias for code expecting Y_obs.
  full_index$Y_obs <- full_index$Z_obs
  
  data_full <- full_index
  
  data_obs <- data_full[
    data_full$observed,
    ,
    drop = FALSE
  ]
  
  
  st_grid <- merge(
    space_grid[, c("space_id", "x", "y", "log_altitude_scaled", "population_density_scaled")],
    time_grid[, c(
      "time_id",
      "t_scaled",
      "sin_year",
      "cos_year",
      "sin2_year",
      "cos2_year",
      "year"
    )],
    by = NULL,
    sort = FALSE
  )
  
  st_grid$st_bau_id <- seq_len(nrow(st_grid))
  
  # Match the simulation interface expected by the model functions.
  
  sigma_noise <- settings$sigma_noise[
    process_names
  ]
  
  if (anyNA(sigma_noise) || length(sigma_noise) != n_process) {
    stop(
      "settings$sigma_noise deve contenere un valore per ",
      "NO2, O3, PM10 e PM25."
    )
  }
  
  real_data <- list(
    scenario_name = "real_EEA_weekly",
    covariate_scaling = list(
      pop_mean = pop_mean,
      pop_sd = pop_sd,
      logalt_mean = logalt_mean,
      logalt_sd = logalt_sd
    ),
    config = list(
      scenario_id = "real_EEA_weekly",
      scenario_label = "EEA Northern Italy - weekly BAUs"
    ),
    
    dimensions = list(
      n_process = n_process,
      n_space_side = NA_integer_,
      n_space = nrow(space_grid),
      n_time = n_time
    ),
    
    basis_settings = list(
      n_basis_space_side = settings$n_basis_space_side,
      n_basis_time = settings$n_basis_time
    ),
    
    precision_parameters = list(
      kappa_space = settings$kappa_space,
      rho_time = settings$rho_time,
      gamma_interaction = settings$gamma_interaction
    ),
    
    sigma_noise = as.numeric(sigma_noise),
    
    process_lookup = process_lookup,
    
    spatial_bau_size_m = settings$spatial_bau_size_m,
    temporal_bau_days = settings$temporal_bau_days,
    
    x_values = sort(unique(space_grid$x)),
    y_values = sort(unique(space_grid$y)),
    
    space_grid = space_grid,
    time_grid = time_grid,
    st_grid = st_grid,
    
    data_full = data_full,
    data_obs = data_obs,
    
    weekly_observed = weekly_observed,
    settings = settings
  )
  
  class(real_data) <- c("cofrk_real_data", "list")
  
  if (verbose) {
    cat("\n============================================\n")
    cat("REAL DATA PREPARATI\n")
    cat("============================================\n")
    cat("Processi:", n_process, "\n")
    cat("BAU spaziali occupate:", nrow(space_grid), "\n")
    cat("Dimensione BAU spaziale:",
        settings$spatial_bau_size_m / 1000, "km\n")
    cat("BAU temporali:", n_time, "\n")
    cat("Dimensione BAU temporale:",
        settings$temporal_bau_days, "giorni\n")
    cat("Osservazioni settimanali:", nrow(data_obs), "\n")
    cat("Dimensione data_full:", nrow(data_full), "\n")
    cat(
      "Basi spaziali:",
      settings$n_basis_space_side, "x",
      settings$n_basis_space_side, "=",
      settings$n_basis_space_side^2, "\n"
    )
    cat("Basi temporali Gaussian:",
        settings$n_basis_time, "\n")
  }
  
  real_data
}


prepare_real_data_from_csv <- function(
    path = file.path("data", "dataset_clean.csv"),
    settings = real_data_settings,
    pollutants = c("NO2", "O3", "PM10", "PM25"),
    verbose = TRUE
) {
  
  if (!file.exists(path)) {
    stop("File non trovato: ", path)
  }
  
  dataset_clean <- utils::read.csv(
    path,
    stringsAsFactors = FALSE
  )
  #dataset_clean <- dataset_clean[dataset_clean$year %in% c(2019, 2020, 2021), ]
  prepare_real_data(
    dataset_clean = dataset_clean,
    settings = settings,
    pollutants = pollutants,
    verbose = verbose
  )
}

