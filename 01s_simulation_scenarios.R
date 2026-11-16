############################################################
# 01s_simulation_scenarios.R
#
# Defines synthetic multivariate spatio-temporal scenarios.
# simulate_scenario() returns the scenario configuration,
# spatial and temporal grids, and full and observed data.
############################################################

get_simulation_scenarios <- function() {
  list(
    S0_independent = list(
      scenario_id = "S0_independent",
      scenario_label = "Independent processes: no borrowing expected",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.0,
      rho_time = 0.2,
      gamma_interaction = 1.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.10, 0.10),
      create_hole_process2 = FALSE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 25,
      hole_t_max = 35,
      use_event_1 = FALSE,
      use_event_2 = FALSE,
      observation_pattern = "random"
    ),
    
    S1_weak_correlation = list(
      scenario_id = "S1_weak_correlation",
      scenario_label = "Weak cross-process correlation",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.4,
      rho_time = 0.2,
      gamma_interaction = 1.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.50, 0.20),
      create_hole_process2 = TRUE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 25,
      hole_t_max = 35,
      use_event_1 = TRUE,
      use_event_2 = TRUE,
      observation_pattern = "random"
    ),
    
    S2_strong_correlation = list(
      scenario_id = "S2_strong_correlation",
      scenario_label = "Strong cross-process correlation, no ST interaction",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.8,
      rho_time = 0.2,
      gamma_interaction = 0.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.20, 0.20),
      create_hole_process2 = FALSE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 25,
      hole_t_max = 35,
      use_event_1 = FALSE,
      use_event_2 = TRUE,
      observation_pattern = "random"
    ),
    
    S3_interaction = list(
      scenario_id = "S3_interaction",
      scenario_label = "Strong correlation with ST interaction",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.8,
      rho_time = 0.8,
      gamma_interaction = 1.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.20, 0.20),
      create_hole_process2 = FALSE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 25,
      hole_t_max = 35,
      use_event_1 = TRUE,
      use_event_2 = TRUE,
      observation_pattern = "random"
    ),
    
    S4_intermittent = list(
      scenario_id = "S4_intermittent",
      scenario_label = "Y2 intermittently observed with final nowcasting gap",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.8,
      rho_time = 0.2,
      gamma_interaction = 1.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.10, 0.10),
      create_hole_process2 = FALSE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 25,
      hole_t_max = 35,
      use_event_1 = TRUE,
      use_event_2 = TRUE,
      observation_pattern = "intermittent",
      y2_frequency = 10,
      nowcasting_gap_start = 42,
      nowcasting_gap_end = 50
    ),
    
    S5_spatial_hole = list(
      scenario_id = "S5_spatial_hole",
      scenario_label = "Spatial hole in Y2 over all times",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.8,
      rho_time = 0.2,
      gamma_interaction = 1.0,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.20, 0.20),
      create_hole_process2 = TRUE,
      hole_x_min = 0.50,
      hole_x_max = 1.00,
      hole_y_min = 0.50,
      hole_y_max = 1.00,
      hole_t_min = 1,
      hole_t_max = 50,
      use_event_1 = FALSE,
      use_event_2 = TRUE,
      observation_pattern = "random"
    ),
    
    S6_ST_hole = list(
      scenario_id = "S6_ST_hole",
      scenario_label = "Spatio-temporal hole in Y2",
      seed_global = 123,
      seed_data = 456,
      alpha_shared = 0.8,
      rho_time = 0.2,
      gamma_interaction = 0.1,
      sigma_noise = c(0.10, 0.10),
      obs_prob_process = c(0.10, 0.10),
      create_hole_process2 = TRUE,
      hole_x_min = 0.00,
      hole_x_max = 1.00,
      hole_y_min = 0.00,
      hole_y_max = 1.00,
      hole_t_min = 20,
      hole_t_max = 35,
      use_event_1 = TRUE,
      use_event_2 = TRUE,
      observation_pattern = "random"
    )
  )
}

validate_scenario_config <- function(config) {
  required_names <- c(
    "scenario_id",
    "scenario_label",
    "seed_global",
    "seed_data",
    "alpha_shared",
    "rho_time",
    "gamma_interaction",
    "sigma_noise",
    "obs_prob_process",
    "create_hole_process2",
    "hole_x_min",
    "hole_x_max",
    "hole_y_min",
    "hole_y_max",
    "hole_t_min",
    "hole_t_max",
    "use_event_1",
    "use_event_2",
    "observation_pattern"
  )
  
  missing_names <- setdiff(required_names, names(config))
  
  if (length(missing_names) > 0L) {
    stop(
      "Configurazione scenario incompleta. Mancano: ",
      paste(missing_names, collapse = ", ")
    )
  }
  
  if (length(config$sigma_noise) != 2L) {
    stop("sigma_noise deve contenere un valore per ciascun processo.")
  }
  
  if (length(config$obs_prob_process) != 2L) {
    stop("obs_prob_process deve contenere un valore per ciascun processo.")
  }
  
  invisible(TRUE)
}

build_simulation_grid <- function(n_space_side, n_time) {
  x_values <- seq(0, 1, length.out = n_space_side)
  y_values <- seq(0, 1, length.out = n_space_side)
  
  space_grid <- expand.grid(
    x = x_values,
    y = y_values
  )
  space_grid$space_id <- seq_len(nrow(space_grid))
  
  time_grid <- data.frame(
    time_id = seq_len(n_time),
    t_scaled = seq(0, 1, length.out = n_time)
  )
  
  st_grid <- merge(
    space_grid,
    time_grid,
    by = NULL,
    sort = FALSE
  )
  st_grid$st_bau_id <- seq_len(nrow(st_grid))
  
  list(
    x_values = x_values,
    y_values = y_values,
    space_grid = space_grid,
    time_grid = time_grid,
    st_grid = st_grid
  )
}

simulate_latent_fields <- function(st_grid, config, beta) {
  set.seed(config$seed_data)
  
  shared_background <- with(
    st_grid,
    0.55 * sin(pi * x) +
      0.45 * cos(pi * y) +
      0.35 * cos(2 * pi * t_scaled)
  )
  
  shared_traffic <- with(
    st_grid,
    exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.10)
  )
  
  event_1 <- with(
    st_grid,
    0.9 * as.numeric(time_id >= 18 & time_id <= 25) *
      exp(-((x - 0.80)^2 + (y - 0.50)^2) / 0.030)
  )
  
  if (!config$use_event_1) {
    event_1[] <- 0
  }
  
  event_2_start <- 10
  event_2_end <- 30
  
  event_2_time <- pmin(
    pmax(
      (st_grid$time_id - event_2_start) /
        (event_2_end - event_2_start),
      0
    ),
    1
  )
  
  event_2_active <- as.numeric(
    st_grid$time_id >= event_2_start &
      st_grid$time_id <= event_2_end
  )
  
  front_center_x <- 0.18 + 0.62 * event_2_time
  front_width_x <- 0.055 + 0.035 * event_2_time
  front_shape_y <- 1 + 0.25 * sin(
    2 * pi * st_grid$y + 0.8 * pi * event_2_time
  )
  front_intensity <- sin(pi * event_2_time)
  
  event_2 <- event_2_active *
    1.35 *
    front_intensity *
    exp(-((st_grid$x - front_center_x)^2) / front_width_x) *
    front_shape_y
  
  if (!config$use_event_2) {
    event_2[] <- 0
  }
  
  shared_signal <- shared_background +
    0.35 * shared_traffic +
    0.55 * event_1 +
    1.10 * event_2
  
  shared_signal <- as.numeric(scale(shared_signal))
  
  specific_y1 <- with(
    st_grid,
    0.35 * sin(2 * pi * x) * cos(2 * pi * t_scaled) +
      0.20 * cos(2 * pi * y)
  )
  
  specific_y2 <- with(
    st_grid,
    0.35 * cos(2 * pi * y) * sin(2 * pi * t_scaled) +
      0.20 * sin(2 * pi * x)
  )
  
  specific_y1 <- as.numeric(scale(specific_y1))
  specific_y2 <- as.numeric(scale(specific_y2))
  
  alpha_shared <- config$alpha_shared
  specific_weight <- sqrt(1 - alpha_shared^2)
  
  y1_latent <- alpha_shared * shared_signal +
    specific_weight * specific_y1
  
  y2_latent <- alpha_shared * shared_signal +
    specific_weight * specific_y2
  
  y1_latent <- y1_latent + rnorm(length(y1_latent), sd = 0.05)
  y2_latent <- y2_latent + rnorm(length(y2_latent), sd = 0.05)
  
  y1_latent <- as.numeric(scale(y1_latent))
  y2_latent <- as.numeric(scale(y2_latent))
  
  trend_y1 <- with(
    st_grid,
    beta$intercept[1] +
      beta$x[1] * x +
      beta$y[1] * y +
      beta$time[1] * t_scaled
  )
  
  trend_y2 <- with(
    st_grid,
    beta$intercept[2] +
      beta$x[2] * x +
      beta$y[2] * y +
      beta$time[2] * t_scaled
  )
  
  y1_true <- y1_latent + trend_y1
  y2_true <- y2_latent + trend_y2
  
  list(
    y1_true = y1_true,
    y2_true = y2_true,
    trend_y1 = trend_y1,
    trend_y2 = trend_y2,
    shared_signal = shared_signal,
    event_1 = event_1,
    event_2 = event_2
  )
}

build_full_data <- function(st_grid, latent_fields, sigma_noise) {
  n_st <- nrow(st_grid)
  
  y1_data <- data.frame(
    process = "Y1",
    process_id = 1L,
    x = st_grid$x,
    y = st_grid$y,
    space_id = st_grid$space_id,
    time_id = st_grid$time_id,
    t_scaled = st_grid$t_scaled,
    st_bau_id = st_grid$st_bau_id,
    mean_true = latent_fields$trend_y1,
    w_true = latent_fields$y1_true - latent_fields$trend_y1,
    Y_true = latent_fields$y1_true,
    Z_obs = latent_fields$y1_true + rnorm(n_st, sd = sigma_noise[1])
  )
  
  y2_data <- data.frame(
    process = "Y2",
    process_id = 2L,
    x = st_grid$x,
    y = st_grid$y,
    space_id = st_grid$space_id,
    time_id = st_grid$time_id,
    t_scaled = st_grid$t_scaled,
    st_bau_id = st_grid$st_bau_id,
    mean_true = latent_fields$trend_y2,
    w_true = latent_fields$y2_true - latent_fields$trend_y2,
    Y_true = latent_fields$y2_true,
    Z_obs = latent_fields$y2_true + rnorm(n_st, sd = sigma_noise[2])
  )
  
  data_full <- rbind(y1_data, y2_data)
  
  # Retain the legacy observation alias for downstream compatibility.
  data_full$Y_obs <- data_full$Z_obs
  
  data_full
}

apply_observation_pattern <- function(
    data_full,
    config,
    n_time,
    observation_seed = NULL
) {
  # Sample monitoring locations for each process; retain them across time.
  data_full$observed <- FALSE
  
  if (!is.null(observation_seed)) {
    set.seed(observation_seed)
  }
  
  selected_space_ids <- vector("list", 2L)
  
  for (process_id in 1:2) {
    available_space_ids <- sort(unique(
      data_full$space_id[
        data_full$process_id == process_id
      ]
    ))
    
    n_observed_space <- as.integer(round(
      config$obs_prob_process[process_id] *
        length(available_space_ids)
    ))
    
    n_observed_space <- max(
      0L,
      min(length(available_space_ids), n_observed_space)
    )
    
    if (n_observed_space == 0L) {
      selected_space_ids[[process_id]] <- integer(0)
    } else {
      selected_space_ids[[process_id]] <- sample(
        available_space_ids,
        size = n_observed_space,
        replace = FALSE
      )
    }
  }
  
  # Intermittent sampling restricts Y2 to selected observation times.
  if (identical(config$observation_pattern, "intermittent")) {
    
    idx_y1 <- which(
      data_full$process_id == 1L &
        data_full$space_id %in% selected_space_ids[[1L]]
    )
    data_full$observed[idx_y1] <- TRUE
    
    y2_available_times <- seq(
      from = 1,
      to = n_time,
      by = config$y2_frequency
    )
    
    idx_y2_available <- which(
      data_full$process_id == 2L &
        data_full$space_id %in% selected_space_ids[[2L]] &
        data_full$time_id %in% y2_available_times
    )
    data_full$observed[idx_y2_available] <- TRUE
    
    idx_nowcasting_gap <- which(
      data_full$process_id == 2L &
        data_full$time_id >= config$nowcasting_gap_start &
        data_full$time_id <= config$nowcasting_gap_end
    )
    data_full$observed[idx_nowcasting_gap] <- FALSE
    
  } else if (identical(config$observation_pattern, "random")) {
    
    for (process_id in 1:2) {
      idx_process <- which(
        data_full$process_id == process_id &
          data_full$space_id %in% selected_space_ids[[process_id]]
      )
      
      data_full$observed[idx_process] <- TRUE
    }
    
  } else {
    stop(
      "Pattern osservativo non riconosciuto: ",
      config$observation_pattern
    )
  }
  
  # Remove Y2 observations inside the configured spatio-temporal gap.
  if (isTRUE(config$create_hole_process2)) {
    idx_hole <- which(
      data_full$process_id == 2L &
        data_full$x >= config$hole_x_min &
        data_full$x <= config$hole_x_max &
        data_full$y >= config$hole_y_min &
        data_full$y <= config$hole_y_max &
        data_full$time_id >= config$hole_t_min &
        data_full$time_id <= config$hole_t_max
    )
    
    data_full$observed[idx_hole] <- FALSE
  }
  
  attr(data_full, "selected_space_ids") <- selected_space_ids
  data_full
}

resample_observations <- function(
    simulation,
    observation_seed
) {
  
  simulation_new <- simulation
  
  data_full <- apply_observation_pattern(
    data_full = simulation$data_full,
    config = simulation$config,
    n_time = simulation$dimensions$n_time,
    observation_seed = observation_seed
  )
  
  simulation_new$data_full <- data_full
  
  simulation_new$data_obs <- data_full[
    data_full$observed,
    ,
    drop = FALSE
  ]
  
  simulation_new$observation_seed <- observation_seed
  
  simulation_new
}

print_simulation_summary <- function(simulation) {
  config <- simulation$config
  data_obs <- simulation$data_obs
  
  cat("\n============================================\n")
  cat("Scenario:", config$scenario_id, "\n")
  cat(config$scenario_label, "\n")
  cat("============================================\n")
  cat("BAU spaziali:", nrow(simulation$space_grid), "\n")
  cat("Tempi:", nrow(simulation$time_grid), "\n")
  cat("BAU spazio-temporali:", nrow(simulation$st_grid), "\n")
  cat("Correlazione empirica Y1-Y2:",
      round(simulation$empirical_process_correlation, 4), "\n")
  cat("Osservazioni per processo:\n")
  print(table(data_obs$process))
  
  invisible(simulation)
}

simulate_scenario <- function(
    scenario_name = "S6_ST_hole",
    n_space_side = 20,
    n_time = 50,
    n_basis_space_side = 5,
    n_basis_time = 10,
    kappa_space = 0.2,
    beta = list(
      intercept = c(0, 0),
      x = c(0, 0),
      y = c(0, 0),
      time = c(0, 0)
    ),
    observation_seed = NULL,
    verbose = TRUE
) {
  scenarios <- get_simulation_scenarios()
  config <- scenarios[[scenario_name]]
  
  if (is.null(config)) {
    stop(
      "Scenario non trovato: ", scenario_name,
      ". Scenari disponibili: ",
      paste(names(scenarios), collapse = ", ")
    )
  }
  
  validate_scenario_config(config)
  set.seed(config$seed_global)
  
  grids <- build_simulation_grid(
    n_space_side = n_space_side,
    n_time = n_time
  )
  
  latent_fields <- simulate_latent_fields(
    st_grid = grids$st_grid,
    config = config,
    beta = beta
  )
  
  data_full <- build_full_data(
    st_grid = grids$st_grid,
    latent_fields = latent_fields,
    sigma_noise = config$sigma_noise
  )
  
  data_full <- apply_observation_pattern(
    data_full = data_full,
    config = config,
    n_time = n_time,
    observation_seed = observation_seed
  )
  
  selected_space_ids <- attr(
    data_full,
    "selected_space_ids"
  )
  
  data_obs <- data_full[data_full$observed, , drop = FALSE]
  
  simulation <- list(
    scenario_name = scenario_name,
    config = config,
    dimensions = list(
      n_process = 2L,
      n_space_side = n_space_side,
      n_space = nrow(grids$space_grid),
      n_time = n_time
    ),
    basis_settings = list(
      n_basis_space_side = n_basis_space_side,
      n_basis_time = n_basis_time
    ),
    precision_parameters = list(
      kappa_space = kappa_space,
      rho_time = config$rho_time,
      gamma_interaction = config$gamma_interaction
    ),
    sigma_noise = config$sigma_noise,
    rho_process_true_approx = config$alpha_shared^2,
    empirical_process_correlation = cor(
      latent_fields$y1_true,
      latent_fields$y2_true
    ),
    observation_seed = observation_seed,
    selected_space_ids = selected_space_ids,
    x_values = grids$x_values,
    y_values = grids$y_values,
    space_grid = grids$space_grid,
    time_grid = grids$time_grid,
    st_grid = grids$st_grid,
    latent_fields = latent_fields,
    data_full = data_full,
    data_obs = data_obs
  )
  
  if (verbose) {
    print_simulation_summary(simulation)
  }
  
  simulation
}
