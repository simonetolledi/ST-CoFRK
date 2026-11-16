
############################################################
# 11s_run_example.R
#
# Monte Carlo simulation pipeline for ST-CoFRK.
# The latent fields, measurement noise, and scenario
# parameters remain fixed across runs. Only the
# observation pattern changes.
############################################################

source("00_load_all.R")

# Select one of the available simulation scenarios.
scenario_name <- "S0_independent"

#"S0_independent"
#"S1_weak_correlation"
#"S2_strong_correlation"
#"S3_interaction"
#"S4_intermittent"
#"S5_spatial_hole"
#"S6_ST_hole"

n_runs <- 20L

observation_seeds <- 1000L + seq_len(n_runs)

results_root <- file.path(
  getwd(),
  "results_scenarios"
)


# Build the latent fields, basis functions, and initial
# precision matrices once, before resampling observations.

simulation_base <- simulate_scenario(
  scenario_name = scenario_name,
  observation_seed = observation_seeds[1L],
  verbose = FALSE
)

basis <- generate_multires_st_basis(
  simulation = simulation_base,
  verbose = FALSE
)

precision <- build_precision_matrices(
  basis = basis,
  kappa_space =
    simulation_base$precision_parameters$kappa_space,
  rho_time =
    simulation_base$precision_parameters$rho_time,
  gamma_interaction =
    simulation_base$precision_parameters$gamma_interaction,
  verbose = FALSE
)


run_single_simulation <- function(
    simulation_base,
    basis,
    precision,
    observation_seed,
    verbose = FALSE
) {
  
  simulation <- resample_observations(
    simulation = simulation_base,
    observation_seed = observation_seed
  )
  
  simulation_y2 <- simulation
  
  # Optional univariate configuration: retain only Y2
  # and reindex it as process 1.
  #-----------------------------------------------
  #simulation_y2$data_full <- simulation$data_full[
  #  simulation$data_full$process_id == 2L, ,
  #  drop = FALSE
  #]
  
  #simulation_y2$data_obs <- simulation$data_obs[
  #  simulation$data_obs$process_id == 2L, ,
  #  drop = FALSE
  #]
  
  #simulation_y2$data_full$process_id <- 1L
  #simulation_y2$data_obs$process_id <- 1L
  
  #simulation_y2$dimensions$n_process <- 1L
  #simulation_y2$sigma_noise <- simulation$sigma_noise[2L]
  #--------------------------------------------------------------
  
  model <- build_model_object(
    simulation = simulation_y2,
    basis = basis,
    mean_formula = ~ 1,
    verbose = verbose
  )
  
  fit <- em_fit(
    model = model,
    precision = precision,
    basis = basis,
    n_em = 20,
    tolerance = 1e-4,
    update_beta = TRUE,
    update_precision = TRUE,
    update_precision_every = 1L,
    verbose = verbose
  )
  
  predictions <- predict_cofrk(
    model = model,
    fit = fit,
    basis = basis,
    verbose = verbose
  )
  
  benchmark <- fit_FRK_benchmark(
    simulation = simulation_y2,
    basis = basis,
    process_ids = 1L,
    verbose = verbose
  )
  
  comparison_data <- predictions$data_full
  
  comparison_data$Y_hat_FRK <-
    benchmark$data_full$Y_hat_FRK
  
  metrics <- compute_model_metrics(
    data_full = comparison_data
  )
  
  summary_row <- build_scenario_summary(
    simulation = simulation_y2,
    fit = fit,
    metrics = metrics
  )
  
  list(
    observation_seed = observation_seed,
    simulation = simulation_y2,
    fit = fit,
    metrics = metrics,
    summary = summary_row
  )
}



# Run the Monte Carlo experiment.

run_results <- vector(
  "list",
  length = n_runs
)

for (run_id in seq_len(n_runs)) {
  
  cat(
    "\n============================================\n",
    "RUN ",
    run_id,
    "/",
    n_runs,
    " | observation_seed = ",
    observation_seeds[run_id],
    "\n============================================\n",
    sep = ""
  )
  
  run_results[[run_id]] <- run_single_simulation(
    simulation_base = simulation_base,
    basis = basis,
    precision = precision,
    observation_seed = observation_seeds[run_id],
    verbose = TRUE
  )
}


# Collect run-level estimates and aggregate Monte Carlo results.

summary_runs <- do.call(
  rbind,
  lapply(
    seq_along(run_results),
    function(run_id) {
      out <- run_results[[run_id]]$summary
      out$run_id <- run_id
      out
    }
  )
)

summary_runs <- summary_runs[
  ,
  c(
    "run_id",
    setdiff(names(summary_runs), "run_id")
  ),
  drop = FALSE
]


monte_carlo_summary <- summarize_monte_carlo(
  summary_runs = summary_runs
)


aggregated <- aggregate_monte_carlo_results(
  run_results = run_results
)


scenario_label <-
  run_results[[1L]]$simulation$config$scenario_label

plots <- make_monte_carlo_plots(
  aggregated = aggregated,
  scenario_label = scenario_label
)


scenario_directory <- file.path(
  results_root,
  scenario_name,
  paste0("monte_carlo_", n_runs, "_runs")
)

save_monte_carlo_results(
  aggregated = aggregated,
  plots = plots,
  summary_runs = summary_runs,
  monte_carlo_summary = monte_carlo_summary,
  output_directory = scenario_directory
)


cat("\n============================================\n")
cat("RIEPILOGO MONTE CARLO\n")
cat("============================================\n")

print(monte_carlo_summary)

cat("\nPrime righe del riepilogo per run:\n")
print(utils::head(summary_runs))

cat("\nRisultati salvati in:\n")
cat(scenario_directory, "\n")


plots[["comparison_map"]]
plots[["prediction_variance_map"]]
plots[["mse_time"]]
plots[["mse_variance"]]
plots[["time_locations"]]


plot_obs <- plot_observed_points_times(
  simulation = run_results[[1]]$simulation,
  process_id = 1L,
  times = c(1L, 11L, 21L, 31L, 41L, 50L)
)

print(plot_obs)



