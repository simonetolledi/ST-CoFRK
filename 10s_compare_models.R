############################################################
# 10_compare_models.R
#
# Costruzione e salvataggio del riepilogo dello scenario.
############################################################

build_scenario_summary <- function(simulation, fit, metrics) {
  Sigma <- as.matrix(fit$Sigma_process)
  n_process <- nrow(Sigma)
  if (!n_process %in% c(1L, 2L)) {
    stop("build_scenario_summary supporta uno o due processi.")
  }
  
  # Nel caso univariato, il processo interno 1 corrisponde a Y2.
  get_sigma <- function(i, j) {
    if (n_process == 1L) {
      if (i == 2L && j == 2L) return(as.numeric(Sigma[1L, 1L]))
      return(NA_real_)
    }
    as.numeric(Sigma[i, j])
  }
  get_fine <- function(i) {
    if (n_process == 1L) {
      if (i == 2L) return(as.numeric(fit$sigma_fine[1L])^2)
      return(NA_real_)
    }
    as.numeric(fit$sigma_fine[i])^2
  }
  get_config <- function(name, index = NULL) {
    x <- simulation$config[[name]]
    if (!is.null(index) && !is.null(x)) x <- x[index]
    if (length(x) != 1L || is.na(x)) return(NA)
    x
  }
  get_metric <- function(tab, name) {
    row <- tab[as.character(tab$process) == "Y2", , drop = FALSE]
    if (nrow(row) != 1L) {
      stop("Attesa esattamente una riga Y2 in metrics$global/missing.")
    }
    as.numeric(row[[name]][1L])
  }
  
  mse_frk_global <- get_metric(metrics$global, "MSE_FRK_univariato")
  mse_cofrk_global <- get_metric(metrics$global, "MSE_COFRK_multivariato")
  mse_frk_missing <- get_metric(metrics$missing, "MSE_FRK_univariato")
  mse_cofrk_missing <- get_metric(metrics$missing, "MSE_COFRK_multivariato")
  improvement_missing <- if (is.finite(mse_frk_missing) && mse_frk_missing > 0) {
    100 * (mse_frk_missing - mse_cofrk_missing) / mse_frk_missing
  } else {
    NA_real_
  }
  
  config <- simulation$config
  rho_hat <- if (n_process == 2L && length(fit$rho_process_hat) == 1L) {
    as.numeric(fit$rho_process_hat)
  } else {
    NA_real_
  }
  rho_true <- simulation$rho_process_true_approx
  if (length(rho_true) != 1L) rho_true <- NA_real_
  seed <- simulation$observation_seed
  if (length(seed) != 1L) seed <- NA_integer_
  
  data.frame(
    scenario_id = get_config("scenario_id"),
    scenario_label = get_config("scenario_label"),
    observation_seed = seed,
    alpha_shared = get_config("alpha_shared"),
    rho_process_true_approx = rho_true,
    rho_process_hat = rho_hat,
    sigma_process_11_hat = get_sigma(1L, 1L),
    sigma_process_22_hat = get_sigma(2L, 2L),
    sigma_process_12_hat = get_sigma(1L, 2L),
    sigma_fine_Y1_variance_hat = get_fine(1L),
    sigma_fine_Y2_variance_hat = get_fine(2L),
    kappa_space_hat = fit$precision_parameters$kappa_space,
    rho_time_hat = fit$precision_parameters$rho_time,
    gamma_interaction_hat = fit$precision_parameters$gamma_interaction,
    obs_prob_Y1 = get_config("obs_prob_process", 1L),
    obs_prob_Y2 = get_config("obs_prob_process", 2L),
    use_event_1 = get_config("use_event_1"),
    use_event_2 = get_config("use_event_2"),
    create_hole_process2 = get_config("create_hole_process2"),
    hole_x_min = get_config("hole_x_min"),
    hole_x_max = get_config("hole_x_max"),
    hole_y_min = get_config("hole_y_min"),
    hole_y_max = get_config("hole_y_max"),
    hole_t_min = get_config("hole_t_min"),
    hole_t_max = get_config("hole_t_max"),
    mse_global_FRK_Y2 = mse_frk_global,
    mse_global_COFRK_Y2 = mse_cofrk_global,
    mse_missing_FRK_Y2 = mse_frk_missing,
    mse_missing_COFRK_Y2 = mse_cofrk_missing,
    improvement_missing_Y2 = improvement_missing
  )
}

append_scenario_summary <- function(
    summary_row,
    summary_path
) {
  dir.create(
    dirname(summary_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  utils::write.table(
    summary_row,
    file = summary_path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(summary_path),
    append = file.exists(summary_path),
    quote = TRUE
  )
  
  invisible(summary_path)
}


summarize_monte_carlo <- function(summary_runs) {
  numeric_columns <- c(
    "rho_process_hat",
    "sigma_process_11_hat",
    "sigma_process_22_hat",
    "sigma_process_12_hat",
    "sigma_fine_Y1_variance_hat",
    "sigma_fine_Y2_variance_hat",
    "kappa_space_hat",
    "rho_time_hat",
    "gamma_interaction_hat",
    "mse_global_FRK_Y2",
    "mse_global_COFRK_Y2",
    "mse_missing_FRK_Y2",
    "mse_missing_COFRK_Y2",
    "improvement_missing_Y2"
  )
  
  available_columns <- intersect(
    numeric_columns,
    names(summary_runs)
  )
  
  if (length(available_columns) == 0L) {
    stop(
      "Nessuna colonna numerica disponibile per ",
      "il riepilogo Monte Carlo."
    )
  }
  
  out <- do.call(
    rbind,
    lapply(
      available_columns,
      function(variable) {
        values <- summary_runs[[variable]]
        valid_values <- values[is.finite(values)]
        
        if (length(valid_values) == 0L) {
          return(data.frame(
            parameter = variable,
            n = 0L,
            mean = NA_real_,
            variance = NA_real_,
            sd = NA_real_,
            min = NA_real_,
            max = NA_real_
          ))
        }
        
        data.frame(
          parameter = variable,
          n = length(valid_values),
          mean = mean(valid_values),
          variance = if (length(valid_values) > 1L) {
            stats::var(valid_values)
          } else {
            0
          },
          sd = if (length(valid_values) > 1L) {
            stats::sd(valid_values)
          } else {
            0
          },
          min = min(valid_values),
          max = max(valid_values)
        )
      }
    )
  )
  
  rownames(out) <- NULL
  out
}