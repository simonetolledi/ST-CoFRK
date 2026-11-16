############################################################
# 06_update_precision_parameters.R
# Update ECM dei parametri di Q_ST con scala identificata in Sigma_process.
############################################################

safe_sparse_logdet <- function(Q) {
  Q <- methods::as((Q + Matrix::t(Q)) / 2, "dgCMatrix")
  out <- Matrix::determinant(Q, logarithm = TRUE)
  value <- as.numeric(out$modulus)
  if (out$sign <= 0 || !is.finite(value)) return(-Inf)
  value
}

# Rimuove la scala globale di Q_ST: la scala resta in Sigma_process.
normalize_precision_scale <- function(Q) {
  Q <- methods::as((Q + Matrix::t(Q)) / 2, "dgCMatrix")
  logdet <- safe_sparse_logdet(Q)
  if (!is.finite(logdet)) stop("Q_ST non SPD durante la normalizzazione.")
  scale_factor <- exp(logdet / nrow(Q))
  list(Q = methods::as(Q / scale_factor, "dgCMatrix"),
       scale_factor = scale_factor)
}

generate_rademacher_probes <- function(n, n_probes = 64L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(sample(c(-1, 1), n * n_probes, replace = TRUE),
         nrow = n, ncol = n_probes)
}

apply_separable_precision <- function(
    x,
    Sigma_inverse,
    Q_ST,
    n_process,
    n_basis_st
) {
  
  input_is_vector <- is.null(dim(x))
  
  if (input_is_vector) {
    x <- matrix(x, ncol = 1L)
  }
  
  if (nrow(x) != n_process * n_basis_st) {
    stop("Dimensioni incompatibili.")
  }
  
  result <- matrix(
    0,
    nrow(x),
    ncol(x)
  )
  
  Q_ST <- methods::as(
    Q_ST,
    "dgCMatrix"
  )
  
  Sigma_inverse <- as.matrix(
    Sigma_inverse
  )
  
  for (j in seq_len(ncol(x))) {
    
    X <- matrix(
      x[, j],
      nrow = n_process,
      ncol = n_basis_st,
      byrow = TRUE
    )
    
    tmp <- as.matrix(
      Sigma_inverse %*%
        X %*%
        Q_ST
    )
    
    result[, j] <- as.vector(
      t(tmp)
    )
  }
  
  if (input_is_vector) {
    as.numeric(result[, 1L])
  } else {
    result
  }
}

solve_posterior_probes <- function(posterior_chol, trace_probes) {
  as.matrix(Matrix::solve(posterior_chol, trace_probes, system = "A"))
}

estimate_hutchinson_trace <- function(trace_probes, posterior_probe_solutions,
                                      Sigma_inverse, Q_ST,
                                      n_process, n_basis_st) {
  QVZ <- apply_separable_precision(
    posterior_probe_solutions, Sigma_inverse, Q_ST,
    n_process, n_basis_st
  )
  mean(colSums(trace_probes * QVZ))
}

expected_precision_objective <- function(eta_hat, Sigma_process, Q_ST,
                                         trace_probes,
                                         posterior_probe_solutions,
                                         return_details = FALSE) {
  n_process <- nrow(Sigma_process)
  if (length(eta_hat) %% n_process != 0L) stop("eta_hat incompatibile.")
  n_basis_st <- length(eta_hat) / n_process
  Sigma_inverse <- solve((Sigma_process + t(Sigma_process)) / 2)
  
  logdet_term <- 0.5 * n_process * safe_sparse_logdet(Q_ST)
  if (!is.finite(logdet_term)) {
    details <- list(value = -Inf, logdet_term = -Inf,
                    mean_term = Inf, trace_term = Inf)
    if (return_details) return(details)
    return(-Inf)
  }
  
  Qmean <- apply_separable_precision(
    eta_hat, Sigma_inverse, Q_ST, n_process, n_basis_st
  )
  mean_term <- as.numeric(crossprod(eta_hat, Qmean))
  trace_term <- estimate_hutchinson_trace(
    trace_probes, posterior_probe_solutions, Sigma_inverse,
    Q_ST, n_process, n_basis_st
  )
  value <- logdet_term - 0.5 * (mean_term + trace_term)
  details <- list(
    value = if (is.finite(value)) value else -Inf,
    logdet_term = logdet_term,
    mean_term = mean_term,
    trace_term = trace_term
  )
  if (return_details) return(details)
  details$value
}

precision_parameter_penalty <- function(parameters, penalty) {
  if (is.null(penalty)) return(0)
  kappa <- parameters$kappa_space
  rho <- parameters$rho_time
  gamma <- parameters$gamma_interaction
  
  0.5 * penalty$lambda_kappa *
    (log(kappa) - log(penalty$kappa_reference))^2 +
    0.5 * penalty$lambda_rho *
    (log(rho) - log(penalty$rho_reference))^2 +
    0.5 * penalty$lambda_gamma *
    (log(gamma) - log(penalty$gamma_reference))^2
}

evaluate_precision_parameters <- function(candidate_parameters, eta_hat,
                                          Sigma_process, basis,
                                          trace_probes,
                                          posterior_probe_solutions,
                                          ridge,
                                          penalty = NULL,
                                          return_details = FALSE) {
  candidate <- build_precision_matrices(
    basis = basis,
    kappa_space = candidate_parameters$kappa_space,
    rho_time = candidate_parameters$rho_time,
    gamma_interaction = candidate_parameters$gamma_interaction,
    ridge = ridge,
    verbose = FALSE
  )
  
  objective_details <- expected_precision_objective(
    eta_hat, Sigma_process, candidate$Q_ST,
    trace_probes, posterior_probe_solutions,
    return_details = TRUE
  )
  raw <- objective_details$value
  pen <- precision_parameter_penalty(candidate_parameters, penalty)
  value <- raw - pen
  
  if (return_details) {
    return(list(
      value = value,
      raw = raw,
      penalty = pen,
      Q_ST = candidate$Q_ST,
      scale_factor = 1,
      logdet_term = objective_details$logdet_term,
      mean_term = objective_details$mean_term,
      trace_term = objective_details$trace_term
    ))
  }
  value
}

.damp_parameters <- function(old, proposed, damping) {
  out <- old
  out$kappa_space <- exp((1 - damping) * log(old$kappa_space) +
                           damping * log(proposed$kappa_space))
  out$rho_time <- exp((1 - damping) * log(old$rho_time) +
                           damping * log(proposed$rho_time))
  out$gamma_interaction <- exp((1 - damping) * log(old$gamma_interaction) +
                                 damping * log(proposed$gamma_interaction))
  out
}

update_precision_jointly <- function(
    current_parameters,
    eta_hat,
    Sigma_process,
    basis,
    trace_probes,
    posterior_probe_solutions,
    update_kappa = TRUE,
    update_rho = TRUE,
    update_gamma = TRUE,
    kappa_bounds,
    rho_bounds,
    gamma_bounds,
    ridge,
    penalty,
    optim_control
) {
  active_parameters <- c(
    if (update_kappa) "kappa_space",
    if (update_rho) "rho_time",
    if (update_gamma) "gamma_interaction"
  )
  
  if (length(active_parameters) == 0L) {
    return(current_parameters)
  }
  
  if (update_gamma && gamma_bounds[1L] <= 0) {
    stop("Con la trasformazione log(gamma), gamma_bounds[1] deve essere > 0.")
  }
  
  transform_parameter <- function(name, value) {
    switch(
      name,
      kappa_space = log(value),
      rho_time = log(value),
      gamma_interaction = log(value),
      stop("Parametro non riconosciuto.")
    )
  }
  
  inverse_transform_parameter <- function(name, value) {
    switch(
      name,
      kappa_space = exp(value),
      rho_time = exp(value),
      gamma_interaction = exp(value),
      stop("Parametro non riconosciuto.")
    )
  }
  
  transformed_bounds <- list(
    kappa_space = log(kappa_bounds),
    rho_time = log(rho_bounds),
    gamma_interaction = log(gamma_bounds)
  )
  
  initial_values <- vapply(
    active_parameters,
    function(name) {
      transform_parameter(name, current_parameters[[name]])
    },
    numeric(1L)
  )
  
  lower_bounds <- vapply(
    active_parameters,
    function(name) transformed_bounds[[name]][1L],
    numeric(1L)
  )
  
  upper_bounds <- vapply(
    active_parameters,
    function(name) transformed_bounds[[name]][2L],
    numeric(1L)
  )
  
  objective_function <- function(transformed_values) {
    candidate_parameters <- current_parameters
    
    for (j in seq_along(active_parameters)) {
      parameter_name <- active_parameters[j]
      candidate_parameters[[parameter_name]] <-
        inverse_transform_parameter(
          parameter_name,
          transformed_values[j]
        )
    }
    
    objective_value <- evaluate_precision_parameters(
      candidate_parameters = candidate_parameters,
      eta_hat = eta_hat,
      Sigma_process = Sigma_process,
      basis = basis,
      trace_probes = trace_probes,
      posterior_probe_solutions = posterior_probe_solutions,
      ridge = ridge,
      penalty = penalty
    )
    
    if (!is.finite(objective_value)) {
      return(.Machine$double.xmax / 100)
    }
    
    -objective_value
  }
  
  optimization_result <- stats::optim(
    par = initial_values,
    fn = objective_function,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = optim_control
  )
  
  proposed_parameters <- current_parameters
  
  for (j in seq_along(active_parameters)) {
    parameter_name <- active_parameters[j]
    proposed_parameters[[parameter_name]] <-
      inverse_transform_parameter(
        parameter_name,
        optimization_result$par[j]
      )
  }
  
  attr(proposed_parameters, "joint_precision_optim") <-
    optimization_result
  
  proposed_parameters
}

update_precision_parameters <- function(
    current_parameters, eta_hat, Sigma_process, basis,
    trace_probes, posterior_probe_solutions,
    update_kappa = TRUE, update_rho = TRUE, update_gamma = TRUE,
    kappa_bounds, rho_bounds,
    gamma_bounds, ridge = 1e-8,
    penalty = NULL, damping = 0.5,
    optim_control = list(maxit = 100L, factr = 1e7, pgtol = 1e-6),
    verbose = TRUE) {
  
  n_eta <- length(eta_hat)
  if (n_eta %% nrow(Sigma_process) != 0L) stop("eta_hat incompatibile.")
  trace_probes <- as.matrix(trace_probes)
  posterior_probe_solutions <- as.matrix(posterior_probe_solutions)
  if (nrow(trace_probes) != n_eta ||
      nrow(posterior_probe_solutions) != n_eta ||
      ncol(trace_probes) != ncol(posterior_probe_solutions)) {
    stop("Dimensioni Hutchinson incompatibili.")
  }
  
  if (update_kappa && kappa_bounds[1L] <= 0) {
    stop("kappa_bounds[1] deve essere > 0.")
  }
  if (update_rho && rho_bounds[1L] <= 0) {
    stop("rho_bounds[1] deve essere > 0.")
  }
  if (update_gamma && gamma_bounds[1L] <= 0) {
    stop("gamma_bounds[1] deve essere > 0 con la trasformazione log(gamma).")
  }
  if (update_gamma && current_parameters$gamma_interaction <= 0) {
    stop("Il valore iniziale di gamma_interaction deve essere > 0.")
  }
  
  # penalty = NULL significa realmente nessuna penalizzazione.
  
  before <- evaluate_precision_parameters(
    current_parameters, eta_hat, Sigma_process, basis,
    trace_probes, posterior_probe_solutions, ridge, penalty,
    return_details = TRUE
  )
  proposed <- update_precision_jointly(
    current_parameters = current_parameters,
    eta_hat = eta_hat,
    Sigma_process = Sigma_process,
    basis = basis,
    trace_probes = trace_probes,
    posterior_probe_solutions = posterior_probe_solutions,
    update_kappa = update_kappa,
    update_rho = update_rho,
    update_gamma = update_gamma,
    kappa_bounds = kappa_bounds,
    rho_bounds = rho_bounds,
    gamma_bounds = gamma_bounds,
    ridge = ridge,
    penalty = penalty,
    optim_control = optim_control
  )
  
  proposed <- .damp_parameters(current_parameters, proposed, damping)
  after <- evaluate_precision_parameters(
    proposed, eta_hat, Sigma_process, basis,
    trace_probes, posterior_probe_solutions, ridge, penalty,
    return_details = TRUE
  )
  
  min_gain <- 1e-8 * max(1, abs(before$value))
  accepted <- is.finite(after$value) && after$value >= before$value + min_gain
  updated <- if (accepted) proposed else current_parameters
  final <- if (accepted) after else before
  
  attr(updated, "expected_precision_objective") <- final$value
  attr(updated, "raw_precision_objective") <- final$raw
  attr(updated, "precision_penalty") <- final$penalty
  attr(updated, "precision_update_accepted") <- accepted
  attr(updated, "Q_ST_evaluated") <- final$Q_ST
  attr(updated, "Q_ST_normalized") <- final$Q_ST
  attr(updated, "Q_ST_scale_factor") <- final$scale_factor
  attr(updated, "precision_logdet_term") <- final$logdet_term
  attr(updated, "precision_mean_term") <- final$mean_term
  attr(updated, "precision_trace_term") <- final$trace_term
  boundary_tolerance <- 1e-6
  attr(updated, "precision_at_boundary") <- c(
    kappa = min(abs(updated$kappa_space - kappa_bounds)) < boundary_tolerance,
    rho = min(abs(updated$rho_time - rho_bounds)) < boundary_tolerance,
    gamma = min(abs(updated$gamma_interaction - gamma_bounds)) < boundary_tolerance
  )
  
  if (verbose) {
    cat("\n---------------------------------------------\n")
    cat(" Precision parameter update\n")
    cat("---------------------------------------------\n")
    
    cat(
      sprintf(
        "%-12s %-14s\n",
        "Parameter",
        "Value"
      )
    )
    
    cat(
      sprintf(
        "%-12s %-14.6f\n",
        "kappa",
        updated$kappa_space
      )
    )
    
    cat(
      sprintf(
        "%-12s %-14.6f\n",
        "rho",
        updated$rho_time
      )
    )
    
    cat(
      sprintf(
        "%-12s %-14.6f\n",
        "gamma",
        updated$gamma_interaction
      )
    )
    
    cat(
      "objective:", round(final$value, 4),
      "| logdet:", round(final$logdet_term, 4),
      "| mean:", round(final$mean_term, 4),
      "| trace:", round(final$trace_term, 4), "\n"
    )
    boundary <- attr(updated, "precision_at_boundary")
    if (any(boundary)) {
      cat("boundary warning:",
          paste(names(boundary)[boundary], collapse = ", "), "\n")
    }
    cat("---------------------------------------------\n\n")
  }
  updated
}

