############################################################
# 12_parameter_recovery.R
# Parameter recovery under the correctly specified ST-CoFRK
# data-generating model. Run from the project directory.
# Requires the existing 00_load_all.R and modules 01--06.
# No FRK benchmark is fitted here.
############################################################

source("00_load_all.R")

# -------------------- USER SETTINGS -----------------------
n_runs <- 20L                 # PILOT: 1; after inspection set to 20
seed_start <- 20261110L
output_dir <- file.path(getwd(), "results_parameter_recovery")

n_space_side <- 20L
n_time <- 50L
n_basis_space_side <- 5L
n_basis_time <- 8L
obs_prob_process <- c(0.8, 0.8)  # spatial monitoring sites, observed over all times

truth <- list(
  kappa_space = 0.2,
  rho_time = 0.2,
  gamma_interaction = 0.5,
  Sigma_process = matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2),
  sigma_fine = c(0.05, 0.05),       # standard deviations, NOT variances
  sigma_noise = c(0.10, 0.10)       # known measurement-error SDs
)

n_em <- 25L
n_probes <- 32L
# -------------------- VALIDATION --------------------------
stopifnot(
  n_runs >= 1L,
  all(obs_prob_process > 0 & obs_prob_process <= 1),
  all(truth$sigma_fine >= 0), all(truth$sigma_noise > 0),
  all(eigen(truth$Sigma_process, symmetric = TRUE)$values > 0),
  abs(sum(diag(truth$Sigma_process)) - 2) < 1e-10
)

# Build the SAME grid, bases and precision as the fitting code.
# The old simulation is used ONLY as a grid/basis scaffold;
# its synthetic fields are discarded completely.
scaffold <- simulate_scenario(
  scenario_name = "S0_independent",
  n_space_side = n_space_side,
  n_time = n_time,
  n_basis_space_side = n_basis_space_side,
  n_basis_time = n_basis_time,
  observation_seed = seed_start,
  verbose = FALSE
)
scaffold$config$obs_prob_process <- obs_prob_process
scaffold$config$observation_pattern <- "random"
scaffold$config$create_hole_process2 <- FALSE
scaffold$sigma_noise <- truth$sigma_noise

basis <- generate_multires_st_basis(scaffold, verbose = FALSE)
precision_true <- build_precision_matrices(
  basis = basis,
  kappa_space = truth$kappa_space,
  rho_time = truth$rho_time,
  gamma_interaction = truth$gamma_interaction,
  verbose = FALSE
)

# Matrix convention: eta is stored process-by-process, with spatial
# basis varying fastest within each temporal basis. If U'U = Q_ST,
# solve(U, z) has covariance Q_ST^{-1} for z ~ N(0,I).
# Multiplying on the left by chol(Sigma_process)' introduces
# the process covariance: Cov(vec(t(E))) = Sigma_process %x% Q_ST^{-1}.
U <- chol(as.matrix(precision_true$Q_ST))
L_process <- t(chol(truth$Sigma_process))
K <- basis$dimensions$n_basis_st

# Precompute full-grid tensor-product design; process blocks are
# assembled by build_multivariate_design() (same routine as the fit).
full_template <- scaffold$data_full
A_full <- build_multivariate_design(
  data = full_template,
  Phi = basis$Phi,
  Psi = basis$Psi,
  n_process = 2L
)

simulate_from_model <- function(run_seed) {
  set.seed(run_seed)
  independent_coefficients <- matrix(rnorm(2L * K), nrow = 2L)
  # Rows = processes, columns = ST coefficients.
  spatial_temporal_coefficients <- t(backsolve(U, t(independent_coefficients)))
  eta_matrix <- L_process %*% spatial_temporal_coefficients
  eta <- as.vector(t(eta_matrix))

  structured <- as.numeric(A_full %*% eta)
  process_id <- full_template$process_id
  fine <- rnorm(length(structured), sd = truth$sigma_fine[process_id])
  measurement <- rnorm(length(structured), sd = truth$sigma_noise[process_id])

  data_full <- full_template
  data_full$mean_true <- 0
  data_full$w_true <- structured + fine
  data_full$Y_true <- structured + fine
  data_full$Z_obs <- data_full$Y_true + measurement
  data_full$Y_obs <- data_full$Z_obs

  # Sample monitoring sites independently for EACH run.
  data_full <- apply_observation_pattern(
    data_full = data_full,
    config = scaffold$config,
    n_time = n_time,
    observation_seed = run_seed + 100000L
  )

  sim <- scaffold
  sim$data_full <- data_full
  sim$data_obs <- data_full[data_full$observed, , drop = FALSE]
  sim$observation_seed <- run_seed + 100000L
  sim$precision_parameters <- precision_true$parameters
  sim$rho_process_true_approx <- truth$Sigma_process[1, 2] /
    sqrt(prod(diag(truth$Sigma_process)))
  sim$empirical_process_correlation <- cor(
    data_full$Y_true[data_full$process_id == 1L],
    data_full$Y_true[data_full$process_id == 2L]
  )
  list(simulation = sim, eta_true = eta)
}








extract_estimates <- function(fit, run_id, n_obs, seconds) {
  S <- fit$Sigma_process
  data.frame(
    run_id = run_id,
    n_obs = n_obs,
    seconds = seconds,
    converged = isTRUE(fit$converged),
    kappa_space = fit$precision_parameters$kappa_space,
    rho_time = fit$precision_parameters$rho_time,
    gamma_interaction = fit$precision_parameters$gamma_interaction,
    Sigma_11 = S[1, 1],
    Sigma_12 = S[1, 2],
    Sigma_22 = S[2, 2],
    correlation = S[1, 2] / sqrt(S[1, 1] * S[2, 2]),
    sigma_fine_1 = fit$sigma_fine[1],
    sigma_fine_2 = fit$sigma_fine[2]
  )
}

true_values <- c(
  kappa_space = truth$kappa_space,
  rho_time = truth$rho_time,
  gamma_interaction = truth$gamma_interaction,
  Sigma_11 = truth$Sigma_process[1, 1],
  Sigma_12 = truth$Sigma_process[1, 2],
  Sigma_22 = truth$Sigma_process[2, 2],
  correlation = truth$Sigma_process[1, 2] /
    sqrt(prod(diag(truth$Sigma_process))),
  sigma_fine_1 = truth$sigma_fine[1],
  sigma_fine_2 = truth$sigma_fine[2]
)

# -------------------- RUN MONTE CARLO ---------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
results <- vector("list", n_runs)
for (i in seq_len(n_runs)) {
  cat(sprintf("\nParameter recovery: run %d/%d\n", i, n_runs))
  generated <- simulate_from_model(seed_start + i - 1L)
  model <- build_model_object(
    simulation = generated$simulation,
    basis = basis,
    mean_formula = ~ 1,
    verbose = FALSE
  )
  # update_beta=FALSE means NO fitted fixed effects; no centering
  # constraint is imposed on eta. X_obs is ignored by em_fit.
  # Start AWAY from the generating parameters to avoid a trivial test.
  precision_initial <- build_precision_matrices(
    basis = basis,
    kappa_space = 0.1,
    rho_time = 0.1,
    gamma_interaction = 0.1,
    verbose = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- em_fit(
    model = model,
    precision = precision_initial,
    basis = basis,
    n_em = n_em,
    tolerance = 1e-4,
    initial_Sigma_process = diag(2),
    initial_sigma_fine = c(0.05, 0.05),
    update_beta = FALSE,
    constrain_eta_mean = FALSE,
    update_sigma_fine = TRUE,
    update_precision = TRUE,
    update_precision_every = 1L,
    precision_burnin = 1L,
    gamma_burnin = 1L,
    n_probes = n_probes,
    process_cov_shrinkage = 0,
    verbose = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  results[[i]] <- extract_estimates(fit, i, nrow(model$data_obs), elapsed)
  print(results[[i]])
  write.csv(do.call(rbind, results[seq_len(i)]),
            file.path(output_dir, "parameter_estimates.csv"), row.names = FALSE)
}

estimates <- do.call(rbind, results)
summary <- do.call(rbind, lapply(names(true_values), function(nm) {
  x <- estimates[[nm]]
  n <- length(x)
  avg <- mean(x)
  sd_mc <- if (n > 1L) sd(x) else NA_real_
  data.frame(
    parameter = nm,
    true = unname(true_values[nm]),
    mean_estimate = avg,
    bias = avg - unname(true_values[nm]),
    sd_mc = sd_mc,
    rmse = sqrt(mean((x - unname(true_values[nm]))^2)),
    ci_mean_low = if (n > 1L) avg - qt(0.975, n - 1) * sd_mc / sqrt(n) else NA_real_,
    ci_mean_high = if (n > 1L) avg + qt(0.975, n - 1) * sd_mc / sqrt(n) else NA_real_
  )
}))
write.csv(summary, file.path(output_dir, "parameter_recovery_summary.csv"), row.names = FALSE)
saveRDS(list(truth = truth, estimates = estimates, summary = summary,
             basis_dimensions = basis$dimensions),
        file.path(output_dir, "parameter_recovery_results.rds"))
cat("\nParameter recovery summary:\n")
print(summary, row.names = FALSE, digits = 4)
cat("\nSaved in:", normalizePath(output_dir), "\n")
cat("Reminder: CI for the MC mean is NOT a confidence interval for a parameter estimate.\n")



















# Profile della log-verosimiglianza completa
# rispetto a kappa, usando eta_true e gli altri parametri veri.

eta_true <- generated$eta_true
Sigma_inv <- solve(truth$Sigma_process)

kappa_grid <- seq(0.005, 0.50, length.out = 60)

profile_kappa <- sapply(kappa_grid, function(k) {
  
  Q <- build_precision_matrices(
    basis = basis,
    kappa_space = k,
    rho_time = truth$rho_time,
    gamma_interaction = truth$gamma_interaction,
    verbose = FALSE
  )$Q_ST
  
  Q0 <- Matrix::kronecker(Sigma_inv, Q)
  
  logdet <- as.numeric(
    Matrix::determinant(
      Matrix::forceSymmetric(Q0),
      logarithm = TRUE
    )$modulus
  )
  
  quadratic <- as.numeric(
    crossprod(eta_true, Q0 %*% eta_true)
  )
  
  0.5 * (logdet - quadratic)
})

plot(
  kappa_grid, profile_kappa,
  type = "l", lwd = 2,
  xlab = "kappa",
  ylab = "Complete-data log-likelihood"
)

abline(v = truth$kappa_space, col = "blue", lty = 2)
abline(v = fit$precision_parameters$kappa_space,
       col = "red", lty = 2)

cat("Maximum profile at kappa =",
    kappa_grid[which.max(profile_kappa)], "\n")

