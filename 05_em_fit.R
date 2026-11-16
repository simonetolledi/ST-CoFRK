############################################################
# 05_em_fit.R
# Fit EM approssimato con Cholesky sparse e vincolo di scala esplicito.
############################################################

build_process_prior_precision <- function(Sigma_process, Q_ST) {
  Sigma_process <- (Sigma_process + t(Sigma_process)) / 2
  methods::as(
    Matrix::kronecker(solve(Sigma_process), Q_ST),
    "CsparseMatrix"
  )
}


normalize_process_covariance_scale <- function(Sigma_process, n_process) {
  Sigma_process <- (Sigma_process + t(Sigma_process)) / 2
  current_trace <- sum(diag(Sigma_process))
  if (!is.finite(current_trace) || current_trace <= 0) {
    stop("Sigma_process non ha una traccia positiva finita.")
  }
  Sigma_process * n_process / current_trace
}

solve_posterior_mode <- function(A_obs, z_residual, Q0,
                                 observation_precision) {
  weighted_A <- Matrix::Diagonal(x = observation_precision) %*% A_obs
  Q_post <- methods::as(Q0 + crossprod(A_obs, weighted_A), "CsparseMatrix")
  rhs <- crossprod(A_obs, observation_precision * z_residual)
  chol_Q_post <- Matrix::Cholesky(Q_post, LDL = FALSE, perm = TRUE)
  list(
    eta_hat = as.numeric(Matrix::solve(chol_Q_post, as.numeric(rhs))),
    Q_post = Q_post,
    chol_Q_post = chol_Q_post,
    rhs = rhs
  )
}

update_process_covariance <- function(
    eta_hat, Q_ST, n_process, n_basis_st,
    trace_probes, posterior_probe_solutions,
    ridge = 1e-8, shrinkage = 0.02) {
  
  eta_matrix <- matrix(
    eta_hat,
    nrow = n_process,
    ncol = n_basis_st,
    byrow = TRUE
  )
  
  mean_contribution <- as.matrix(
    eta_matrix %*% Q_ST %*% t(eta_matrix)
  )
  
  covariance_contribution <- matrix(
    0,
    n_process,
    n_process
  )
  
  for (m in seq_len(ncol(trace_probes))) {
    
    Z <- matrix(
      trace_probes[, m],
      nrow = n_process,
      ncol = n_basis_st,
      byrow = TRUE
    )
    
    X <- matrix(
      posterior_probe_solutions[, m],
      nrow = n_process,
      ncol = n_basis_st,
      byrow = TRUE
    )
    
    covariance_contribution <-
      covariance_contribution +
      as.matrix(
        Z %*% Q_ST %*% t(X)
      )
  }
  
  covariance_contribution <-
    covariance_contribution /
    ncol(trace_probes)
  
  Sigma_new <- as.matrix(
    (
      mean_contribution +
        covariance_contribution
    ) / n_basis_st
  )
  
  Sigma_new <- (
    Sigma_new +
      t(Sigma_new)
  ) / 2
  
  Sigma_new <-
    (1 - shrinkage) * Sigma_new +
    shrinkage * diag(
      diag(Sigma_new),
      n_process
    )
  
  min_eigenvalue <- min(
    eigen(
      Sigma_new,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  if (min_eigenvalue < ridge) {
    Sigma_new <- Sigma_new +
      diag(
        ridge - min_eigenvalue,
        n_process
      )
  }
  
  Sigma_new
}

# M-step corretta: include diag(A V_post A').
update_fine_scale_sd <- function(
    residual, A_obs, process_id, sigma_noise, n_process,
    trace_probes, posterior_probe_solutions,
    minimum_variance = 1e-8) {
  
  AZ <- as.matrix(A_obs %*% trace_probes)
  AVZ <- as.matrix(A_obs %*% posterior_probe_solutions)
  posterior_variance <- pmax(rowMeans(AZ * AVZ), 0)
  expected_squared_residual <- residual^2 + posterior_variance
  
  sigma_fine_squared <- numeric(n_process)
  for (process in seq_len(n_process)) {
    idx <- which(process_id == process)
    sigma_fine_squared[process] <- max(
      minimum_variance,
      mean(expected_squared_residual[idx]) - sigma_noise[process]^2
    )
  }
  sqrt(sigma_fine_squared)
}


update_beta_wls <- function(
    X_obs, z_obs, A_obs, eta_hat, observation_precision,
    ridge = 1e-8) {
  if (is.null(X_obs) || ncol(X_obs) == 0L) {
    return(numeric(0))
  }
  
  adjusted_response <- z_obs - as.numeric(A_obs %*% eta_hat)
  weighted_X <- Matrix::Diagonal(x = observation_precision) %*% X_obs
  lhs <- crossprod(X_obs, weighted_X)
  rhs <- crossprod(X_obs, observation_precision * adjusted_response)
  
  lhs <- methods::as(
    lhs + Matrix::Diagonal(n = ncol(X_obs), x = ridge),
    "CsparseMatrix"
  )
  
  as.numeric(Matrix::solve(lhs, rhs))
}


# Vincolo: media del campo latente sulle osservazioni = 0, per processo.
# C ha una riga per processo e agisce sui coefficienti eta.
build_latent_mean_constraints <- function(A_obs, process_id, n_process) {
  C <- matrix(0, nrow = n_process, ncol = ncol(A_obs))
  for (j in seq_len(n_process)) {
    idx <- which(process_id == j)
    if (!length(idx)) stop("Processo senza osservazioni: ", j)
    C[j, ] <- as.numeric(Matrix::colMeans(A_obs[idx, , drop = FALSE]))
  }
  if (qr(t(C))$rank < n_process) stop("Vincoli latenti linearmente dipendenti.")
  C
}

# Precalcola U = V C' e K = C V C', con V = Q_post^{-1}.
prepare_latent_constraint <- function(chol_Q_post, C) {
  U <- as.matrix(Matrix::solve(chol_Q_post, t(C), system = "A"))
  K <- C %*% U
  K <- (K + t(K)) / 2
  K_chol <- chol(K)
  list(C = C, U = U, K_chol = K_chol)
}

apply_latent_constraint <- function(B, constraint) {
  if (is.null(constraint)) return(B)
  was_vector <- is.null(dim(B))
  B <- as.matrix(B)
  adjustment <- constraint$U %*% backsolve(
    constraint$K_chol,
    forwardsolve(t(constraint$K_chol), constraint$C %*% B)
  )
  result <- B - adjustment
  if (was_vector) as.numeric(result) else result
}

update_beta_marginal_gls <- function(
    X_obs, z_obs, A_obs,
    observation_precision,
    chol_Q_post,
    ridge = 1e-8,
    latent_constraint = NULL) {
  
  if (is.null(X_obs) || ncol(X_obs) == 0L) {
    return(numeric(0))
  }
  
  # Applica V^{-1} a una matrice/vector B usando Woodbury:
  #
  # V^{-1} B =
  # R^{-1} B -
  # R^{-1} A Q_post^{-1} A' R^{-1} B
  #
  # con:
  # Q_post = Q0 + A' R^{-1} A
  
  apply_Vinv <- function(B) {
    
    B <- as.matrix(B)
    
    Rinv_B <- observation_precision * B
    
    At_Rinv_B <- crossprod(A_obs, Rinv_B)
    
    Qpost_inv_At_Rinv_B <- Matrix::solve(
      chol_Q_post,
      At_Rinv_B
    )
    
    correction <- observation_precision *
      as.matrix(A_obs %*% Qpost_inv_At_Rinv_B)
    
    result <- Rinv_B - correction
    if (!is.null(latent_constraint)) {
      # Vincolo sulla prior: correzione al GLS marginale.
      D <- observation_precision *
        as.matrix(A_obs %*% latent_constraint$U)
      result <- result + D %*% backsolve(
        latent_constraint$K_chol,
        forwardsolve(t(latent_constraint$K_chol), crossprod(D, B))
      )
    }
    result
  }
  
  # V^{-1} X
  Vinv_X <- apply_Vinv(X_obs)
  
  # V^{-1} z
  Vinv_z <- as.numeric(
    apply_Vinv(matrix(z_obs, ncol = 1L))
  )
  
  lhs <- crossprod(X_obs, Vinv_X)
  rhs <- crossprod(X_obs, Vinv_z)
  
  lhs <- as.matrix(lhs)
  diag(lhs) <- diag(lhs) + ridge
  
  as.numeric(solve(lhs, rhs))
}




relative_change <- function(new_value, old_value) {
  norm(as.matrix(new_value - old_value), type = "F") /
    max(1e-8, norm(as.matrix(old_value), type = "F"))
}

precision_convergence_vector <- function(parameters) {
  c(
    log(parameters$kappa_space),
    log(parameters$rho_time),
    log(parameters$gamma_interaction)
  )
}

em_fit <- function(
    model, precision, basis = NULL,
    n_em = 40, tolerance = 1e-4,
    initial_Sigma_process = NULL,
    initial_sigma_fine = NULL,
    initial_beta = NULL,
    beta_ridge = 1e-8,
    update_beta = FALSE,
    constrain_eta_mean = update_beta,
    update_sigma_fine = TRUE,
    update_precision = FALSE,
    update_precision_every = 2L,
    precision_burnin = 3L,
    gamma_burnin = 3L,
    update_kappa = TRUE,
    update_rho = TRUE,
    update_gamma = TRUE,
    kappa_bounds = c(1e-4, 3),
    rho_bounds = c(1e-4, 3),
    gamma_bounds = c(1e-6, 10),
    n_probes = 32L,
    precision_penalty = NULL,
    precision_damping = 0.5,
    process_cov_shrinkage = 0.1,
    verbose = TRUE) {
  
  n_process <- model$dimensions$n_process
  n_basis_st <- model$dimensions$n_basis_st
  if (is.null(initial_Sigma_process)) initial_Sigma_process <- diag(n_process)
  if (is.null(initial_sigma_fine)) initial_sigma_fine <- rep(0.10, n_process)
  if (update_precision && is.null(basis)) stop("basis necessario.")
  
  Sigma_process <- normalize_process_covariance_scale(
    initial_Sigma_process, n_process
  )
  sigma_fine <- initial_sigma_fine
  
  has_fixed_effects <- update_beta && !is.null(model$X_obs) && ncol(model$X_obs) > 0L
  if (constrain_eta_mean && !has_fixed_effects) {
    stop("constrain_eta_mean richiede update_beta = TRUE e X_obs non vuota.")
  }
  latent_C <- if (constrain_eta_mean) {
    build_latent_mean_constraints(
      model$A_obs, model$data_obs$process_id, n_process
    )
  } else NULL
  if (has_fixed_effects) {
    if (is.null(model$z_obs)) {
      stop("model$z_obs necessario quando model$X_obs e' presente.")
    }
    if (nrow(model$X_obs) != length(model$z_obs)) {
      stop("Numero di righe di model$X_obs diverso dalla lunghezza di model$z_obs.")
    }
    if (is.null(initial_beta)) initial_beta <- rep(0, ncol(model$X_obs))
    if (length(initial_beta) != ncol(model$X_obs)) {
      stop("initial_beta deve avere lunghezza pari a ncol(model$X_obs).")
    }
    beta <- as.numeric(initial_beta)
  } else {
    beta <- numeric(0)
  }
  
  precision_parameters <- precision$parameters
  numerical_ridge <- if (is.null(precision$ridge)) 1e-8 else precision$ridge
  
  # Q_ST mantiene la propria scala; l'identificabilità globale
  # è imposta da trace(Sigma_process) = n_process.
  Q_ST <- precision$Q_ST
  precision$Q_ST <- Q_ST
  precision$Q_ST_scale_factor <- 1
  
  # Sonde comuni a tutte le iterazioni: meno rumore nel criterio EM.
  trace_probes <- generate_rademacher_probes(
    n_process * n_basis_st, n_probes = n_probes, seed = 1000L
  )
  history <- vector("list", n_em)
  
  for (em_iteration in seq_len(n_em)) {
    noise_variance <- model$sigma_noise[model$data_obs$process_id]^2 +
      sigma_fine[model$data_obs$process_id]^2
    observation_precision <- 1 / noise_variance
    
    z_residual <- if (has_fixed_effects) {
      model$z_obs - as.numeric(model$X_obs %*% beta)
    } else {
      model$z_residual
    }
    
    Q0 <- build_process_prior_precision(Sigma_process, Q_ST)
    
    # Prima soluzione del posterior: costruisce e fattorizza Q_post.
    posterior <- solve_posterior_mode(
      model$A_obs, z_residual, Q0, observation_precision
    )
    
    latent_constraint <- if (constrain_eta_mean) {
      prepare_latent_constraint(posterior$chol_Q_post, latent_C)
    } else NULL
    eta_hat <- apply_latent_constraint(posterior$eta_hat, latent_constraint)
    
    beta_new <- if (has_fixed_effects) {
      update_beta_marginal_gls(
        X_obs = model$X_obs,
        z_obs = model$z_obs,
        A_obs = model$A_obs,
        observation_precision = observation_precision,
        chol_Q_post = posterior$chol_Q_post,
        ridge = beta_ridge,
        latent_constraint = latent_constraint
      )
    } else {
      beta
    }
    
    if (has_fixed_effects) {
      
      # Beta cambia il residual, ma NON Q_post.
      z_residual <- model$z_obs -
        as.numeric(model$X_obs %*% beta_new)
      
      # Cambia soltanto il RHS.
      rhs_new <- crossprod(
        model$A_obs,
        observation_precision * z_residual
      )
      
      # Riutilizza la Cholesky già calcolata.
      eta_hat <- as.numeric(
        Matrix::solve(
          posterior$chol_Q_post,
          as.numeric(rhs_new)
        )
      )
      
      eta_hat <- apply_latent_constraint(eta_hat, latent_constraint)
      
      # Aggiorniamo eta_hat anche dentro posterior perché il codice
      # successivo utilizza l'oggetto posterior.
      posterior$eta_hat <- eta_hat
      posterior$rhs <- rhs_new
    }
    
    posterior_probe_solutions <- solve_posterior_probes(
      posterior$chol_Q_post, trace_probes
    )
    posterior_probe_solutions <- apply_latent_constraint(
      posterior_probe_solutions, latent_constraint
    )
    
    Sigma_new <- update_process_covariance(
      eta_hat, Q_ST, n_process, n_basis_st,
      trace_probes, posterior_probe_solutions,
      ridge = 1e-8,
      shrinkage = process_cov_shrinkage
    )
    Sigma_new <- normalize_process_covariance_scale(
      Sigma_new, n_process
    )
    
    residual_obs <- z_residual -
      as.numeric(model$A_obs %*% eta_hat)
    
    sigma_fine_new <- if (update_sigma_fine) {
      
      # Update EM non penalizzato
      sigma_fine_mle <- update_fine_scale_sd(
        residual_obs, model$A_obs, model$data_obs$process_id,
        model$sigma_noise, n_process,
        trace_probes, posterior_probe_solutions
      )
      
      # Shrinkage della VARIANZA verso sigma_fine_target^2
      lambda_fine <- 0
      sigma_fine_target <- 0
      
      sigma_fine_new <- sqrt(
        (1 - lambda_fine) * sigma_fine_mle^2 +
          lambda_fine * sigma_fine_target^2
      )
      
      sigma_fine_new
      
    } else {
      sigma_fine
    }
    
    precision_parameters_old <- precision_parameters
    precision_updated <- update_precision &&
      em_iteration > precision_burnin &&
      em_iteration %% update_precision_every == 0L
    
    if (precision_updated) {
      precision_parameters <- update_precision_parameters(
        current_parameters = precision_parameters,
        eta_hat = eta_hat,
        Sigma_process = Sigma_process,
        basis = basis,
        trace_probes = trace_probes,
        posterior_probe_solutions = posterior_probe_solutions,
        update_kappa = update_kappa,
        update_rho = update_rho,
        update_gamma = update_gamma && em_iteration > gamma_burnin,
        kappa_bounds = kappa_bounds,
        rho_bounds = rho_bounds,
        gamma_bounds = gamma_bounds,
        ridge = numerical_ridge,
        penalty = precision_penalty,
        damping = precision_damping,
        verbose = verbose
      )
      
      # Usa esattamente la Q_ST già valutata e accettata.
      Q_candidate <- attr(precision_parameters, "Q_ST_evaluated")
      if (!is.null(Q_candidate)) {
        Q_ST <- Q_candidate
      } else {
        rebuilt <- build_precision_matrices(
          basis = basis,
          kappa_space = precision_parameters$kappa_space,
          rho_time = precision_parameters$rho_time,
          gamma_interaction = precision_parameters$gamma_interaction,
          ridge = numerical_ridge,
          verbose = FALSE
        )
        Q_ST <- rebuilt$Q_ST
      }
      precision$Q_ST <- Q_ST
      precision$parameters <- precision_parameters
    }
    
    change_beta <- if (has_fixed_effects) {
      relative_change(beta_new, beta)
    } else 0
    change_Sigma <- relative_change(Sigma_new, Sigma_process)
    change_fine <- relative_change(sigma_fine_new, sigma_fine)
    change_precision <- if (precision_updated) {
      relative_change(
        precision_convergence_vector(precision_parameters),
        precision_convergence_vector(precision_parameters_old)
      )
    } else 0
    change_em <- max(change_beta, change_Sigma, change_fine, change_precision)
    
    history[[em_iteration]] <- list(
      iteration = em_iteration,
      beta = beta_new,
      Sigma_process = Sigma_new,
      sigma_fine = sigma_fine_new,
      precision_parameters = precision_parameters,
      precision_update_accepted =
        attr(precision_parameters, "precision_update_accepted"),
      precision_objective =
        attr(precision_parameters, "expected_precision_objective"),
      precision_logdet_term =
        attr(precision_parameters, "precision_logdet_term"),
      precision_mean_term =
        attr(precision_parameters, "precision_mean_term"),
      precision_trace_term =
        attr(precision_parameters, "precision_trace_term"),
      precision_at_boundary =
        attr(precision_parameters, "precision_at_boundary"),
      gamma_updated = precision_updated && update_gamma &&
        em_iteration > gamma_burnin,
      change_beta = change_beta,
      change_Sigma = change_Sigma,
      change_fine = change_fine,
      change_precision = change_precision,
      change_em = change_em
    )
    
    if (verbose) {
      cat("\nEM iteration", em_iteration, "\n")
      cat("Sigma_process:\n"); print(round(Sigma_new, 4))
      if (has_fixed_effects) {
        cat("beta:", round(beta_new, 4), "\n")
      }
      if (constrain_eta_mean) {
        cat("latent means:", signif(as.numeric(latent_C %*% eta_hat), 4), "\n")
      }
      cat("sigma_fine:", round(sigma_fine_new, 4), "\n")
      cat("relative change:", signif(change_em, 4),
          "| beta:", signif(change_beta, 4),
          "| Sigma:", signif(change_Sigma, 4),
          "| fine:", signif(change_fine, 4),
          "| precision:", signif(change_precision, 4), "\n")
    }
    
    beta <- beta_new
    Sigma_process <- Sigma_new
    sigma_fine <- sigma_fine_new
    if (change_em < tolerance && em_iteration > precision_burnin) {
      history <- history[seq_len(em_iteration)]
      break
    }
  }
  
  noise_variance <- model$sigma_noise[model$data_obs$process_id]^2 +
    sigma_fine[model$data_obs$process_id]^2
  observation_precision <- 1 / noise_variance
  z_residual <- if (has_fixed_effects) {
    model$z_obs - as.numeric(model$X_obs %*% beta)
  } else {
    model$z_residual
  }
  Q0 <- build_process_prior_precision(Sigma_process, Q_ST)
  posterior <- solve_posterior_mode(
    model$A_obs, z_residual, Q0, observation_precision
  )
  final_constraint <- if (constrain_eta_mean) {
    prepare_latent_constraint(posterior$chol_Q_post, latent_C)
  } else NULL
  posterior$eta_hat <- apply_latent_constraint(
    posterior$eta_hat, final_constraint
  )
  
  rho_process_hat <- if (n_process >= 2L) {
    Sigma_process[1L, 2L] /
      sqrt(Sigma_process[1L, 1L] * Sigma_process[2L, 2L])
  } else NA_real_
  
  fit <- list(
    beta_hat = beta,
    eta_hat = posterior$eta_hat,
    constrain_eta_mean = constrain_eta_mean,
    latent_mean_constraint = latent_C,
    Sigma_process = Sigma_process,
    sigma_fine = sigma_fine,
    rho_process_hat = rho_process_hat,
    Q_ST = Q_ST,
    Q0 = Q0,
    Q_post = posterior$Q_post,
    chol_Q_post = posterior$chol_Q_post,
    observation_precision = observation_precision,
    precision_parameters = precision_parameters,
    history = history,
    converged = length(history) > 0L &&
      tail(history, 1L)[[1L]]$change_em < tolerance
  )
  class(fit) <- c("cofrk_fit", "list")
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Fit EM completato\n")
    cat("============================================\n")
    cat("rho_process_hat:", round(rho_process_hat, 4), "\n")
    if (has_fixed_effects) {
      cat("beta_hat:", round(beta, 4), "\n")
    }
    cat("sigma_fine:", round(sigma_fine, 4), "\n")
    
    cat(
      "precision parameters:",
      "kappa =", round(
        precision_parameters$kappa_space,
        4
      ),
      "| rho =", round(
        precision_parameters$rho_time,
        4
      ),
      "| gamma =", round(
        precision_parameters$gamma_interaction,
        6
      ),
      "\n"
    )
    
    cat("converged:", fit$converged, "\n")
  }
  fit
}
