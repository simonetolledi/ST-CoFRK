############################################################
# 04_build_model_object.R
#
# Costruzione delle matrici di progetto e dell'oggetto
# contenente i dati necessari al fit.
############################################################

build_mean_design <- function(data, mean_formula, n_process) {
  X_base <- stats::model.matrix(
    mean_formula,
    data = data
  )
  
  X <- do.call(
    cbind,
    lapply(seq_len(n_process), function(process_id) {
      X_base * as.numeric(data$process_id == process_id)
    })
  )
  
  colnames(X) <- unlist(
    lapply(seq_len(n_process), function(process_id) {
      paste0(colnames(X_base), "_Y", process_id)
    })
  )
  
  X
}

build_multivariate_design <- function(
    data,
    Phi,
    Psi,
    n_process
) {
  n_basis_space <- ncol(Phi)
  n_basis_time <- ncol(Psi)
  n_basis_st <- n_basis_space * n_basis_time
  n_data <- nrow(data)
  
  if (any(data$process_id < 1L | data$process_id > n_process)) {
    stop("process_id contiene valori non validi.")
  }
  
  # Spatial basis is compactly supported:
  # keep only spatial basis functions that are actually non-zero.
  spatial_nonzero <- lapply(
    seq_len(nrow(Phi)),
    function(space_id) which(Phi[space_id, ] != 0)
  )
  
  n_nonzero <- vapply(
    seq_len(n_data),
    function(row_id) {
      length(spatial_nonzero[[data$space_id[row_id]]]) * n_basis_time
    },
    integer(1)
  )
  
  total_nonzero <- sum(n_nonzero)
  
  row_index <- integer(total_nonzero)
  column_index <- integer(total_nonzero)
  values <- numeric(total_nonzero)
  
  position <- 1L
  
  for (row_id in seq_len(n_data)) {
    process_id <- data$process_id[row_id]
    space_id <- data$space_id[row_id]
    time_id <- data$time_id[row_id]
    
    active_space <- spatial_nonzero[[space_id]]
    n_active_space <- length(active_space)
    
    if (n_active_space == 0L) {
      next
    }
    
    # Ordering: time external, space internal.
    weights <-
      rep(
        Psi[time_id, ],
        each = length(active_space)
      ) *
      rep(
        Phi[space_id, active_space],
        times = n_basis_time
      )
    
    # Within-process ST indices corresponding only to active spatial bases.
    st_columns <-
      rep(
        (seq_len(n_basis_time) - 1L) * n_basis_space,
        each = n_active_space
      ) +
      rep(
        active_space,
        times = n_basis_time
      )
    
    process_columns <- (
      (process_id - 1L) * n_basis_st
    ) + st_columns
    
    n_entries <- length(weights)
    idx <- position:(position + n_entries - 1L)
    
    row_index[idx] <- row_id
    column_index[idx] <- process_columns
    values[idx] <- weights
    
    position <- position + n_entries
  }
  
  Matrix::sparseMatrix(
    i = row_index,
    j = column_index,
    x = values,
    dims = c(n_data, n_process * n_basis_st)
  )
}
build_model_object <- function(
    simulation,
    basis,
    mean_formula = ~ x + y + t_scaled,
    verbose = TRUE
) {
  data_full <- simulation$data_full
  data_obs <- simulation$data_obs
  n_process <- simulation$dimensions$n_process
  
  X_full <- build_mean_design(
    data = data_full,
    mean_formula = mean_formula,
    n_process = n_process
  )
  
  X_obs <- build_mean_design(
    data = data_obs,
    mean_formula = mean_formula,
    n_process = n_process
  )
  
  A_obs <- build_multivariate_design(
    data = data_obs,
    Phi = basis$Phi,
    Psi = basis$Psi,
    n_process = n_process
  )
  
  z_obs <- data_obs$Z_obs
  
  data_full$mean_hat <- 0
  data_obs$mean_hat <- 0
  
  model <- list(
    data_full = data_full,
    data_obs = data_obs,
    X_full = X_full,
    X_obs = X_obs,
    A_obs = A_obs,
    z_obs = z_obs,
    z_residual = z_obs,
    sigma_noise = simulation$sigma_noise,
    dimensions = list(
      n_process = n_process,
      n_observed = nrow(data_obs),
      n_full = nrow(data_full),
      n_beta = ncol(X_obs),
      n_basis_space = basis$dimensions$n_basis_space,
      n_basis_time = basis$dimensions$n_basis_time,
      n_basis_st = basis$dimensions$n_basis_st,
      n_eta = n_process * basis$dimensions$n_basis_st
    )
  )
  
  class(model) <- c("cofrk_model", "list")
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Oggetto modello\n")
    cat("============================================\n")
    cat("Osservazioni:", model$dimensions$n_observed, "\n")
    cat("Coefficienti beta:", model$dimensions$n_beta, "\n")
    cat("Coefficienti eta:", model$dimensions$n_eta, "\n")
    cat("Dimensione X_obs:",
        nrow(X_obs), "x", ncol(X_obs), "\n")
    cat("Dimensione A_obs:",
        nrow(A_obs), "x", ncol(A_obs), "\n")
  }
  
  model
}