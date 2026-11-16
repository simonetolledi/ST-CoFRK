############################################################
# 03_build_precision_matrices.R
#
# Costruzione delle matrici di precisione:
#   Q_space
#   Q_time
#   Q_ST
############################################################

build_spatial_precision <- function(
    n_basis_space_side,
    kappa_space,
    ridge = 1e-6
) {
  if (n_basis_space_side < 2L) {
    stop("n_basis_space_side deve essere almeno 2.")
  }
  if (kappa_space < 0) {
    stop("kappa_space deve essere non negativo.")
  }

  n_basis_space <- n_basis_space_side^2
  Q_space <- Matrix::Matrix(
    0,
    nrow = n_basis_space,
    ncol = n_basis_space,
    sparse = TRUE
  )

  basis_index <- matrix(
    seq_len(n_basis_space),
    nrow = n_basis_space_side,
    ncol = n_basis_space_side,
    byrow = FALSE
  )

  for (row_id in seq_len(n_basis_space_side)) {
    for (col_id in seq_len(n_basis_space_side)) {
      current <- basis_index[row_id, col_id]
      neighbors <- integer(0)

      if (row_id > 1L) {
        neighbors <- c(neighbors, basis_index[row_id - 1L, col_id])
      }
      if (row_id < n_basis_space_side) {
        neighbors <- c(neighbors, basis_index[row_id + 1L, col_id])
      }
      if (col_id > 1L) {
        neighbors <- c(neighbors, basis_index[row_id, col_id - 1L])
      }
      if (col_id < n_basis_space_side) {
        neighbors <- c(neighbors, basis_index[row_id, col_id + 1L])
      }

      Q_space[current, current] <- length(neighbors) + kappa_space^2
      Q_space[current, neighbors] <- -1
    }
  }

  Q_space <- Q_space +
    Matrix::Diagonal(n_basis_space, x = rep(ridge, n_basis_space))

  methods::as(Q_space, "dgCMatrix")
}

build_temporal_precision <- function(
    n_basis_time,
    rho_time,
    ridge = 1e-6
) {
  if (n_basis_time < 1L) {
    stop("n_basis_time deve essere almeno 1.")
  }
  
  if (rho_time < 0) {
    stop("rho_time deve essere non negativo.")
  }
  
  Q_time <- Matrix::Matrix(
    0,
    nrow = n_basis_time,
    ncol = n_basis_time,
    sparse = TRUE
  )
  
  if (n_basis_time == 1L) {
    
    Q_time[1L, 1L] <- rho_time^2
    
  } else {
    
    for (time_basis_id in seq_len(n_basis_time)) {
      
      neighbors <- integer(0)
      
      if (time_basis_id > 1L) {
        neighbors <- c(
          neighbors,
          time_basis_id - 1L
        )
      }
      
      if (time_basis_id < n_basis_time) {
        neighbors <- c(
          neighbors,
          time_basis_id + 1L
        )
      }
      
      Q_time[
        time_basis_id,
        time_basis_id
      ] <- length(neighbors) + rho_time^2
      
      Q_time[
        time_basis_id,
        neighbors
      ] <- -1
    }
  }
  
  Q_time <- Q_time +
    Matrix::Diagonal(
      n_basis_time,
      x = rep(ridge, n_basis_time)
    )
  
  methods::as(Q_time, "dgCMatrix")
}

build_st_precision <- function(
    Q_space,
    Q_time,
    gamma_interaction,
    ridge = 1e-6
) {
  if (gamma_interaction < 0) {
    stop("gamma_interaction deve essere non negativo.")
  }

  n_basis_space <- nrow(Q_space)
  n_basis_time <- nrow(Q_time)

  I_space <- Matrix::Diagonal(n_basis_space)
  I_time <- Matrix::Diagonal(n_basis_time)

  Q_ST <- Matrix::kronecker(I_time, Q_space) +
    Matrix::kronecker(Q_time, I_space) +
    gamma_interaction * Matrix::kronecker(Q_time, Q_space)

  n_basis_st <- n_basis_space * n_basis_time
  Q_ST <- Q_ST +
    Matrix::Diagonal(n_basis_st, x = rep(ridge, n_basis_st))

  methods::as(Q_ST, "dgCMatrix")
}

build_precision_matrices <- function(
    basis,
    kappa_space,
    rho_time,
    gamma_interaction,
    ridge = 1e-6,
    verbose = TRUE
) {
  n_basis_space_side <- basis$dimensions$n_basis_space_side
  n_basis_space <- basis$dimensions$n_basis_space
  n_basis_time <- basis$dimensions$n_basis_time

  if (n_basis_space != n_basis_space_side^2) {
    stop(
      "La precisione spaziale attuale richiede una griglia quadrata ",
      "di basi."
    )
  }

  Q_space <- build_spatial_precision(
    n_basis_space_side = n_basis_space_side,
    kappa_space = kappa_space,
    ridge = ridge
  )

  Q_time <- build_temporal_precision(
    n_basis_time = n_basis_time,
    rho_time = rho_time,
    ridge = ridge
  )

  Q_ST <- build_st_precision(
    Q_space = Q_space,
    Q_time = Q_time,
    gamma_interaction = gamma_interaction,
    ridge = ridge
  )

  precision <- list(
    Q_space = Q_space,
    Q_time = Q_time,
    Q_ST = Q_ST,
    parameters = list(
      kappa_space = kappa_space,
      rho_time = rho_time,
      gamma_interaction = gamma_interaction
    ),
    ridge = ridge
  )

  class(precision) <- c("cofrk_precision", "list")

  if (verbose) {
    cat("\n============================================\n")
    cat("Matrici di precisione\n")
    cat("============================================\n")
    cat("Q_space:", nrow(Q_space), "x", ncol(Q_space), "\n")
    cat("Q_time:", nrow(Q_time), "x", ncol(Q_time), "\n")
    cat("Q_ST:", nrow(Q_ST), "x", ncol(Q_ST), "\n")
    cat("kappa_space:", kappa_space, "\n")
    cat("rho_time:", rho_time, "\n")
    cat("gamma_interaction:", gamma_interaction, "\n")
  }

  precision
}
