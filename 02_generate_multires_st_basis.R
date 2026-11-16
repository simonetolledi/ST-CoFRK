############################################################
# 02_generate_multires_st_basis.R
#
# Costruzione delle basi spazio-temporali del modello.
#
# In questa prima versione viene generata una sola risoluzione:
#   - basi spaziali bisquare;
#   - basi temporali gaussiane.
#
# La struttura restituita è già predisposta per aggiungere
# successivamente più livelli di risoluzione.
############################################################

validate_basis_inputs <- function(
    space_grid,
    time_grid,
    n_basis_space_side,
    n_basis_time
) {
  required_space_columns <- c("x", "y", "space_id")
  required_time_columns <- c("time_id", "t_scaled")

  missing_space_columns <- setdiff(
    required_space_columns,
    names(space_grid)
  )

  missing_time_columns <- setdiff(
    required_time_columns,
    names(time_grid)
  )

  if (length(missing_space_columns) > 0L) {
    stop(
      "space_grid non contiene le colonne richieste: ",
      paste(missing_space_columns, collapse = ", ")
    )
  }

  if (length(missing_time_columns) > 0L) {
    stop(
      "time_grid non contiene le colonne richieste: ",
      paste(missing_time_columns, collapse = ", ")
    )
  }

  if (n_basis_space_side < 2L) {
    stop("n_basis_space_side deve essere almeno 2.")
  }

  if (n_basis_time < 1L) {
    stop("n_basis_time deve essere almeno 1.")
  }

  if (anyDuplicated(space_grid$space_id)) {
    stop("space_grid$space_id deve identificare univocamente le righe.")
  }

  if (anyDuplicated(time_grid$time_id)) {
    stop("time_grid$time_id deve identificare univocamente le righe.")
  }

  invisible(TRUE)
}

build_spatial_basis <- function(
    space_grid,
    n_basis_space_side,
    range_multiplier = 1.5
) {
  if (!requireNamespace("FRK", quietly = TRUE)) {
    stop(
      "Il pacchetto 'FRK' è necessario per costruire le basi spaziali."
    )
  }

  center_x <- seq(0, 1, length.out = n_basis_space_side)
  center_y <- seq(0, 1, length.out = n_basis_space_side)

  spatial_centers <- expand.grid(
    center_x = center_x,
    center_y = center_y
  )

  center_spacing <- diff(center_x)[1]
  spatial_range <- range_multiplier * center_spacing
  n_basis_space <- nrow(spatial_centers)

  spatial_basis_object <- FRK::local_basis(
    manifold = FRK::plane(),
    loc = as.matrix(
      spatial_centers[, c("center_x", "center_y")]
    ),
    scale = rep(spatial_range, n_basis_space),
    type = "bisquare"
  )

  Phi <- as.matrix(
    FRK::eval_basis(
      spatial_basis_object,
      as.matrix(space_grid[, c("x", "y")])
    )
  )

  Phi[is.na(Phi)] <- 0
  
  basis_sum <- rowSums(Phi)
  
  if (any(basis_sum <= 1e-10)) {
    stop("Alcuni BAU non sono coperti dalle basi spaziali.")
  }
  
  Phi <- sweep(
    Phi,
    1L,
    basis_sum,
    FUN = "/"
  )

  colnames(Phi) <- paste0("space_basis_", seq_len(n_basis_space))
  rownames(Phi) <- paste0("space_", space_grid$space_id)

  list(
    matrix = Phi,
    object = spatial_basis_object,
    centers = spatial_centers,
    range = spatial_range,
    center_spacing = center_spacing,
    n_basis = n_basis_space
  )
}

build_temporal_basis <- function(
    time_grid,
    n_basis_time,
    scale_time = NULL,
    normalize = TRUE
) {
  time_min <- min(time_grid$time_id)
  time_max <- max(time_grid$time_id)
  n_time <- nrow(time_grid)

  temporal_centers <- seq(
    time_min,
    time_max,
    length.out = n_basis_time
  )

  if (is.null(scale_time)) {
    scale_time <- 0.55*(n_time / n_basis_time)
  }

  Psi <- vapply(
    temporal_centers,
    FUN = function(center) {
      exp(
        -0.5 *
          ((time_grid$time_id - center) / scale_time)^2
      )
    },
    FUN.VALUE = numeric(n_time)
  )

  Psi <- as.matrix(Psi)


  Psi[is.na(Psi)] <- 0

  colnames(Psi) <- paste0("time_basis_", seq_len(n_basis_time))
  rownames(Psi) <- paste0("time_", time_grid$time_id)

  list(
    matrix = Psi,
    centers = temporal_centers,
    scale = scale_time,
    normalized = normalize,
    n_basis = n_basis_time
  )
}

generate_multires_st_basis <- function(
    simulation,
    spatial_range_multiplier = 1.5,
    temporal_scale = NULL,
    normalize_temporal = TRUE,
    verbose = TRUE
) {
  if (is.null(simulation$space_grid) ||
      is.null(simulation$time_grid) ||
      is.null(simulation$basis_settings)) {
    stop(
      "simulation deve contenere space_grid, time_grid e basis_settings."
    )
  }

  n_basis_space_side <- simulation$basis_settings$n_basis_space_side
  n_basis_time <- simulation$basis_settings$n_basis_time

  validate_basis_inputs(
    space_grid = simulation$space_grid,
    time_grid = simulation$time_grid,
    n_basis_space_side = n_basis_space_side,
    n_basis_time = n_basis_time
  )

  spatial_basis <- build_spatial_basis(
    space_grid = simulation$space_grid,
    n_basis_space_side = n_basis_space_side,
    range_multiplier = spatial_range_multiplier
  )

  temporal_basis <- build_temporal_basis(
    time_grid = simulation$time_grid,
    n_basis_time = n_basis_time,
    scale_time = temporal_scale,
    normalize = normalize_temporal
  )

  basis <- list(
    resolution_names = "level_1",
    Phi = spatial_basis$matrix,
    Psi = temporal_basis$matrix,
    spatial = spatial_basis,
    temporal = temporal_basis,
    dimensions = list(
      n_basis_space_side = n_basis_space_side,
      n_basis_space = spatial_basis$n_basis,
      n_basis_time = temporal_basis$n_basis,
      n_basis_st = spatial_basis$n_basis * temporal_basis$n_basis
    ),
    ordering = list(
      coefficient_order = "time_then_space",
      description = paste(
        "I coefficienti spazio-temporali sono ordinati per tempo;",
        "all'interno di ogni base temporale variano le basi spaziali."
      )
    )
  )

  class(basis) <- c("multires_st_basis", "list")

  if (verbose) {
    print_basis_summary(basis)
  }

  basis
}

print_basis_summary <- function(basis) {
  cat("\n============================================\n")
  cat("Basi spazio-temporali\n")
  cat("============================================\n")
  cat(
    "Basi spaziali:",
    basis$dimensions$n_basis_space,
    "(", basis$dimensions$n_basis_space_side,
    "x", basis$dimensions$n_basis_space_side, ")\n"
  )
  cat("Basi temporali:", basis$dimensions$n_basis_time, "\n")
  cat("Basi spazio-temporali:", basis$dimensions$n_basis_st, "\n")
  cat("Range spaziale:", round(basis$spatial$range, 4), "\n")
  cat("Scala temporale:", round(basis$temporal$scale, 4), "\n")

  invisible(basis)
}
