############################################################
# 08_fit_FRK_benchmark.R
#
# Benchmark FRK spazio-temporale univariato ufficiale.
#
# Modifica minima:
# - il FIT usa sempre e solo le osservazioni della run;
# - auto_BAUs() costruisce invece il dominio usando l'intera
#   griglia spazio-temporale della simulazione, così il supporto
#   FRK non cambia tra una run e l'altra.
############################################################

find_first_column <- function(column_names, patterns) {
  for (pattern in patterns) {
    matches <- grep(
      pattern,
      column_names,
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(matches) > 0L) {
      return(matches[1L])
    }
  }
  NA_character_
}

fit_FRK_benchmark <- function(
    simulation,
    basis,
    process_ids = 2L,
    n_em = 50,
    tolerance = 0.01,
    verbose = TRUE
) {
  data_full <- simulation$data_full
  data_obs <- simulation$data_obs
  space_grid <- simulation$space_grid
  n_time <- simulation$dimensions$n_time
  n_space <- simulation$dimensions$n_space
  
  data_full$Y_hat_FRK <- NA_real_
  fitted_models <- vector("list", length(process_ids))
  names(fitted_models) <- paste0("process_", process_ids)
  
  for (benchmark_index in seq_along(process_ids)) {
    process_id <- process_ids[benchmark_index]
    
    if (verbose) {
      cat("\nFit FRK benchmark per processo", process_id, "\n")
    }
    
    obs_process <- data_obs[
      data_obs$process_id == process_id,
      ,
      drop = FALSE
    ]
    
    df_obs_frk <- data.frame(
      x = obs_process$x,
      y = obs_process$y,
      time = as.Date("2000-01-01") +
        (obs_process$time_id - 1L),
      z = obs_process$Z_obs,
      std = rep(
        simulation$sigma_noise[process_id],
        nrow(obs_process)
      )
    )
    
    ST_object <- spacetime::stConstruct(
      x = df_obs_frk,
      space = c("x", "y"),
      time = "time",
      crs = sp::CRS(NA_character_),
      interval = TRUE
    )
    
    ST_object@data <- data.frame(
      z = df_obs_frk$z,
      std = df_obs_frk$std
    )
    
    ########################################################
    # MODIFICA MINIMA:
    # supporto geometrico completo e fisso tra le run
    ########################################################
    
    support_process <- data_full[
      data_full$process_id == process_id,
      ,
      drop = FALSE
    ]
    
    df_support <- data.frame(
      x = support_process$x,
      y = support_process$y,
      time = as.Date("2000-01-01") +
        (support_process$time_id - 1L)
    )
    
    ST_support <- spacetime::stConstruct(
      x = df_support,
      space = c("x", "y"),
      time = "time",
      crs = sp::CRS(NA_character_),
      interval = TRUE
    )
    
    grid_BAUs <- FRK::auto_BAUs(
      manifold = FRK::STplane(),
      data = ST_support,
      cellsize = c(
        diff(simulation$x_values)[1L],
        diff(simulation$y_values)[1L],
        1
      ),
      type = "grid",
      convex = -0.05,
      tunit = "days",
      nonconvex_hull = FALSE
    )
    
    grid_BAUs$fs <- 1
    
    bau_data <- as(grid_BAUs, "data.frame")
    
    x_column <- find_first_column(
      names(bau_data),
      c("^x$", "coords.x1", "^s1$")
    )
    y_column <- find_first_column(
      names(bau_data),
      c("^y$", "coords.x2", "^s2$")
    )
    time_column <- find_first_column(
      names(bau_data),
      c("time", "date")
    )
    
    if (is.na(x_column) || is.na(y_column)) {
      stop("Non trovo le coordinate spaziali nelle BAU FRK.")
    }
    
    grid_BAUs$cov_x <- bau_data[[x_column]]
    grid_BAUs$cov_y <- bau_data[[y_column]]
    
    if (!is.na(time_column)) {
      grid_BAUs$cov_t <- as.numeric(
        as.Date(bau_data[[time_column]]) -
          as.Date("2000-01-01")
      ) / (n_time - 1L)
    } else {
      grid_BAUs$cov_t <- rep(
        seq(0, 1, length.out = n_time),
        each = n_space
      )[seq_len(nrow(bau_data))]
    }
    
    spatial_basis_frk <- basis$spatial$object
    
    temporal_basis_frk <- FRK::local_basis(
      manifold = FRK::real_line(),
      type = "Gaussian",
      loc = matrix(
        basis$temporal$centers,
        ncol = 1L
      ),
      scale = rep(
        basis$temporal$scale,
        basis$dimensions$n_basis_time
      )
    )
    
    st_basis_frk <- FRK::TensorP(
      spatial_basis_frk,
      temporal_basis_frk
    )
    
    S_frk <- FRK::SRE(
      f = z ~ 1,
      data = list(ST_object),
      basis = st_basis_frk,
      BAUs = grid_BAUs,
      est_error = FALSE,
      average_in_BAU = TRUE
    )
    
    S_frk <- FRK::SRE.fit(
      S_frk,
      n_EM = n_em,
      tol = tolerance,
      print_lik = FALSE
    )
    
    prediction_grid <- methods::selectMethod(
      "predict",
      signature = "SRE"
    )(
      S_frk,
      obs_fs = FALSE
    )
    
    prediction_data <- as(
      prediction_grid,
      "data.frame"
    )
    
    prediction_column <- find_first_column(
      names(prediction_grid@data),
      c("^mu$", "pred", "mean", "fit")
    )
    
    if (is.na(prediction_column)) {
      stop("Non trovo la colonna di predizione FRK.")
    }
    
    x_prediction_column <- find_first_column(
      names(prediction_data),
      c("^x$", "coords.x1", "^s1$")
    )
    y_prediction_column <- find_first_column(
      names(prediction_data),
      c("^y$", "coords.x2", "^s2$")
    )
    time_prediction_column <- find_first_column(
      names(prediction_data),
      c("time", "date")
    )
    
    if (
      is.na(x_prediction_column) ||
      is.na(y_prediction_column)
    ) {
      stop("Non trovo le coordinate nelle predizioni FRK.")
    }
    
    if (!is.na(time_prediction_column)) {
      prediction_time <- as.Date(
        prediction_data[[time_prediction_column]]
      )
      prediction_time_id <- as.numeric(
        prediction_time - as.Date("2000-01-01")
      ) + 1L
    } else {
      rounded_location <- round(
        prediction_data[[x_prediction_column]],
        10
      ) + 1000 * round(
        prediction_data[[y_prediction_column]],
        10
      )
      n_space_prediction <- length(unique(rounded_location))
      prediction_time_id <- rep(
        seq_len(
          ceiling(
            nrow(prediction_data) /
              n_space_prediction
          )
        ),
        each = n_space_prediction
      )[seq_len(nrow(prediction_data))]
    }
    
    prediction_space_id <- apply(
      prediction_data[
        ,
        c(x_prediction_column, y_prediction_column),
        drop = FALSE
      ],
      1L,
      function(coordinates) {
        which.min(
          (space_grid$x - coordinates[1L])^2 +
            (space_grid$y - coordinates[2L])^2
        )
      }
    )
    
    prediction_key <- data.frame(
      space_id = prediction_space_id,
      time_id = prediction_time_id,
      prediction = prediction_data[[prediction_column]]
    )
    
    prediction_key <- prediction_key[
      prediction_key$time_id >= 1L &
        prediction_key$time_id <= n_time,
      ,
      drop = FALSE
    ]
    
    prediction_key <- stats::aggregate(
      prediction ~ space_id + time_id,
      data = prediction_key,
      FUN = mean
    )
    
    idx_full_process <- which(
      data_full$process_id == process_id
    )
    
    full_key <- data_full[
      idx_full_process,
      c("space_id", "time_id"),
      drop = FALSE
    ]
    full_key$row_id <- idx_full_process
    
    merged_prediction <- merge(
      full_key,
      prediction_key,
      by = c("space_id", "time_id"),
      all.x = TRUE,
      sort = FALSE
    )
    
    data_full$Y_hat_FRK[
      merged_prediction$row_id
    ] <- merged_prediction$prediction
    
    fitted_models[[benchmark_index]] <- S_frk
    
    if (verbose) {
      cat(
        "Predizioni assegnate:",
        sum(!is.na(data_full$Y_hat_FRK[idx_full_process])),
        "su", length(idx_full_process), "\n"
      )
    }
  }
  
  list(
    data_full = data_full,
    fitted_models = fitted_models,
    process_ids = process_ids
  )
}
