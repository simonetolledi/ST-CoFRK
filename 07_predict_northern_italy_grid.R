############################################################
# 07_predict_northern_italy_grid.R
#
# Prediction del fit COFRK GIA' STIMATO su una griglia completa
# di BAU da 1.6 km, ritagliata sul vero confine del Nord Italia.
#
# IMPORTANTE:
# - NON rifitta il modello;
# - usa fit$eta_hat gia' stimato;
# - usa ESATTAMENTE le stesse basi spaziali del fit
#   (basis$spatial$object);
# - usa le stesse basi temporali Gaussiane (basis$Psi);
# - il Nord Italia e' definito tramite NUTS-1 Eurostat:
#       ITC = Nord-Ovest
#       ITH = Nord-Est
#
# Dipendenze:
#   sf, Matrix, FRK, giscoR
#
# Una tantum, se serve:
# install.packages("giscoR")
############################################################


get_northern_italy_boundary <- function(
    year = "2021",
    epsg = 3035
) {
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Serve il pacchetto 'sf'.")
  }
  
  if (!requireNamespace("giscoR", quietly = TRUE)) {
    stop(
      "Serve il pacchetto 'giscoR'. Installalo una volta con:\n",
      "install.packages('giscoR')"
    )
  }
  
  italy_nuts1 <- giscoR::gisco_get_nuts(
    year = year,
    epsg = epsg,
    nuts_level = 1,
    country = "IT",
    resolution = "03"
  )
  
  north <- italy_nuts1[
    italy_nuts1$NUTS_ID %in% c("ITC", "ITH"),
  ]
  
  if (nrow(north) == 0L) {
    stop("Non trovo ITC/ITH nei confini NUTS Eurostat.")
  }
  
  north <- sf::st_make_valid(north)
  
  # Un unico poligono per il Nord Italia.
  north_union <- sf::st_union(north)
  
  sf::st_sf(
    region = "Northern Italy",
    geometry = north_union
  )
}


build_northern_italy_prediction_grid <- function(
    real_data,
    cellsize_m = NULL,
    nuts_year = "2021",
    epsg = NULL
) {
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Serve il pacchetto 'sf'.")
  }
  
  if (is.null(cellsize_m)) {
    cellsize_m <- real_data$settings$spatial_bau_size_m
  }
  
  if (is.null(cellsize_m)) {
    cellsize_m <- 1600
  }
  
  if (is.null(epsg)) {
    epsg <- real_data$settings$projected_crs
  }
  
  if (is.null(epsg)) {
    epsg <- 3035
  }
  
  north <- get_northern_italy_boundary(
    year = nuts_year,
    epsg = epsg
  )
  
  ##########################################################
  # Griglia regolare da 1.6 km
  #
  # what = "centers": generiamo i centroidi della griglia.
  # Conserviamo solo quelli che cadono nel Nord Italia.
  ##########################################################
  
  bbox <- sf::st_bbox(north)
  
  # ----------------------------------------------------------
  # Griglia 1.6 km allineata ESATTAMENTE ai BAU del fit
  # ----------------------------------------------------------
  
  old_space <- real_data$space_grid
  
  x_ref <- min(old_space$bau_easting, na.rm = TRUE)
  y_ref <- min(old_space$bau_northing, na.rm = TRUE)
  
  # Estendiamo la stessa lattice del fit fino a coprire il Nord Italia
  kx_min <- floor((bbox["xmin"] - x_ref) / cellsize_m)
  kx_max <- ceiling((bbox["xmax"] - x_ref) / cellsize_m)
  
  ky_min <- floor((bbox["ymin"] - y_ref) / cellsize_m)
  ky_max <- ceiling((bbox["ymax"] - y_ref) / cellsize_m)
  
  x_seq <- x_ref + (kx_min:kx_max) * cellsize_m
  y_seq <- y_ref + (ky_min:ky_max) * cellsize_m
  
  xy_grid <- expand.grid(
    bau_easting = x_seq,
    bau_northing = y_seq
  )
  
  grid_centers <- sf::st_as_sf(
    xy_grid,
    coords = c("bau_easting", "bau_northing"),
    crs = epsg,
    remove = FALSE
  )
  
  inside <- lengths(
    sf::st_intersects(grid_centers, north)
  ) > 0L
  
  grid_centers <- grid_centers[inside, ]
  
  xy <- sf::st_coordinates(grid_centers)
  
  prediction_grid <- data.frame(
    grid_id = seq_len(nrow(grid_centers)),
    bau_easting = xy[, 1L],
    bau_northing = xy[, 2L]
  )
  
  ##########################################################
  # Coordinate normalizzate DEL FIT
  #
  # Il fit originale ha costruito le basi su x,y in [0,1].
  # Per NON cambiare il significato di eta_hat dobbiamo usare
  # la stessa trasformazione con cui erano state create le
  # coordinate del vecchio space_grid.
  ##########################################################
  
  old_space <- real_data$space_grid
  
  required <- c(
    "x", "y",
    "bau_easting", "bau_northing"
  )
  
  missing <- setdiff(required, names(old_space))
  
  if (length(missing) > 0L) {
    stop(
      "real_data$space_grid non contiene: ",
      paste(missing, collapse = ", "),
      ". Servono per ricostruire la normalizzazione usata nel fit."
    )
  }
  
  # Ricaviamo la trasformazione affine direttamente dai dati
  # che sono stati realmente usati nel fit.
  fit_x_min <- min(old_space$bau_easting, na.rm = TRUE)
  fit_x_max <- max(old_space$bau_easting, na.rm = TRUE)
  fit_y_min <- min(old_space$bau_northing, na.rm = TRUE)
  fit_y_max <- max(old_space$bau_northing, na.rm = TRUE)
  
  if (fit_x_max <= fit_x_min || fit_y_max <= fit_y_min) {
    stop("Range spaziale del fit non valido.")
  }
  
  prediction_grid$x <- (
    prediction_grid$bau_easting - fit_x_min
  ) / (fit_x_max - fit_x_min)
  
  prediction_grid$y <- (
    prediction_grid$bau_northing - fit_y_min
  ) / (fit_y_max - fit_y_min)
  
  ##########################################################
  # Diagnostica extrapolation
  ##########################################################
  
  prediction_grid$inside_fit_bbox <-
    prediction_grid$x >= 0 &
    prediction_grid$x <= 1 &
    prediction_grid$y >= 0 &
    prediction_grid$y <= 1
  
  
  
  ##########################################################
  # Population density - Eurostat Census Grid
  ##########################################################
  
  eurostat_gpkg <- file.path(
    "eurostat_tmp",
    "Eurostat_Census-GRID_2021_V3",
    "ESTAT_Census_2021_V3.gpkg"
  )
  
  if (!file.exists(eurostat_gpkg)) {
    stop("GeoPackage Eurostat non trovato: ", eurostat_gpkg)
  }
  
  prediction_sf <- sf::st_as_sf(
    prediction_grid,
    coords = c("bau_easting", "bau_northing"),
    crs = epsg,
    remove = FALSE
  )
  
  population_grid <- sf::st_read(
    eurostat_gpkg,
    layer = "ESTAT_Census_2021_V3",
    wkt_filter = sf::st_as_text(
      sf::st_as_sfc(sf::st_bbox(prediction_sf))
    ),
    quiet = TRUE
  )
  
  population_grid <- population_grid[, "T", drop = FALSE]
  population_grid$T <- as.numeric(population_grid$T)
  
  # Per dati puntuali vogliamo UNA sola cella Eurostat per BAU.
  # st_intersects può restituire più celle quando il punto cade
  # esattamente su un bordo della griglia Eurostat.
  
  pop_idx <- sf::st_intersects(
    prediction_sf,
    population_grid
  )
  
  population_density <- vapply(
    pop_idx,
    function(idx) {
      if (length(idx) == 0L) {
        return(NA_real_)
      }
      
      # Se il punto cade sul bordo di più celle, usa la prima.
      # I valori duplicati derivano solo dall'ambiguità geometrica
      # del confine, non da BAU differenti.
      population_grid$T[idx[1L]]
    },
    numeric(1)
  )
  
  prediction_grid$population_density <- population_density
  
  prediction_grid$population_density[
    is.na(prediction_grid$population_density)
  ] <- 0
  prediction_grid$population_density[
    is.na(prediction_grid$population_density)
  ] <- 0
  
  prediction_grid$population_density_scaled <-
    (
      prediction_grid$population_density -
        real_data$covariate_scaling$pop_mean
    ) /
    real_data$covariate_scaling$pop_sd
  
  # -------------------------------------------------------------------------
  # Altitude from Copernicus DEM
  # -------------------------------------------------------------------------
  
  dem_path <- file.path(
    "data",
    "copernicus_dem_90_northern_italy.tiff"
  )
  
  if (!file.exists(dem_path)) {
    stop("Copernicus DEM not found: ", dem_path)
  }
  
  dem <- terra::rast(dem_path)
  
  prediction_points <- sf::st_as_sf(
    prediction_grid,
    coords = c("bau_easting", "bau_northing"),
    crs = epsg,
    remove = FALSE
  )
  
  prediction_points_dem <- sf::st_transform(
    prediction_points,
    terra::crs(dem)
  )
  
  dem_values <- terra::extract(
    dem,
    terra::vect(prediction_points_dem)
  )
  
  prediction_grid$altitude <- dem_values[[2]]
  
  prediction_grid$log_altitude <-
    log1p(pmax(prediction_grid$altitude, 0))
  
  prediction_grid$log_altitude_scaled <-
    (
      prediction_grid$log_altitude -
        real_data$covariate_scaling$logalt_mean
    ) /
    real_data$covariate_scaling$logalt_sd
  
  rm(dem, prediction_points, prediction_points_dem, dem_values)
  gc()
  
  
  
  
  attr(prediction_grid, "north_boundary") <- north
  attr(prediction_grid, "cellsize_m") <- cellsize_m
  attr(prediction_grid, "fit_bounds") <- list(
    easting_min = fit_x_min,
    easting_max = fit_x_max,
    northing_min = fit_y_min,
    northing_max = fit_y_max
  )
  
  prediction_grid
}


evaluate_fitted_spatial_basis_on_grid <- function(
    prediction_grid,
    basis
) {
  
  if (!requireNamespace("FRK", quietly = TRUE)) {
    stop("Serve il pacchetto 'FRK'.")
  }
  
  if (is.null(basis$spatial$object)) {
    stop(
      "basis$spatial$object non trovato. ",
      "Serve l'oggetto delle basi usato nel fit."
    )
  }
  
  Phi_new <- as.matrix(
    FRK::eval_basis(
      basis$spatial$object,
      as.matrix(
        prediction_grid[, c("x", "y")]
      )
    )
  )
  
  Phi_new[is.na(Phi_new)] <- 0
  
  basis_sum <- rowSums(Phi_new)
  
  valid <- basis_sum > 1e-10
  
  Phi_new[valid, ] <- sweep(
    Phi_new[valid, , drop = FALSE],
    1L,
    basis_sum[valid],
    FUN = "/"
  )
  
  colnames(Phi_new) <- colnames(basis$Phi)
  
  Phi_new
}


predict_cofrk_northern_italy <- function(
    real_data,
    fit,
    basis,
    time_ids = NULL,
    pollutants = NULL,
    cellsize_m = NULL,
    nuts_year = "2021",
    epsg = NULL,
    drop_zero_basis = TRUE,
    include_fine_scale_at_observed = TRUE,
    verbose = TRUE
) {
  
  ##########################################################
  # 1. Griglia geografica completa
  ##########################################################
  
  grid <- build_northern_italy_prediction_grid(
    real_data = real_data,
    cellsize_m = cellsize_m,
    nuts_year = nuts_year,
    epsg = epsg
  )
  
  north_boundary <- attr(grid, "north_boundary")
  
  ##########################################################
  # 2. Stesse basi spaziali del fit
  ##########################################################
  
  Phi_new <- evaluate_fitted_spatial_basis_on_grid(
    prediction_grid = grid,
    basis = basis
  )
  
  basis_coverage <- rowSums(abs(Phi_new)) > 0
  grid$basis_covered <- basis_coverage
  
  if (drop_zero_basis) {
    grid <- grid[basis_coverage, , drop = FALSE]
    Phi_new <- Phi_new[basis_coverage, , drop = FALSE]
  }
  
  grid$space_id <- seq_len(nrow(grid))
  
  ##########################################################
  # 3. Tempi
  ##########################################################
  
  time_grid <- real_data$time_grid
  
  if (is.null(time_ids)) {
    idx <- unique(
      as.integer(
        round(
          seq(
            1,
            nrow(time_grid),
            length.out = min(4L, nrow(time_grid))
          )
        )
      )
    )
    time_ids <- time_grid$time_id[idx]
  }
  
  if (!all(time_ids %in% time_grid$time_id)) {
    stop("Almeno un time_id richiesto non esiste nel fit.")
  }
  
  ##########################################################
  # 4. Processi
  ##########################################################
  
  process_lookup <- real_data$process_lookup
  
  if (is.null(pollutants)) {
    pollutants <- names(process_lookup)
  }
  
  if (!all(pollutants %in% names(process_lookup))) {
    stop(
      "Pollutant non riconosciuto. Disponibili: ",
      paste(names(process_lookup), collapse = ", ")
    )
  }
  
  process_ids <- unname(process_lookup[pollutants])
  
  ##########################################################
  # 5. Prediction SENZA costruire A_new
  #
  # build_multivariate_design() usa, per ogni riga:
  #
  #   as.vector(outer(Psi[time_id, ], Phi[space_id, ]))
  #
  # Qui calcoliamo ESATTAMENTE lo stesso prodotto scalare,
  # ma sfruttiamo la struttura tensoriale:
  #
  #   eta_process -> matrice M x R
  #   Psi_t %*% eta_matrix -> coefficienti spaziali al tempo t
  #   Phi_new %*% coeff_spatial -> prediction su tutti i BAU
  #
  # Quindi non viene mai allocata la gigantesca A_new.
  ##########################################################
  
  n_basis_space <- ncol(Phi_new)
  n_basis_time <- ncol(basis$Psi)
  n_basis_st <- n_basis_space * n_basis_time
  n_process <- real_data$dimensions$n_process
  
  expected_eta_length <- n_process * n_basis_st
  
  if (length(fit$eta_hat) != expected_eta_length) {
    stop(
      "Dimensione di fit$eta_hat incompatibile con le basi: ",
      length(fit$eta_hat), " invece di ", expected_eta_length, "."
    )
  }
  
  prediction_parts <- vector(
    "list",
    length(process_ids) * length(time_ids)
  )
  
  part_id <- 0L
  
  for (process_id in process_ids) {
    
    pollutant_name <- names(process_lookup)[
      match(process_id, process_lookup)
    ]
    
    eta_start <- (process_id - 1L) * n_basis_st + 1L
    eta_end <- process_id * n_basis_st
    
    eta_process <- fit$eta_hat[eta_start:eta_end]
    
    # Ordinamento coerente con Q_ST e A_obs:
    # per ogni base temporale, scorrono tutte le basi spaziali.
    # eta_matrix: R x M.
    eta_matrix <- matrix(
      eta_process,
      nrow = n_basis_space,
      ncol = n_basis_time,
      byrow = FALSE
    )
    
    for (time_id in time_ids) {
      
      part_id <- part_id + 1L
      
      psi_t <- basis$Psi[time_id, ]
      
      # Prediction della componente strutturata:
      # Phi_new (n_BAU x R) %*%
      # eta_matrix (R x M) %*%
      # psi_t (M)
      w_hat <- as.numeric(
        Phi_new %*%
          eta_matrix %*%
          psi_t
      )
      
      d <- data.frame(
        process = pollutant_name,
        process_id = process_id,
        space_id = grid$space_id,
        time_id = time_id,
        x = grid$x,
        y = grid$y,
        population_density_scaled =
          grid$population_density_scaled,
        log_altitude_scaled =
          grid$log_altitude_scaled,
        bau_easting = grid$bau_easting,
        bau_northing = grid$bau_northing,
        inside_fit_bbox = grid$inside_fit_bbox,
        w_hat_smooth = w_hat,
        stringsAsFactors = FALSE
      )
      
      time_row <- match(time_id, time_grid$time_id)
      
      if ("t_scaled" %in% names(time_grid)) {
        d$t_scaled <- time_grid$t_scaled[time_row]
      }
      
      if ("sin_year" %in% names(time_grid)) {
        d$sin_year <- time_grid$sin_year[time_row]
      }
      
      if ("cos_year" %in% names(time_grid)) {
        d$cos_year <- time_grid$cos_year[time_row]
      }
      
      if ("sin2_year" %in% names(time_grid)) {
        d$sin2_year <- time_grid$sin2_year[time_row]
      }
      
      if ("cos2_year" %in% names(time_grid)) {
        d$cos2_year <- time_grid$cos2_year[time_row]
      }
      
      if ("year" %in% names(time_grid)) {
        d$year <- factor(
          as.character(time_grid$year[time_row]),
          levels = levels(time_grid$year)
        )
      }
      
      if ("week_start" %in% names(time_grid)) {
        d$week_start <- time_grid$week_start[time_row]
      }
      
      if ("week_end" %in% names(time_grid)) {
        d$week_end <- time_grid$week_end[time_row]
      }
      
      prediction_parts[[part_id]] <- d
      
      if (verbose) {
        cat(
          "\rPrediction:",
          part_id, "/",
          length(process_ids) * length(time_ids),
          "|", pollutant_name,
          "| time_id =", time_id,
          "      "
        )
        flush.console()
      }
    }
  }
  
  data_new <- do.call(rbind, prediction_parts)
  rownames(data_new) <- NULL
  
  if (verbose) cat("\n")
  
  ##########################################################
  # 6. Mean term
  #
  # Nel fit corrente update_beta = FALSE -> beta_hat vuoto,
  # quindi mean_hat = 0.
  ##########################################################
  
  data_new$mean_hat <- 0
  
  if (!is.null(fit$beta_hat) &&
      length(fit$beta_hat) > 0L) {
    
    if (!exists("build_mean_design", mode = "function")) {
      stop(
        "build_mean_design() non disponibile per calcolare ",
        "il mean sui nuovi punti."
      )
    }
    
    X_new <- build_mean_design(
      data = data_new,
      mean_formula = mean_formula,
      n_process = n_process
    )
    
    if (ncol(X_new) != length(fit$beta_hat)) {
      stop(
        "Il design del mean sui nuovi punti non coincide ",
        "con beta_hat del fit."
      )
    }
    
    data_new$mean_hat <- as.numeric(
      X_new %*% fit$beta_hat
    )
  }
  
  # Parte strutturata della prediction.
  data_new$Y_hat_smooth <-
    data_new$mean_hat +
    data_new$w_hat_smooth
  
  # Fine-scale: come nelle simulazioni, viene ricostruita SOLO
  # per combinazioni process-space-time realmente osservate.
  # Sui nuovi BAU / tempi non osservati E[xi | dati] = 0.
  data_new$xi_hat <- 0
  data_new$observed_in_fit <- FALSE
  
  if (include_fine_scale_at_observed) {
    old_space <- real_data$space_grid
    req_space <- c("space_id", "bau_easting", "bau_northing")
    if (!all(req_space %in% names(old_space))) {
      stop("real_data$space_grid non contiene le coordinate metriche necessarie per ricostruire xi.")
    }
    
    # Match ESATTO dei centri BAU: niente nearest-neighbour.
    coord_key_old <- paste(old_space$bau_easting, old_space$bau_northing, sep = "_")
    coord_key_new <- paste(data_new$bau_easting, data_new$bau_northing, sep = "_")
    old_space_id <- old_space$space_id[match(coord_key_new, coord_key_old)]
    
    obs <- real_data$data_obs
    obs_key <- paste(obs$process_id, obs$space_id, obs$time_id, sep = "_")
    pred_key <- paste(data_new$process_id, old_space_id, data_new$time_id, sep = "_")
    obs_idx <- match(pred_key, obs_key)
    matched <- !is.na(obs_idx) & !is.na(old_space_id)
    
    if (any(matched)) {
      pid <- data_new$process_id[matched]
      residual_observed <- obs$Z_obs[obs_idx[matched]] -
        data_new$Y_hat_smooth[matched]
      
      sigma_noise <- if (!is.null(real_data$sigma_noise)) {
        real_data$sigma_noise
      } else if (exists("model", inherits = TRUE) && !is.null(get("model", inherits = TRUE)$sigma_noise)) {
        get("model", inherits = TRUE)$sigma_noise
      } else {
        stop("sigma_noise non disponibile: serve real_data$sigma_noise oppure model$sigma_noise.")
      }
      
      shrinkage <- fit$sigma_fine[pid]^2 /
        (fit$sigma_fine[pid]^2 + sigma_noise[pid]^2)
      
      data_new$xi_hat[matched] <- shrinkage * residual_observed
      data_new$observed_in_fit[matched] <- TRUE
    }
  }
  
  data_new$w_hat <- data_new$w_hat_smooth + data_new$xi_hat
  data_new$Y_hat <- data_new$Y_hat_smooth + data_new$xi_hat
  
  # Conserva prediction raw per diagnostica
  data_new$Z_hat_raw <- data_new$Y_hat
  
  # Prediction fisicamente ammissibile: concentrazioni >= 0
  data_new$Z_hat <- pmax(data_new$Y_hat, 0)
  
  ##########################################################
  # Output
  ##########################################################
  
  out <- list(
    data_full = data_new,
    grid = grid,
    boundary = north_boundary,
    Phi_new = Phi_new,
    time_ids = time_ids,
    pollutants = pollutants,
    include_fine_scale_at_observed = include_fine_scale_at_observed,
    reused_existing_fit = TRUE
  )
  
  class(out) <- c(
    "cofrk_grid_predictions",
    "cofrk_predictions",
    "list"
  )
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Predizioni COFRK - griglia Nord Italia\n")
    cat("============================================\n")
    cat("BAU geografici predetti:", nrow(grid), "\n")
    cat("Righe di prediction:", nrow(data_new), "\n")
    cat("Settimane:", paste(time_ids, collapse = ", "), "\n")
    cat("Pollutant:", paste(pollutants, collapse = ", "), "\n")
    cat(
      "BAU dentro il bounding box del fit:",
      sum(grid$inside_fit_bbox), "/", nrow(grid), "\n"
    )
    cat("Fine-scale ricostruita su:", sum(data_new$observed_in_fit), "righe osservate\n")
    cat("Fit EM riutilizzato: SI\n")
    cat("A_full costruita: NO (prediction tensoriale memory-safe)\n")
  }
  
  out
}

plot_cofrk_northern_italy <- function(
    grid_predictions,
    pollutant = NULL,
    time_id = NULL
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Serve il pacchetto 'ggplot2'.")
  }
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Serve il pacchetto 'sf'.")
  }
  
  d <- grid_predictions$data_full
  
  if (is.null(pollutant)) {
    pollutant <- unique(d$process)[1L]
  }
  
  if (is.null(time_id)) {
    time_id <- unique(d$time_id)[1L]
  }
  
  dd <- d[
    d$process == pollutant &
      d$time_id == time_id,
    ,
    drop = FALSE
  ]
  
  if (nrow(dd) == 0L) {
    stop("Nessuna prediction per pollutant/time_id richiesti.")
  }
  
  subtitle <- paste0(
    "time_id = ",
    time_id
  )
  
  if ("week_start" %in% names(dd)) {
    subtitle <- paste0(
      "Week starting ",
      format(as.Date(dd$week_start[1L]), "%Y-%m-%d")
    )
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = dd,
      ggplot2::aes(
        x = bau_easting,
        y = bau_northing,
        fill = Y_hat
      )
    ) +
    ggplot2::geom_sf(
      data = grid_predictions$boundary,
      fill = NA,
      linewidth = 0.35
    ) +
    ggplot2::coord_sf(
      datum = NA,
      expand = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C"
    ) +
    ggplot2::labs(
      title = paste(
        "COFRK prediction -",
        pollutant
      ),
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = "Prediction"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}


plot_cofrk_northern_italy_facets <- function(
    grid_predictions
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Serve il pacchetto 'ggplot2'.")
  }
  
  d <- grid_predictions$data_full
  
  if ("week_start" %in% names(d)) {
    d$week_label <- format(
      as.Date(d$week_start),
      "%Y-%m-%d"
    )
  } else {
    d$week_label <- paste(
      "Week",
      d$time_id
    )
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = d,
      ggplot2::aes(
        x = bau_easting,
        y = bau_northing,
        fill = Y_hat
      )
    ) +
    ggplot2::geom_sf(
      data = grid_predictions$boundary,
      fill = NA,
      linewidth = 0.25
    ) +
    ggplot2::facet_grid(
      process ~ week_label
    ) +
    ggplot2::coord_sf(
      datum = NA,
      expand = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C"
    ) +
    ggplot2::labs(
      title = "COFRK predictions over Northern Italy",
      subtitle = "Existing EM fit evaluated on the complete 1.6 km BAU grid",
      x = NULL,
      y = NULL,
      fill = "Prediction"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}


############################################################
# USO - NON RIFARE IL FIT
############################################################
#
# source("07_predict_northern_italy_grid.R")
#
# # Default: 4 settimane rappresentative, tutti i pollutant.
# pred_north <- predict_cofrk_northern_italy(
#   real_data = real_data,
#   fit = fit,
#   basis = basis
# )
#
# # Tutte le 4 settimane in un'unica figura:
# plot_cofrk_northern_italy_facets(pred_north)
#
# # Una singola mappa:
# plot_cofrk_northern_italy(
#   pred_north,
#   pollutant = "PM10",
#   time_id = pred_north$time_ids[1]
# )
#
# # Settimane specifiche:
# pred_north <- predict_cofrk_northern_italy(
#   real_data = real_data,
#   fit = fit,
#   basis = basis,
#   time_ids = c(20, 100, 180, 250)
# )
#
# ATTENZIONE:
# Per predire TUTTE le 261 settimane sulla griglia 1.6 km:
#
# pred_north_all <- predict_cofrk_northern_italy(
#   real_data = real_data,
#   fit = fit,
#   basis = basis,
#   time_ids = real_data$time_grid$time_id
# )
#
# Questo oggetto puo' diventare molto grande. Per le figure
# della tesi conviene partire da poche settimane.
############################################################

