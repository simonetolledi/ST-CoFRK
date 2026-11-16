############################################################
# 09_real_data_plots.R
#
# Richiede:
#   real_data
#   fit
#   predictions <- predict_cofrk(...)
#
# Output:
#   - mappe osservato / predetto per settimane selezionate
#   - serie temporali aggregate osservato / predetto
#   - scatter observed vs fitted
#   - matrice di correlazione stimata tra processi
############################################################

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Serve il pacchetto 'ggplot2'.")
}

.get_prediction_data <- function(predictions) {
  # predict_cofrk() nelle versioni del progetto puo' restituire
  # direttamente un data.frame oppure una lista che lo contiene.
  if (is.data.frame(predictions)) return(predictions)
  
  candidates <- c(
    "predictions", "data", "data_full",
    "prediction_data", "pred_data"
  )
  
  for (nm in candidates) {
    if (!is.null(predictions[[nm]]) &&
        is.data.frame(predictions[[nm]])) {
      return(predictions[[nm]])
    }
  }
  
  stop(
    "Non trovo il data.frame delle predizioni. ",
    "Esegui str(predictions) e controlla l'output."
  )
}


.prepare_plot_data <- function(real_data, predictions) {
  
  pred <- .get_prediction_data(predictions)
  
  required <- c("process_id", "space_id", "time_id")
  missing <- setdiff(required, names(pred))
  
  if (length(missing) > 0L) {
    stop(
      "Nelle predizioni mancano: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # Nome processo
  process_names <- names(real_data$process_lookup)
  
  if (!"process" %in% names(pred)) {
    pred$process <- process_names[pred$process_id]
  }
  
  # Coordinate BAU
  sg <- real_data$space_grid
  
  for (nm in c("x", "y", "bau_easting", "bau_northing")) {
    if (nm %in% names(sg) && !nm %in% names(pred)) {
      pred[[nm]] <- sg[[nm]][
        match(pred$space_id, sg$space_id)
      ]
    }
  }
  
  # Informazioni temporali
  tg <- real_data$time_grid
  
  for (nm in c("week_start", "week_end", "t_scaled")) {
    if (nm %in% names(tg) && !nm %in% names(pred)) {
      pred[[nm]] <- tg[[nm]][
        match(pred$time_id, tg$time_id)
      ]
    }
  }
  
  pred
}


prediction_value_column <- function(pred) {
  # Per le mappe principali preferiamo la prediction del processo Y.
  candidates <- c("Z_hat", "Y_hat", "Y_hat_smooth")
  
  found <- candidates[candidates %in% names(pred)]
  
  if (length(found) == 0L) {
    stop(
      "Non trovo Y_hat, Y_hat_smooth o Z_hat nelle predizioni."
    )
  }
  
  found[1L]
}


select_representative_weeks <- function(real_data, n = 4L) {
  
  ids <- real_data$time_grid$time_id
  
  if (length(ids) <= n) return(ids)
  
  unique(
    as.integer(
      round(seq(1, length(ids), length.out = n))
    )
  )
}


read_northern_italy_municipalities <- function(
    boundary,
    shapefile = file.path(
      "Com01012026_g",
      "Com01012026_g_WGS84.shp"
    )
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Serve il pacchetto 'sf'.")
  }
  
  if (!file.exists(shapefile)) {
    stop(
      "Shapefile dei comuni non trovato: ", shapefile, "\n",
      "La cartella Com01012026_g deve trovarsi nella directory del progetto."
    )
  }
  
  comuni <- sf::st_read(shapefile, quiet = TRUE)
  comuni <- sf::st_make_valid(comuni)
  
  # Gaia/ISTAT: COD_REG da 1 a 8 identifica le regioni del Nord Italia.
  required_istat <- c("PRO_COM_T", "COMUNE", "COD_REG")
  missing_istat <- setdiff(required_istat, names(comuni))
  if (length(missing_istat) > 0L) {
    stop(
      "Nello shapefile ISTAT mancano i campi attesi: ",
      paste(missing_istat, collapse = ", ")
    )
  }
  
  comuni$COD_REG <- as.integer(as.character(comuni$COD_REG))
  comuni <- comuni[comuni$COD_REG %in% 1:8, , drop = FALSE]
  
  # Portiamo i comuni nello stesso CRS usato dalla prediction.
  comuni <- sf::st_transform(comuni, sf::st_crs(boundary))
  
  # Conserviamo esplicitamente identificativo ISTAT, nome e regione.
  comuni[, c("PRO_COM_T", "COMUNE", "COD_REG"), drop = FALSE]
}


aggregate_predictions_to_municipalities <- function(
    real_data,
    predictions,
    municipalities = NULL,
    shapefile = file.path(
      "Com01012026_g",
      "Com01012026_g_WGS84.shp"
    )
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Serve il pacchetto 'sf'.")
  }
  
  pred <- .prepare_plot_data(real_data, predictions)
  value_col <- prediction_value_column(pred)
  
  required <- c("space_id", "bau_easting", "bau_northing", "process", "time_id")
  missing <- setdiff(required, names(pred))
  if (length(missing) > 0L) {
    stop(
      "Per l'aggregazione comunale mancano: ",
      paste(missing, collapse = ", ")
    )
  }
  
  boundary <- predictions$boundary
  if (is.null(boundary)) stop("predictions$boundary non trovato.")
  
  if (is.null(municipalities)) {
    municipalities <- read_northern_italy_municipalities(
      boundary = boundary,
      shapefile = shapefile
    )
  }
  
  # Costruiamo una volta sola i poligoni quadrati dei BAU da 1.6 km.
  # L'aggregazione comunale e' area-weighted: un comune riceve contributo
  # da ogni BAU che lo interseca, proporzionalmente all'area di sovrapposizione.
  bau <- unique(pred[, c("space_id", "bau_easting", "bau_northing"), drop = FALSE])
  
  bau_size <- real_data$settings$spatial_bau_size_m
  if (is.null(bau_size) || !is.finite(bau_size) || bau_size <= 0) {
    stop("real_data$settings$spatial_bau_size_m non valido.")
  }
  half <- bau_size / 2
  
  square_geom <- lapply(seq_len(nrow(bau)), function(i) {
    x <- bau$bau_easting[i]
    y <- bau$bau_northing[i]
    sf::st_polygon(list(matrix(
      c(
        x-half, y-half,
        x+half, y-half,
        x+half, y+half,
        x-half, y+half,
        x-half, y-half
      ),
      ncol = 2,
      byrow = TRUE
    )))
  })
  
  bau_sf <- sf::st_sf(
    space_id = bau$space_id,
    geometry = sf::st_sfc(square_geom, crs = sf::st_crs(boundary))
  )
  
  id_col <- "PRO_COM_T"
  comuni_int <- municipalities[, id_col, drop = FALSE]
  names(comuni_int)[names(comuni_int) == id_col] <- "municipality_id"
  
  message("Calcolo intersezioni area-weighted BAU-comuni (una sola volta)...")
  overlap <- suppressWarnings(sf::st_intersection(bau_sf, comuni_int))
  overlap$overlap_area <- as.numeric(sf::st_area(overlap))
  weights <- sf::st_drop_geometry(overlap)
  weights <- weights[is.finite(weights$overlap_area) & weights$overlap_area > 0, , drop = FALSE]
  weights$municipality_id <- as.character(weights$municipality_id)
  
  # Agganciamo i pesi geometrici alle prediction. Il peso dipende solo dallo spazio,
  # quindi la stessa tabella vale per tutti i pollutant e tutte le settimane.
  pred_small <- pred[, c("space_id", "process", "time_id", value_col), drop = FALSE]
  names(pred_small)[names(pred_small) == value_col] <- "prediction_bau"
  
  pw <- merge(weights, pred_small, by = "space_id", all = FALSE, sort = FALSE)
  pw$weighted_prediction <- pw$overlap_area * pw$prediction_bau
  
  num <- stats::aggregate(
    pw$weighted_prediction,
    by = list(
      municipality_id = pw$municipality_id,
      process = pw$process,
      time_id = pw$time_id
    ),
    FUN = sum,
    na.rm = TRUE
  )
  names(num)[4L] <- "weighted_sum"
  
  den <- stats::aggregate(
    pw$overlap_area,
    by = list(
      municipality_id = pw$municipality_id,
      process = pw$process,
      time_id = pw$time_id
    ),
    FUN = sum,
    na.rm = TRUE
  )
  names(den)[4L] <- "overlap_area"
  
  municipal_mean <- merge(
    num, den,
    by = c("municipality_id", "process", "time_id"),
    all = TRUE,
    sort = FALSE
  )
  municipal_mean$prediction <- municipal_mean$weighted_sum / municipal_mean$overlap_area
  
  municipal_n <- stats::aggregate(
    pw$space_id,
    by = list(
      municipality_id = pw$municipality_id,
      process = pw$process,
      time_id = pw$time_id
    ),
    FUN = function(x) length(unique(x))
  )
  names(municipal_n)[4L] <- "n_bau"
  
  municipal_mean <- merge(
    municipal_mean,
    municipal_n,
    by = c("municipality_id", "process", "time_id"),
    all.x = TRUE,
    sort = FALSE
  )
  
  municipal_mean$week_start <- real_data$time_grid$week_start[
    match(municipal_mean$time_id, real_data$time_grid$time_id)
  ]
  
  list(
    municipalities = municipalities,
    data = municipal_mean,
    value_column = value_col
  )
}


plot_prediction_maps <- function(
    real_data,
    predictions,
    time_ids = NULL,
    pollutants = names(real_data$process_lookup),
    n_weeks = 4L,
    shapefile = file.path(
      "Com01012026_g",
      "Com01012026_g_WGS84.shp"
    )
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Serve il pacchetto 'patchwork'.")
  }
  
  if (is.null(time_ids)) {
    time_ids <- select_representative_weeks(real_data, n = n_weeks)
  }
  
  agg <- aggregate_predictions_to_municipalities(
    real_data = real_data,
    predictions = predictions,
    shapefile = shapefile
  )
  
  comuni <- agg$municipalities

  
  d <- agg$data[
    agg$data$time_id %in% time_ids &
      agg$data$process %in% pollutants,
    ,
    drop = FALSE
  ]
  
  # Livelli delle date presi direttamente dai time_id richiesti: nessun facet NA.
  selected_dates <- real_data$time_grid$week_start[
    match(time_ids, real_data$time_grid$time_id)
  ]
  week_levels <- format(as.Date(selected_dates), "%Y-%m-%d")
  
  id_col <- "PRO_COM_T"
  
  plot_list <- lapply(pollutants, function(p) {
    dp <- d[d$process == p, , drop = FALSE]
    
    map_parts <- lapply(seq_along(time_ids), function(k) {
      tid <- time_ids[k]
      vals <- dp[dp$time_id == tid, , drop = FALSE]
      
      m <- comuni
      m$municipality_id <- as.character(m[[id_col]])
      vals$municipality_id <- as.character(vals$municipality_id)
      
      m <- merge(
        m,
        vals[, c("municipality_id", "prediction", "n_bau")],
        by = "municipality_id",
        all.x = TRUE,
        sort = FALSE
      )
      
      # Assegniamo esplicitamente la data a TUTTI i comuni del pannello,
      # anche a eventuali comuni senza prediction.
      m$week_label <- factor(week_levels[k], levels = week_levels)
      m
    })
    
    map_sf <- do.call(rbind, map_parts)
    lims <- range(dp$prediction, na.rm = TRUE)
    
    ggplot2::ggplot(map_sf) +
      ggplot2::geom_sf(
        ggplot2::aes(fill = prediction),
        colour = NA,
        linewidth = 0
      ) +
      ggplot2::facet_wrap(~ week_label, nrow = 1, drop = FALSE) +
      ggplot2::coord_sf(datum = NA, expand = FALSE) +
      ggplot2::scale_fill_viridis_c(
        option = "C",
        limits = lims,
        na.value = "grey92",
        name = expression(mu*g/m^3)
      ) +
      ggplot2::labs(title = p, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0)
      )
  })
  
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_annotation(
      title = "Predicted pollutant concentrations",
      subtitle = paste0(
        "Weekly COFRK predictions area-weighted from ",
        real_data$settings$spatial_bau_size_m / 1000,
        " km BAUs to municipalities"
      )
    )
}


plot_observed_maps <- function(
    real_data,
    boundary,
    time_ids = NULL,
    pollutants = names(real_data$process_lookup),
    n_weeks = 4L
) {
  
  d <- real_data$data_obs
  
  if (is.null(time_ids)) {
    time_ids <- select_representative_weeks(
      real_data,
      n = n_weeks
    )
  }
  
  d <- d[
    d$time_id %in% time_ids &
      d$process %in% pollutants,
    ,
    drop = FALSE
  ]
  
  sg <- real_data$space_grid
  
  for (nm in c("bau_easting", "bau_northing")) {
    if (nm %in% names(sg) && !nm %in% names(d)) {
      d[[nm]] <- sg[[nm]][
        match(d$space_id, sg$space_id)
      ]
    }
  }
  
  if ("week_start" %in% names(d)) {
    d$week_label <- format(
      as.Date(d$week_start),
      "%Y-%m-%d"
    )
  } else {
    d$week_label <- paste("Week", d$time_id)
  }
  
  # Mantiene l'ordine temporale corretto
  week_levels <- unique(
    d$week_label[order(d$time_id)]
  )
  
  d$week_label <- factor(
    d$week_label,
    levels = week_levels
  )
  
  xvar <- if ("bau_easting" %in% names(d)) {
    "bau_easting"
  } else {
    "x"
  }
  
  yvar <- if ("bau_northing" %in% names(d)) {
    "bau_northing"
  } else {
    "y"
  }
  
  # Un plot separato per pollutant:
  # così ciascuno ha la propria colour scale
  plot_list <- lapply(pollutants, function(p) {
    
    dp <- d[d$process == p, , drop = FALSE]
    
    ggplot2::ggplot(
      dp,
      ggplot2::aes(
        x = .data[[xvar]],
        y = .data[[yvar]]
      )
    ) +
      ggplot2::geom_sf(
        data = boundary,
        inherit.aes = FALSE,
        fill = NA,
        colour = "grey40",
        linewidth = 0.35
      ) +
      ggplot2::geom_point(
        ggplot2::aes(colour = Z_obs),
        size = 2
      ) +
      ggplot2::facet_wrap(
        ~ week_label,
        nrow = 1
      ) +
      ggplot2::coord_sf() +
      ggplot2::scale_colour_viridis_c(
        option = "C",
        na.value = "grey90",
        name = expression(mu*g/m^3)
      ) +
      ggplot2::labs(
        title = p,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0
        )
      )
  })
  
  # Una riga per pollutant
  combined <- patchwork::wrap_plots(
    plot_list,
    ncol = 4
  ) +
    patchwork::plot_annotation(
      title = "Observed pollutant concentrations",
      subtitle = "Weekly BAU averages - pollutant-specific colour scales"
    )
  
  combined
}


plot_observed_vs_fitted <- function(real_data, model, fit) {
  
  obs <- real_data$data_obs
  
  mean_obs <- as.numeric(model$X_obs %*% fit$beta_hat)
  w_obs <- as.numeric(model$A_obs %*% fit$eta_hat)
  
  obs$fitted_raw <- mean_obs + w_obs
  obs$fitted <- pmax(obs$fitted_raw, 0)
  
  ggplot2::ggplot(
    obs,
    ggplot2::aes(
      x = Z_obs,
      y = fitted
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.25,
      size = 1.8
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    ggplot2::facet_wrap(
      ~ process,
      scales = "free"
    ) +
    ggplot2::labs(
      title = "Observed vs fitted",
      subtitle = "In-sample diagnostic",
      x = "Observed weekly concentration",
      y = "Fitted weekly concentration"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}


plot_temporal_response <- function(
    real_data,
    predictions
) {
  
  pred <- .prepare_plot_data(real_data, predictions)
  value_col <- prediction_value_column(pred)
  
  # Media spaziale delle prediction per settimana/processo.
  pred_mean <- stats::aggregate(
    pred[[value_col]],
    by = list(
      process = pred$process,
      time_id = pred$time_id
    ),
    FUN = mean,
    na.rm = TRUE
  )
  
  names(pred_mean)[3L] <- "predicted"
  print(head(pred_mean, 20))
  cat("\nNumero righe pred_mean:", nrow(pred_mean), "\n")
  cat("Processi:", unique(pred_mean$process), "\n")
  cat("Range time_id:", range(pred_mean$time_id), "\n")
  cat("NA predicted:", sum(is.na(pred_mean$predicted)), "\n")
  cat("Range predicted:", range(pred_mean$predicted, na.rm = TRUE), "\n")
  obs <- real_data$data_obs
  
  obs_mean <- stats::aggregate(
    obs$Z_obs,
    by = list(
      process = obs$process,
      time_id = obs$time_id
    ),
    FUN = mean,
    na.rm = TRUE
  )
  
  names(obs_mean)[3L] <- "observed"
  
  d <- merge(
    pred_mean,
    obs_mean,
    by = c("process", "time_id"),
    all = TRUE,
    sort = TRUE
  )
  
  d$week_start <- real_data$time_grid$week_start[
    match(d$time_id, real_data$time_grid$time_id)
  ]
  
  long <- rbind(
    data.frame(
      process = d$process,
      week_start = d$week_start,
      type = "Predicted",
      value = d$predicted
    ),
    data.frame(
      process = d$process,
      week_start = d$week_start,
      type = "Observed",
      value = d$observed
    )
  )
  
  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = week_start,
      y = value,
      colour = type
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ process,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = "Temporal evolution of pollutant concentrations",
      subtitle = "Spatial mean by week",
      x = NULL,
      y = "Concentration",
      colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}


plot_process_correlation <- function(
    fit,
    process_names = c("NO2", "O3", "PM10", "PM25")
) {
  
  if (is.null(fit$Sigma_process)) {
    stop("fit$Sigma_process non trovato.")
  }
  
  R <- stats::cov2cor(
    as.matrix(fit$Sigma_process)
  )
  
  if (nrow(R) != length(process_names)) {
    process_names <- paste0(
      "Process ", seq_len(nrow(R))
    )
  }
  
  rownames(R) <- process_names
  colnames(R) <- process_names
  
  d <- as.data.frame(
    as.table(R),
    stringsAsFactors = FALSE
  )
  
  names(d) <- c(
    "process_1",
    "process_2",
    "correlation"
  )
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = process_1,
      y = process_2,
      fill = correlation
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", correlation)
      ),
      size = 4
    ) +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      midpoint = 0
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Estimated cross-process correlation",
      subtitle = "Derived from Sigma_process",
      x = NULL,
      y = NULL,
      fill = "Correlation"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}


make_real_data_plots <- function(
    real_data,
    model,
    fit,
    predictions,
    time_ids = NULL
) {
  
  list(
    predicted_maps = plot_prediction_maps(
      real_data,
      predictions,
      time_ids = time_ids
    ),
    
    observed_maps = plot_observed_maps(
      real_data,
      boundary = predictions$boundary,
      time_ids = time_ids
    ),
    
    observed_vs_fitted = plot_observed_vs_fitted(
      real_data = real_data,
      model = model,
      fit = fit
    ),
    
    temporal_response = plot_temporal_response(
      real_data,
      predictions
    ),
    
    process_correlation = plot_process_correlation(
      fit,
      process_names = names(real_data$process_lookup)
    )
  )
}


############################################################
# ESEMPIO DOPO IL FIT:
#
# predictions <- predict_cofrk(
#   model = model,
#   fit = fit,
#   basis = basis,
#   verbose = TRUE
# )
#
# plots <- make_real_data_plots(
#   real_data = real_data,
#   fit = fit,
#   predictions = predictions
# )
#
# plots$predicted_maps
# plots$observed_maps
# plots$observed_vs_fitted
# plots$temporal_response
# plots$process_correlation
#
# Settimane specifiche:
# plots <- make_real_data_plots(
#   real_data,
#   fit,
#   predictions,
#   time_ids = c(20, 100, 180, 250)
# )
############################################################

