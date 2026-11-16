############################################################
# 09_metrics_and_plots.R
#
# Metriche, grafici e salvataggio dei risultati.
############################################################

scale_fill_field <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c(
      "purple4", "royalblue3", "white", "gold", "red3"
    )
  )
}

scale_fill_error <- function() {
  ggplot2::scale_fill_gradient2(
    low = "red3",
    mid = "white",
    high = "blue3",
    midpoint = 0
  )
}

scale_color_models <- function() {
  ggplot2::scale_color_manual(
    values = c(
      "FRK" = "royalblue3",
      "ST-CoFRK" = "red3"
    ),
    limits = c("ST-CoFRK", "FRK"),
    drop = FALSE
  )
}

compute_model_metrics <- function(data_full) {
  data_full$error_COFRK <- data_full$Y_hat -
    data_full$Y_true
  data_full$error_COFRK_smooth <-
    data_full$Y_hat_smooth - data_full$Y_true
  data_full$error_FRK <- data_full$Y_hat_FRK -
    data_full$Y_true
  
  metric_formula <- stats::as.formula(
    "cbind(
       MSE_FRK_univariato = error_FRK^2,
       MSE_COFRK_multivariato = error_COFRK^2,
       MSE_COFRK_multivariato_smooth =
         error_COFRK_smooth^2
     ) ~ process"
  )
  
  metrics_global <- stats::aggregate(
    metric_formula,
    data = data_full,
    FUN = mean
  )
  
  metrics_missing <- stats::aggregate(
    metric_formula,
    data = data_full[!data_full$observed, , drop = FALSE],
    FUN = mean
  )
  
  metrics_observed <- stats::aggregate(
    metric_formula,
    data = data_full[data_full$observed, , drop = FALSE],
    FUN = mean
  )
  
  metrics_time_y2 <- stats::aggregate(
    cbind(
      FRK_ST_univariato = error_FRK^2,
      COFRK_multivariato = error_COFRK^2
    ) ~ time_id,
    data = data_full[data_full$process == "Y2", , drop = FALSE],
    FUN = mean
  )
  
  list(
    data_full = data_full,
    global = metrics_global,
    missing = metrics_missing,
    observed = metrics_observed,
    time_y2 = metrics_time_y2
  )
}

make_model_plots <- function(
    metrics,
    scenario_label,
    process_to_plot = "Y2",
    times_to_plot = c(5, 15, 25, 35),
    selected_points = data.frame(
      x_target = c(
        0.25, 0.30, 0.45,
        0.60, 0.75, 0.90,
        0.00, 1.00, 0.50
      ),
      y_target = c(
        0.45, 0.90, 0.25,
        0.30, 0.75, 0.60,
        0.00, 1.00, 0.50
      )
    )
) {
  data_full <- metrics$data_full
  
  ##########################################################
  # 1. Mappe vero vs predetto
  ##########################################################
  
  plot_data <- data_full[
    data_full$process == process_to_plot &
      data_full$time_id %in% times_to_plot,
    ,
    drop = FALSE
  ]
  
  comparison_long <- data.frame(
    x = rep(plot_data$x, 3L),
    y = rep(plot_data$y, 3L),
    time_id = rep(plot_data$time_id, 3L),
    value = c(
      plot_data$Y_true,
      plot_data$Y_hat_FRK,
      plot_data$Y_hat
    ),
    type = rep(
      c(
        "True",
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(plot_data)
    )
  )
  
  comparison_map <- ggplot2::ggplot(
    comparison_long,
    ggplot2::aes(
      x = x,
      y = y,
      fill = value
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      type ~ time_id,
      scales = "fixed"
    ) +
    scale_fill_field() +
    ggplot2::labs(
      
      x = "x",
      y = "y",
      fill = "Value"
    )
  
  ##########################################################
  # 2. Mappe degli errori
  ##########################################################
  
  error_long <- data.frame(
    x = rep(plot_data$x, 2L),
    y = rep(plot_data$y, 2L),
    time_id = rep(plot_data$time_id, 2L),
    error = c(
      plot_data$error_FRK,
      plot_data$error_COFRK
    ),
    model = rep(
      c(
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(plot_data)
    )
  )
  
  error_map <- ggplot2::ggplot(
    error_long,
    ggplot2::aes(
      x = x,
      y = y,
      fill = error
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      model ~ time_id,
      scales = "fixed"
    ) +
    scale_fill_error() +
    ggplot2::labs(
      
      x = "x",
      y = "y",
      fill = "Error"
    )
  
  ##########################################################
  # 3. MSE nel tempo
  ##########################################################
  
  mse_time_long <- reshape(
    metrics$time_y2,
    varying = c(
      "FRK_ST_univariato",
      "COFRK_multivariato"
    ),
    v.names = "MSE",
    timevar = "model",
    times = c(
      "FRK",
      "ST-CoFRK"
    ),
    direction = "long"
  )
  
  mse_time_long$model <- factor(
    mse_time_long$model,
    levels = c("ST-CoFRK", "FRK")
  )
  
  mse_time_plot <- ggplot2::ggplot(
    mse_time_long,
    ggplot2::aes(
      x = time_id,
      y = MSE,
      color = model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    scale_color_models() +
    ggplot2::labs(
      
      x = "Time",
      y = "MSE",
      color = "Model"
    )
  
  ##########################################################
  # 4. Serie temporali in locazioni specifiche
  ##########################################################
  
  space_grid <- unique(
    data_full[
      ,
      c("space_id", "x", "y"),
      drop = FALSE
    ]
  )
  
  space_grid <- space_grid[
    order(space_grid$space_id),
    ,
    drop = FALSE
  ]
  
  selected_space_ids <- integer(
    nrow(selected_points)
  )
  
  for (location_index in seq_len(nrow(selected_points))) {
    
    squared_distance <-
      (
        space_grid$x -
          selected_points$x_target[location_index]
      )^2 +
      (
        space_grid$y -
          selected_points$y_target[location_index]
      )^2
    
    selected_space_ids[location_index] <-
      space_grid$space_id[
        which.min(squared_distance)
      ]
  }
  
  selected_locations <- space_grid[
    match(
      selected_space_ids,
      space_grid$space_id
    ),
    ,
    drop = FALSE
  ]
  
  selected_locations$location <- paste0(
    "Location ",
    seq_len(nrow(selected_locations)),
    " (x=",
    round(selected_locations$x, 2),
    ", y=",
    round(selected_locations$y, 2),
    ")"
  )
  
  time_location_data <- merge(
    data_full[
      data_full$process == process_to_plot &
        data_full$space_id %in% selected_space_ids,
      ,
      drop = FALSE
    ],
    selected_locations[
      ,
      c("space_id", "location"),
      drop = FALSE
    ],
    by = "space_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  time_location_long <- data.frame(
    time_id = rep(
      time_location_data$time_id,
      3L
    ),
    value = c(
      time_location_data$Y_true,
      time_location_data$Y_hat_FRK,
      time_location_data$Y_hat
    ),
    model = rep(
      c(
        "True",
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(time_location_data)
    ),
    location = rep(
      time_location_data$location,
      3L
    )
  )
  
  time_locations_plot <- ggplot2::ggplot(
    time_location_long,
    ggplot2::aes(
      x = time_id,
      y = value,
      color = model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(
      ~ location,
      scales = "fixed"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "True" = "black",
        "FRK" = "royalblue3",
        "ST-CoFRK" = "red3"
      )
    ) +
    ggplot2::labs(
      
      x = "Time",
      y = "Value",
      color = "Series"
    )
  
  ##########################################################
  # Output
  ##########################################################
  
  list(
    comparison_map = comparison_map,
    error_map = error_map,
    mse_time = mse_time_plot,
    time_locations = time_locations_plot,
    selected_locations = selected_locations
  )
}

plot_observed_points_times <- function(
    simulation,
    process_id = NULL,
    times = c(1L, 10L, 25L, 40L, 50L)
) {
  data_full <- simulation$data_full
  available_ids <- sort(unique(data_full$process_id))
  
  if (is.null(process_id)) {
    # If Y2 is still labeled in the simulation, select it directly.
    if ("process" %in% names(data_full) && any(data_full$process == "Y2")) {
      y2_ids <- unique(data_full$process_id[data_full$process == "Y2"])
      if (length(y2_ids) != 1L) stop("Y2 has multiple process IDs.")
      process_id <- y2_ids
    } else if (length(available_ids) == 1L) {
      process_id <- available_ids
    } else if (2L %in% available_ids) {
      process_id <- 2L
    } else {
      stop("Cannot identify Y2 automatically; specify process_id.")
    }
  }
  
  if (length(process_id) != 1L || is.na(process_id) ||
      !process_id %in% available_ids) {
    stop("Invalid process_id. Available IDs: ", paste(available_ids, collapse = ", "))
  }
  
  plot_data <- data_full[
    data_full$process_id == process_id & data_full$time_id %in% times,
    , drop = FALSE
  ]
  if (nrow(plot_data) == 0L) {
    stop("No observations for process_id=", process_id,
         " and selected times. Available times: ",
         paste(range(data_full$time_id), collapse = "-"))
  }
  
  plot_data$time_label <- factor(
    paste0("t = ", plot_data$time_id),
    levels = paste0("t = ", times)
  )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(colour = "grey85", size = 1.5) +
    ggplot2::geom_point(
      data = plot_data[!is.na(plot_data$observed) & plot_data$observed, , drop = FALSE],
      colour = "black", size = 2
    ) +
    ggplot2::facet_wrap(~ time_label, nrow = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "x", y = "y") +
    ggplot2::theme_minimal()
}


save_model_results <- function(
    metrics,
    plots,
    output_directory
) {
  plots_directory <- file.path(
    output_directory,
    "plots"
  )
  dir.create(
    plots_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  utils::write.csv(
    metrics$global,
    file.path(output_directory, "metrics_global.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    metrics$missing,
    file.path(output_directory, "metrics_missing.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    metrics$observed,
    file.path(output_directory, "metrics_observed.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    metrics$time_y2,
    file.path(output_directory, "mse_time_y2.csv"),
    row.names = FALSE
  )
  
  ggplot2::ggsave(
    file.path(plots_directory, "Y2_true_vs_predicted.png"),
    plots$comparison_map,
    width = 11,
    height = 7,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(plots_directory, "Y2_prediction_error.png"),
    plots$error_map,
    width = 11,
    height = 5,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(plots_directory, "Y2_MSE_over_time.png"),
    plots$mse_time,
    width = 9,
    height = 5,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_time_series_selected_locations.png"
    ),
    plots$time_locations,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  invisible(output_directory)
}



############################################################
# Aggregazione Monte Carlo
############################################################

row_variance <- function(x) {
  if (ncol(x) <= 1L) {
    return(rep(0, nrow(x)))
  }
  apply(x, 1L, stats::var)
}


aggregate_monte_carlo_results <- function(run_results) {
  if (length(run_results) < 1L) {
    stop("run_results deve contenere almeno una run.")
  }
  
  data_list <- lapply(
    run_results,
    function(run) run$metrics$data_full
  )
  
  reference_data <- data_list[[1L]]
  
  reference_key <- paste(
    reference_data$process_id,
    reference_data$space_id,
    reference_data$time_id,
    sep = "_"
  )
  
  for (run_id in seq_along(data_list)) {
    current_key <- paste(
      data_list[[run_id]]$process_id,
      data_list[[run_id]]$space_id,
      data_list[[run_id]]$time_id,
      sep = "_"
    )
    
    if (!identical(reference_key, current_key)) {
      stop(
        "L'ordinamento dei dati completi non coincide tra le run."
      )
    }
  }
  
  y_true_matrix <- do.call(
    cbind,
    lapply(data_list, function(x) x$Y_true)
  )
  
  z_obs_matrix <- do.call(
    cbind,
    lapply(data_list, function(x) x$Z_obs)
  )
  
  frk_matrix <- do.call(
    cbind,
    lapply(data_list, function(x) x$Y_hat_FRK)
  )
  
  cofrk_matrix <- do.call(
    cbind,
    lapply(data_list, function(x) x$Y_hat)
  )
  
  cofrk_smooth_matrix <- do.call(
    cbind,
    lapply(data_list, function(x) x$Y_hat_smooth)
  )
  
  observed_matrix <- do.call(
    cbind,
    lapply(
      data_list,
      function(x) as.numeric(x$observed)
    )
  )
  
  aggregated_data <- reference_data[
    ,
    c(
      "process",
      "process_id",
      "x",
      "y",
      "space_id",
      "time_id",
      "t_scaled",
      "st_bau_id"
    ),
    drop = FALSE
  ]
  
  aggregated_data$Y_true <- rowMeans(y_true_matrix)
  aggregated_data$Y_true_variance <- row_variance(y_true_matrix)
  
  aggregated_data$Z_obs <- rowMeans(z_obs_matrix)
  aggregated_data$Z_obs_variance <- row_variance(z_obs_matrix)
  
  aggregated_data$Y_hat_FRK_mean <- rowMeans(frk_matrix)
  aggregated_data$Y_hat_FRK_variance <- row_variance(frk_matrix)
  aggregated_data$Y_hat_FRK_sd <-
    sqrt(aggregated_data$Y_hat_FRK_variance)
  
  aggregated_data$Y_hat_COFRK_mean <- rowMeans(cofrk_matrix)
  aggregated_data$Y_hat_COFRK_variance <- row_variance(cofrk_matrix)
  aggregated_data$Y_hat_COFRK_sd <-
    sqrt(aggregated_data$Y_hat_COFRK_variance)
  
  aggregated_data$Y_hat_COFRK_smooth_mean <-
    rowMeans(cofrk_smooth_matrix)
  aggregated_data$Y_hat_COFRK_smooth_variance <-
    row_variance(cofrk_smooth_matrix)
  
  aggregated_data$error_FRK_mean <-
    aggregated_data$Y_hat_FRK_mean -
    aggregated_data$Y_true
  
  aggregated_data$error_COFRK_mean <-
    aggregated_data$Y_hat_COFRK_mean -
    aggregated_data$Y_true
  
  # Y_true è identico tra le run, quindi la varianza
  # dell'errore coincide con la varianza della predizione.
  aggregated_data$error_FRK_variance <-
    aggregated_data$Y_hat_FRK_variance
  
  aggregated_data$error_COFRK_variance <-
    aggregated_data$Y_hat_COFRK_variance
  
  aggregated_data$observed_frequency <-
    rowMeans(observed_matrix)
  
  time_runs <- do.call(
    rbind,
    lapply(
      seq_along(run_results),
      function(run_id) {
        out <- run_results[[run_id]]$metrics$time_y2
        out$run_id <- run_id
        out
      }
    )
  )
  
  summarize_time_model <- function(
    values,
    time_id,
    model_name
  ) {
    time_values <- sort(unique(time_id))
    
    do.call(
      rbind,
      lapply(
        time_values,
        function(current_time) {
          x <- values[time_id == current_time]
          x <- x[is.finite(x)]
          
          data.frame(
            time_id = current_time,
            model = model_name,
            mean = mean(x),
            variance = if (length(x) > 1L) {
              stats::var(x)
            } else {
              0
            },
            sd = if (length(x) > 1L) {
              stats::sd(x)
            } else {
              0
            }
          )
        }
      )
    )
  }
  
  time_summary <- rbind(
    summarize_time_model(
      values = time_runs$FRK_ST_univariato,
      time_id = time_runs$time_id,
      model_name = "FRK"
    ),
    summarize_time_model(
      values = time_runs$COFRK_multivariato,
      time_id = time_runs$time_id,
      model_name = "ST-CoFRK"
    )
  )
  
  rownames(time_summary) <- NULL
  
  list(
    n_runs = length(run_results),
    data_full = aggregated_data,
    time_y2 = time_summary,
    time_y2_runs = time_runs
  )
}


############################################################
# Grafici Monte Carlo
############################################################

make_monte_carlo_plots <- function(
    aggregated,
    scenario_label,
    process_to_plot = "Y2",
    times_to_plot = c(5, 15, 25, 35, 35),
    selected_points = data.frame(
      x_target = c(
        0.10, 0.30, 0.45,
        0.60, 0.75, 0.90,
        0.00, 1.00, 0.50
      ),
      y_target = c(
        0.50, 0.90, 0.25,
        0.30, 0.75, 0.60,
        0.00, 1.00, 0.50
      )
    )
) {
  data_full <- aggregated$data_full
  
  plot_data <- data_full[
    data_full$process == process_to_plot &
      data_full$time_id %in% times_to_plot,
    ,
    drop = FALSE
  ]
  
  ##########################################################
  # 1. Vero vs predizione media
  ##########################################################
  
  comparison_long <- data.frame(
    x = rep(plot_data$x, 3L),
    y = rep(plot_data$y, 3L),
    time_id = rep(plot_data$time_id, 3L),
    value = c(
      plot_data$Y_true,
      plot_data$Y_hat_FRK_mean,
      plot_data$Y_hat_COFRK_mean
    ),
    type = rep(
      c(
        "True",
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(plot_data)
    )
  )
  
  comparison_map <- ggplot2::ggplot(
    comparison_long,
    ggplot2::aes(
      x = x,
      y = y,
      fill = value
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      type ~ time_id,
      scales = "fixed"
    ) +
    scale_fill_field() +
    ggplot2::labs(
      
      x = "x",
      y = "y",
      fill = "Value"
    )
  
  ##########################################################
  # 2. Errore medio
  ##########################################################
  
  error_long <- data.frame(
    x = rep(plot_data$x, 2L),
    y = rep(plot_data$y, 2L),
    time_id = rep(plot_data$time_id, 2L),
    error = c(
      plot_data$error_FRK_mean,
      plot_data$error_COFRK_mean
    ),
    model = rep(
      c(
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(plot_data)
    )
  )
  
  error_map <- ggplot2::ggplot(
    error_long,
    ggplot2::aes(
      x = x,
      y = y,
      fill = error
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      model ~ time_id,
      scales = "fixed"
    ) +
    scale_fill_error() +
    ggplot2::labs(
      
      x = "x",
      y = "y",
      fill = "Mean error"
    )
  
  ##########################################################
  # 3. Varianza della predizione tra run
  ##########################################################
  
  variance_long <- data.frame(
    x = rep(plot_data$x, 2L),
    y = rep(plot_data$y, 2L),
    time_id = rep(plot_data$time_id, 2L),
    variance = c(
      plot_data$Y_hat_FRK_variance,
      plot_data$Y_hat_COFRK_variance
    ),
    model = rep(
      c(
        "FRK",
        "ST-CoFRK"
      ),
      each = nrow(plot_data)
    )
  )
  
  prediction_variance_map <- ggplot2::ggplot(
    variance_long,
    ggplot2::aes(
      x = x,
      y = y,
      fill = variance
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      model ~ time_id,
      scales = "fixed"
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "red3"
    ) +
    ggplot2::labs(
      
      x = "x",
      y = "y",
      fill = "Variance"
    )
  
  ##########################################################
  # 4. MSE medio nel tempo con banda +/- 1 sd
  ##########################################################
  
  time_plot_data <- aggregated$time_y2
  # Normalize labels also when aggregated results were created with an older script.
  # Never change the underlying MSE or variance values.
  model_labels <- trimws(as.character(time_plot_data$model))
  model_labels[model_labels %in% c(
    "FRK ST univariato", "FRK_ST_univariato", "FRK"
  )] <- "FRK"
  model_labels[model_labels %in% c(
    "COFRK ST multivariato", "COFRK_multivariato", "COFRK ST multivariato - media",
    "ST-CoFRK", "COFRK"
  )] <- "ST-CoFRK"
  unexpected_models <- setdiff(unique(model_labels), c("ST-CoFRK", "FRK"))
  if (length(unexpected_models)) {
    stop("Unknown model labels in aggregated$time_y2: ",
         paste(unexpected_models, collapse = ", "))
  }
  time_plot_data$model <- factor(model_labels, levels = c("ST-CoFRK", "FRK"))
  
  mse_time_plot <- ggplot2::ggplot(
    time_plot_data,
    ggplot2::aes(
      x = time_id,
      y = mean,
      color = model,
      fill = model
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = pmax(mean - sd, 0),
        ymax = mean + sd
      ),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.2) +
    scale_color_models() +
    ggplot2::scale_fill_manual(
      values = c(
        "FRK" = "royalblue3",
        "ST-CoFRK" = "red3"
      ),
      limits = c("ST-CoFRK", "FRK"),
      drop = FALSE
    ) +
    ggplot2::labs(
      
      x = "Time",
      y = "Mean MSE",
      color = "Model",
      fill = "Model"
    )
  
  ##########################################################
  # 5. Varianza del MSE nel tempo
  ##########################################################
  
  mse_variance_plot <- ggplot2::ggplot(
    time_plot_data,
    ggplot2::aes(
      x = time_id,
      y = variance,
      color = model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 1.8) +
    scale_color_models() +
    ggplot2::labs(
      
      x = "Time",
      y = "MSE variance",
      color = "Model"
    )
  
  ##########################################################
  # 6. Serie temporali medie in locazioni selezionate
  ##########################################################
  
  space_grid <- unique(
    data_full[
      ,
      c("space_id", "x", "y"),
      drop = FALSE
    ]
  )
  
  space_grid <- space_grid[
    order(space_grid$space_id),
    ,
    drop = FALSE
  ]
  
  selected_space_ids <- integer(
    nrow(selected_points)
  )
  
  for (location_index in seq_len(nrow(selected_points))) {
    squared_distance <-
      (
        space_grid$x -
          selected_points$x_target[location_index]
      )^2 +
      (
        space_grid$y -
          selected_points$y_target[location_index]
      )^2
    
    selected_space_ids[location_index] <-
      space_grid$space_id[
        which.min(squared_distance)
      ]
  }
  
  selected_locations <- space_grid[
    match(
      selected_space_ids,
      space_grid$space_id
    ),
    ,
    drop = FALSE
  ]
  
  selected_locations$location <- paste0(
    "Location ",
    seq_len(nrow(selected_locations)),
    " (x=",
    round(selected_locations$x, 2),
    ", y=",
    round(selected_locations$y, 2),
    ")"
  )
  
  time_location_data <- merge(
    data_full[
      data_full$process == process_to_plot &
        data_full$space_id %in% selected_space_ids,
      ,
      drop = FALSE
    ],
    selected_locations[
      ,
      c("space_id", "location"),
      drop = FALSE
    ],
    by = "space_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  true_data <- data.frame(
    time_id = time_location_data$time_id,
    location = time_location_data$location,
    mean = time_location_data$Y_true,
    sd = 0,
    model = "True"
  )
  
  frk_data <- data.frame(
    time_id = time_location_data$time_id,
    location = time_location_data$location,
    mean = time_location_data$Y_hat_FRK_mean,
    sd = time_location_data$Y_hat_FRK_sd,
    model = "FRK"
  )
  
  cofrk_data <- data.frame(
    time_id = time_location_data$time_id,
    location = time_location_data$location,
    mean = time_location_data$Y_hat_COFRK_mean,
    sd = time_location_data$Y_hat_COFRK_sd,
    model = "ST-CoFRK"
  )
  
  time_location_long <- rbind(
    true_data,
    frk_data,
    cofrk_data
  )
  
  time_locations_plot <- ggplot2::ggplot(
    time_location_long,
    ggplot2::aes(
      x = time_id,
      y = mean,
      color = model,
      fill = model
    )
  ) +
    ggplot2::geom_ribbon(
      data = time_location_long[
        time_location_long$model != "True",
        ,
        drop = FALSE
      ],
      ggplot2::aes(
        ymin = mean - sd,
        ymax = mean + sd
      ),
      alpha = 0.12,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::facet_wrap(
      ~ location,
      scales = "fixed"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "True" = "black",
        "FRK" = "royalblue3",
        "ST-CoFRK" = "red3"
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "True" = "black",
        "FRK" = "royalblue3",
        "ST-CoFRK" = "red3"
      )
    ) +
    ggplot2::labs(
      
      x = "Time",
      y = "Value",
      color = "Series",
      fill = "Series"
    )
  
  list(
    comparison_map = comparison_map,
    error_map = error_map,
    prediction_variance_map = prediction_variance_map,
    mse_time = mse_time_plot,
    mse_variance = mse_variance_plot,
    time_locations = time_locations_plot,
    selected_locations = selected_locations
  )
}


############################################################
# Salvataggio risultati Monte Carlo
############################################################

save_monte_carlo_results <- function(
    aggregated,
    plots,
    summary_runs,
    monte_carlo_summary,
    output_directory
) {
  plots_directory <- file.path(
    output_directory,
    "plots"
  )
  
  dir.create(
    plots_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  utils::write.csv(
    summary_runs,
    file.path(
      output_directory,
      "all_runs_summary.csv"
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    monte_carlo_summary,
    file.path(
      output_directory,
      "monte_carlo_summary.csv"
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    aggregated$time_y2,
    file.path(
      output_directory,
      "mse_time_y2_mean_variance.csv"
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    aggregated$data_full,
    file.path(
      output_directory,
      "prediction_mean_variance.csv"
    ),
    row.names = FALSE
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_true_vs_mean_prediction.png"
    ),
    plots$comparison_map,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_mean_prediction_error.png"
    ),
    plots$error_map,
    width = 11,
    height = 5,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_prediction_variance.png"
    ),
    plots$prediction_variance_map,
    width = 11,
    height = 5,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_MSE_mean_over_time.png"
    ),
    plots$mse_time,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_MSE_variance_over_time.png"
    ),
    plots$mse_variance,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(
      plots_directory,
      "Y2_time_series_mean_sd.png"
    ),
    plots$time_locations,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  invisible(output_directory)
}