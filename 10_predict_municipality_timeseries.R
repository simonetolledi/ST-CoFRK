# Prediction coFRK: un comune x tutto il periodo
predict_cofrk_municipality_timeseries <- function(real_data, fit, basis, comune, comuni_shp,
                                                  pollutants=NULL, cellsize_m=NULL,
                                                  nuts_year="2021", epsg=NULL,
                                                  mean_formula=NULL,
                                                  drop_zero_basis=TRUE, verbose=TRUE) {
  if (!requireNamespace("sf", quietly=TRUE)) stop("Serve 'sf'.")
  
  comuni <- sf::st_read(comuni_shp, quiet=TRUE)
  req <- c("PRO_COM_T","COMUNE","COD_REG")
  if (!all(req %in% names(comuni))) stop("Shapefile ISTAT incompatibile.")
  comuni$PRO_COM_T <- as.character(comuni$PRO_COM_T)
  comuni$COMUNE <- as.character(comuni$COMUNE)
  comuni <- comuni[suppressWarnings(as.integer(comuni$COD_REG)) %in% 1:8,]
  hit <- which(tolower(trimws(comuni$COMUNE)) == tolower(trimws(comune)))
  if (length(hit)!=1L) stop("Comune non trovato o nome non univoco: ", comune)
  comune_sf <- comuni[hit,]
  
  grid_all <- build_northern_italy_prediction_grid(
    real_data=real_data, cellsize_m=cellsize_m, nuts_year=nuts_year, epsg=epsg)
  boundary <- attr(grid_all,"north_boundary")
  crs_grid <- if (!is.null(boundary)) sf::st_crs(boundary) else sf::st_crs(epsg)
  if (is.na(crs_grid)) stop("CRS griglia non determinabile: passa epsg esplicitamente.")
  comune_sf <- sf::st_transform(comune_sf, crs_grid)
  
  if (is.null(cellsize_m)) {
    dx <- diff(sort(unique(grid_all$bau_easting))); dx <- dx[dx>0]
    dy <- diff(sort(unique(grid_all$bau_northing))); dy <- dy[dy>0]
    cellsize_m <- min(c(dx,dy), na.rm=TRUE)
  }
  
  h <- cellsize_m/2
  square <- function(x,y) sf::st_polygon(list(matrix(c(
    x-h,y-h, x+h,y-h, x+h,y+h, x-h,y+h, x-h,y-h), ncol=2, byrow=TRUE)))
  
  cells <- sf::st_sf(
    grid_row=seq_len(nrow(grid_all)),
    geometry=sf::st_sfc(
      lapply(seq_len(nrow(grid_all)), function(i)
        square(grid_all$bau_easting[i],grid_all$bau_northing[i])),
      crs=crs_grid
    )
  )
  
  keep <- lengths(sf::st_intersects(cells, comune_sf))>0
  inter <- suppressWarnings(sf::st_intersection(
    cells[keep,],
    comune_sf[,"PRO_COM_T",drop=FALSE]
  ))
  inter$area_overlap <- as.numeric(sf::st_area(inter))
  aa <- tapply(inter$area_overlap, inter$grid_row, sum, na.rm=TRUE)
  rr <- as.integer(names(aa))
  grid <- grid_all[rr,,drop=FALSE]
  grid$area_overlap <- as.numeric(aa[as.character(rr)])
  
  Phi_new <- evaluate_fitted_spatial_basis_on_grid(grid,basis)
  covered <- rowSums(abs(Phi_new))>0
  
  if (drop_zero_basis) {
    grid <- grid[covered,,drop=FALSE]
    Phi_new <- Phi_new[covered,,drop=FALSE]
  }
  
  if (!nrow(grid)) stop("Comune fuori dal supporto delle basi.")
  grid$space_id <- seq_len(nrow(grid))
  
  tg <- real_data$time_grid
  time_ids <- tg$time_id
  lookup <- real_data$process_lookup
  
  if (is.null(pollutants)) pollutants <- names(lookup)
  if (!all(pollutants %in% names(lookup))) stop("Pollutant non riconosciuto.")
  
  pids <- unname(lookup[pollutants])
  
  R <- ncol(Phi_new)
  M <- ncol(basis$Psi)
  MR <- M*R
  np <- real_data$dimensions$n_process
  
  if (length(fit$eta_hat)!=np*MR)
    stop("eta_hat incompatibile con le basi.")
  
  # Usa esattamente la stessa formula del mean usata nel fit.
  if (is.null(mean_formula) && exists("model", inherits=TRUE)) {
    mm <- get("model", inherits=TRUE)
    if (!is.null(mm$mean_formula)) mean_formula <- mm$mean_formula
  }
  
  if (is.null(mean_formula)) {
    stop(
      "mean_formula non disponibile. Passa mean_formula = model$mean_formula ",
      "oppure assicurati che l'oggetto 'model' sia caricato."
    )
  }
  
  parts <- vector("list",length(pids)*length(time_ids))
  z <- 0L
  
  for (k in seq_along(pids)) {
    pid <- pids[k]
    eta <- fit$eta_hat[((pid-1L)*MR+1L):(pid*MR)]
    E <- matrix(eta, nrow=R, ncol=M)
    for (tid in time_ids) {
      z <- z+1L
      tr <- match(tid,tg$time_id)
      
      
      wh <- as.numeric(
        Phi_new %*% (E %*% basis$Psi[tid, ])
      )
      
      d <- data.frame(
        process=pollutants[k],
        process_id=pid,
        space_id=grid$space_id,
        time_id=tid,
        x=grid$x,
        y=grid$y,
        population_density_scaled=grid$population_density_scaled,
        log_altitude_scaled=grid$log_altitude_scaled,
        bau_easting=grid$bau_easting,
        bau_northing=grid$bau_northing,
        area_overlap=grid$area_overlap,
        w_hat_smooth=wh
      )
      
      for (v in c(
        "t_scaled","sin_year","cos_year",
        "sin2_year","cos2_year",
        "week_start","week_end"
      )) {
        if (v %in% names(tg))
          d[[v]] <- tg[[v]][tr]
      }
      
      if ("year" %in% names(tg))
        d$year <- factor(
          as.character(tg$year[tr]),
          levels=levels(tg$year)
        )
      
      parts[[z]] <- d
    }
  }
  
  dn <- do.call(rbind,parts)
  rownames(dn) <- NULL
  
  dn$mean_hat <- 0
  
  if (!is.null(fit$beta_hat) && length(fit$beta_hat)>0) {
    
    X <- build_mean_design(
      data=dn,
      mean_formula=mean_formula,
      n_process=np
    )
    
    if (ncol(X) != length(fit$beta_hat)) {
      stop(
        "Mean design incompatibile con beta_hat: ",
        ncol(X), " colonne contro ", length(fit$beta_hat), " beta."
      )
    }
    
    if (!is.null(names(fit$beta_hat)) &&
        all(nzchar(names(fit$beta_hat))) &&
        !identical(colnames(X), names(fit$beta_hat))) {
      stop("Nomi/ordine dei beta non coerenti con il design di prediction.")
    }
    
    dn$mean_hat <- as.numeric(X %*% fit$beta_hat)
  }
  
  # Prediction del campo strutturato
  dn$Y_hat <- dn$mean_hat + dn$w_hat_smooth
  dn$Z_hat_raw <- dn$Y_hat
  dn$Z_hat <- pmax(dn$Z_hat_raw,0)
  
  # Aggregazione area-weighted a livello comunale
  sp <- split(
    dn,
    interaction(dn$process,dn$time_id,drop=TRUE)
  )
  
  municipal <- do.call(rbind,lapply(sp,function(dd) {
    
    w <- dd$area_overlap
    
    Z_hat_raw <- weighted.mean(
      dd$Y_hat,
      w,
      na.rm=TRUE
    )
    
    o <- data.frame(
      COMUNE=comune_sf$COMUNE[1],
      PRO_COM_T=comune_sf$PRO_COM_T[1],
      process=dd$process[1],
      process_id=dd$process_id[1],
      time_id=dd$time_id[1],
      Z_hat_raw=Z_hat_raw,
      Z_hat=pmax(Z_hat_raw,0),
      n_bau=nrow(dd)
    )
    
    if ("week_start" %in% names(dd))
      o$week_start <- dd$week_start[1]
    
    if ("week_end" %in% names(dd))
      o$week_end <- dd$week_end[1]
    
    o
  }))
  
  municipal <- municipal[
    order(municipal$process_id,municipal$time_id),
  ]
  rownames(municipal) <- NULL
  
  if (verbose)
    cat(
      "\nComune:",comune_sf$COMUNE[1],
      "| BAU:",nrow(grid),
      "| settimane:",length(time_ids),
      "| fit riutilizzato: SI\n"
    )
  
  list(
    data=municipal,
    bau_predictions=dn,
    bau_grid=grid,
    comune=comune_sf,
    pollutants=pollutants,
    time_ids=time_ids,
    reused_existing_fit=TRUE
  )
}


plot_municipality_timeseries <- function(x, real_data, pollutants=NULL) {
  
  if (!requireNamespace("ggplot2",quietly=TRUE))
    stop("Serve 'ggplot2'.")
  
  if (!requireNamespace("sf",quietly=TRUE))
    stop("Serve 'sf'.")
  
  if (!requireNamespace("dplyr",quietly=TRUE))
    stop("Serve 'dplyr'.")
  
  d <- x$data
  
  if (is.null(pollutants))
    pollutants <- x$pollutants
  
  d <- d[
    d$process %in% pollutants,
    ,
    drop=FALSE
  ]
  
  d$date <- if ("week_start" %in% names(d)) {
    as.Date(d$week_start)
  } else if ("week_end" %in% names(d)) {
    as.Date(d$week_end)
  } else {
    d$time_id
  }
  
  # Osservazioni appartenenti al comune
  obs <- real_data$data_obs
  
  obs$bau_easting <- real_data$space_grid$bau_easting[
    match(obs$space_id,real_data$space_grid$space_id)
  ]
  
  obs$bau_northing <- real_data$space_grid$bau_northing[
    match(obs$space_id,real_data$space_grid$space_id)
  ]
  
  obs_sf <- sf::st_as_sf(
    obs,
    coords=c("bau_easting","bau_northing"),
    crs=sf::st_crs(x$comune),
    remove=FALSE
  )
  
  inside <- lengths(
    sf::st_intersects(obs_sf,x$comune)
  ) > 0
  
  obs_comune <- obs[
    inside & obs$process %in% pollutants,
    ,
    drop=FALSE
  ]
  
  # Media osservata settimanale
  obs_weekly <- obs_comune |>
    dplyr::group_by(process,time_id) |>
    dplyr::summarise(
      Z_obs=mean(Z_obs,na.rm=TRUE),
      n_obs=sum(!is.na(Z_obs)),
      .groups="drop"
    )
  
  comparison <- d |>
    dplyr::left_join(
      obs_weekly,
      by=c("process","time_id")
    )
  
  ggplot2::ggplot(
    comparison,
    ggplot2::aes(x=date)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y=Z_hat,
        colour="coFRK prediction"
      ),
      linewidth=0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        y=Z_obs,
        colour="Observed"
      ),
      size=1.2,
      alpha=0.65,
      na.rm=TRUE
    ) +
    ggplot2::facet_wrap(
      ~process,
      ncol=1,
      scales="free_y"
    ) +
    ggplot2::scale_colour_manual(
      values=c(
        "coFRK prediction"="black",
        "Observed"="#D55E00"
      )
    ) +
    ggplot2::labs(
      title=paste0(
        "Observed and predicted pollutant concentrations - ",
        unique(d$COMUNE)[1]
      ),
      subtitle="Weekly observed BAUs and area-weighted municipal coFRK prediction",
      x=NULL,
      y=expression(
        paste("Concentration (",mu,"g/",m^3,")")
      ),
      colour=NULL
    ) +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(
      panel.grid.minor=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(face="bold"),
      legend.position="top"
    )
}