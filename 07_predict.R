############################################################
# 07_predict.R
#
# Predizione del processo:
#   Y = mean + w + xi
# e dell'osservazione:
#   Z = Y + epsilon
############################################################

match_observed_rows <- function(data_full, data_obs) {
  full_key <- paste(
    data_full$process_id,
    data_full$space_id,
    data_full$time_id,
    sep = "_"
  )
  
  obs_key <- paste(
    data_obs$process_id,
    data_obs$space_id,
    data_obs$time_id,
    sep = "_"
  )
  
  match(obs_key, full_key)
}

predict_cofrk <- function(
    model,
    fit,
    basis,
    include_fine_scale_at_observed = TRUE,
    verbose = TRUE
) {
  data_full <- model$data_full
  data_obs <- model$data_obs
  n_process <- model$dimensions$n_process
  
  A_full <- build_multivariate_design(
    data = data_full,
    Phi = basis$Phi,
    Psi = basis$Psi,
    n_process = n_process
  )
  
  data_full$w_hat_smooth <- as.numeric(
    A_full %*% fit$eta_hat
  )
  
  if (!is.null(fit$beta_hat) && length(fit$beta_hat) > 0L) {
    data_full$mean_hat <- as.numeric(
      model$X_full %*% fit$beta_hat
    )
  }
  
  data_full$Y_hat_smooth <- data_full$mean_hat +
    data_full$w_hat_smooth
  
  data_full$xi_hat <- 0
  
  observed_rows <- match_observed_rows(
    data_full = data_full,
    data_obs = data_obs
  )
  
  if (anyNA(observed_rows)) {
    stop("Non è stato possibile associare tutti i dati osservati.")
  }
  
  if (include_fine_scale_at_observed) {
    residual_observed <- data_obs$Z_obs -
      data_full$Y_hat_smooth[observed_rows]
    
    shrinkage <- fit$sigma_fine[
      data_obs$process_id
    ]^2 / (
      fit$sigma_fine[data_obs$process_id]^2 +
        model$sigma_noise[data_obs$process_id]^2
    )
    
    data_full$xi_hat[observed_rows] <-
      shrinkage * residual_observed
  }
  
  # Processo latente stimato.
  data_full$w_hat <- data_full$w_hat_smooth +
    data_full$xi_hat
  data_full$Y_hat <- data_full$Y_hat_smooth +
    data_full$xi_hat
  
  # Media predittiva dell'osservazione:
  # E[Z | dati] = E[Y | dati], perché E[epsilon] = 0.
  data_full$Z_hat <- data_full$Y_hat
  
  predictions <- list(
    data_full = data_full,
    A_full = A_full,
    observed_rows = observed_rows,
    include_fine_scale_at_observed =
      include_fine_scale_at_observed
  )
  
  class(predictions) <- c("cofrk_predictions", "list")
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Predizioni COFRK\n")
    cat("============================================\n")
    cat("Predizioni totali:", nrow(data_full), "\n")
    cat(
      "Fine-scale ricostruita su:",
      sum(data_full$xi_hat != 0),
      "punti osservati\n"
    )
  }
  
  predictions
}
