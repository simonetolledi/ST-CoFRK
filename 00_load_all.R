############################################################
# 00_load_all.R
#
# Pacchetti e caricamento di tutti i moduli.
############################################################

required_packages <- c(
  "MASS",
  "Matrix",
  "ggplot2",
  "FRK",
  "sp",
  "splines",
  "spacetime",
  "zoo"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Installa prima i pacchetti mancanti: ",
    paste(missing_packages, collapse = ", ")
  )
}

project_files <- c(
  "01s_simulation_scenarios.R",
  "02_generate_multires_st_basis.R",
  "03_build_precision_matrices.R",
  "04_build_model_object.R",
  "05_em_fit.R",
  "06_update_precision_parameters.R",
  "07_predict.R",
  "08s_fit_FRK_benchmark.R",
  "09s_metrics_and_plots.R",
  "10s_compare_models.R"
)

missing_files <- project_files[
  !file.exists(project_files)
]

if (length(missing_files) > 0L) {
  stop(
    "File di progetto mancanti: ",
    paste(missing_files, collapse = ", "),
    ". Imposta la working directory sulla cartella del progetto."
  )
}

invisible(
  lapply(project_files, source)
)
