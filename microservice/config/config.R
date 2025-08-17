# Basic service and model configuration

safe_get_env <- function(key, default = "") {
  val <- Sys.getenv(key, unset = NA)
  if (is.na(val) || identical(val, "")) {
    return(default)
  }
  val
}

# Try to resolve repository root in both local and container
resolve_app_root <- function() {
  # Priority: APP_ROOT env
  app_root <- safe_get_env("APP_ROOT", "")
  if (!identical(app_root, "") && dir.exists(app_root)) {
    return(normalizePath(app_root, mustWork = FALSE))
  }

  # If running with WORKDIR=/app in Docker
  candidate <- "."
  if (file.exists(file.path(candidate, "RF-AmPEP30.R"))) {
    return(normalizePath(candidate, mustWork = FALSE))
  }

  # If running from microservice/api directory locally
  candidate <- file.path("..", "..")
  if (file.exists(file.path(candidate, "RF-AmPEP30.R"))) {
    return(normalizePath(candidate, mustWork = FALSE))
  }

  # Fallback: current working directory
  normalizePath(getwd(), mustWork = FALSE)
}

APP_ROOT <- resolve_app_root()

SERVICE_CONFIG <- list(
  name = safe_get_env("SERVICE_NAME", "deep-ampep30-microservice"),
  version = safe_get_env("SERVICE_VERSION", "1.0.0"),
  host = safe_get_env("PLUMBER_HOST", "0.0.0.0"),
  port = as.integer(safe_get_env("PLUMBER_PORT", "8002")),
  log_level = safe_get_env("R_LOG_LEVEL", "INFO"),
  default_method = tolower(safe_get_env("DEFAULT_METHOD", "rf"))
)

MODEL_CONFIG <- list(
  rf_model_path = safe_get_env("MODEL_RF_PATH", file.path(APP_ROOT, "AmPEP30-RF-1200tree.mdl")),
  cnn_model_path = safe_get_env("MODEL_CNN_PATH", file.path(APP_ROOT, "AmPEP30-CNN.mdl")),
  enable_rf = tolower(safe_get_env("ENABLE_RF", "true")) %in% c("1", "true", "yes"),
  enable_cnn = tolower(safe_get_env("ENABLE_CNN", "false")) %in% c("1", "true", "yes")
)

API_CONFIG <- list(
  max_request_size_bytes = as.integer(safe_get_env("MAX_REQUEST_SIZE_BYTES", as.character(10 * 1024 * 1024))),
  timeout_seconds = as.integer(safe_get_env("API_TIMEOUT_SECONDS", "3600")),
  cors_enabled = tolower(safe_get_env("CORS_ENABLED", "true")) %in% c("1", "true", "yes")
)

VALIDATION_CONFIG <- list(
  min_sequence_length = as.integer(safe_get_env("MIN_SEQUENCE_LENGTH", "5")),
  max_sequence_length = as.integer(safe_get_env("MAX_SEQUENCE_LENGTH", "30")),
  allowed_amino_acids = safe_get_env("ALLOWED_AMINO_ACIDS", "ACDEFGHIKLMNPQRSTVWY"),
  max_sequences_per_request = as.integer(safe_get_env("MAX_SEQUENCES_PER_REQUEST", "100"))
)

SCRIPT_CONFIG <- list(
  rf_script_path = file.path(APP_ROOT, "RF-AmPEP30.R"),
  cnn_script_path = file.path(APP_ROOT, "Deep-AmPEP30.R")
)
