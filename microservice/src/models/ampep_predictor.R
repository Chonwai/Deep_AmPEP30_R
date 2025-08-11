# Enhanced AmPEP Predictor with industry-standard practices
# This file is sourced from microservice/api/, so paths are relative to that directory

# Ensure dependencies are loaded
if (!exists("SERVICE_CONFIG")) {
  source(file.path("..", "config", "config.R"))
}
if (!exists("log_message")) {
  source(file.path("..", "src", "utils", "helpers.R"))
}

# Source enhanced predictors
source(file.path("..", "src", "models", "rf_predictor_enhanced.R"))
source(file.path("..", "src", "models", "deep_predictor_enhanced.R"))

# Legacy model loading functions (kept for backward compatibility)
load_rf_model <- function() {
  load_rf_model_cached()
}

load_cnn_model <- function() {
  load_cnn_model_cached()
}

# Enhanced prediction functions using new predictors
predict_rf_from_fasta_string <- function(fasta_content, temp_dir = tempdir()) {
  log_message("Using enhanced RF predictor", level = "INFO")

  tryCatch(
    {
      result <- predict_rf_enhanced(fasta_content)
      # Filter out invalid sequences for backward compatibility
      valid_results <- result[result$status == "success", ]
      if (nrow(valid_results) == 0) {
        stop("No valid sequences could be processed")
      }
      return(valid_results[, c("sequence_name", "prediction", "probability")])
    },
    error = function(e) {
      log_message("Enhanced RF prediction failed, falling back to legacy method", level = "WARNING")

      # Legacy fallback
      predict_rf_legacy(fasta_content, temp_dir)
    }
  )
}

predict_cnn_from_fasta_string <- function(fasta_content, temp_dir = tempdir()) {
  log_message("Using enhanced CNN predictor", level = "INFO")

  tryCatch(
    {
      result <- predict_cnn_enhanced(fasta_content)
      # Filter out invalid sequences for backward compatibility
      valid_results <- result[result$status == "success", ]
      if (nrow(valid_results) == 0) {
        stop("No valid sequences could be processed")
      }
      return(valid_results[, c("sequence_name", "prediction", "probability")])
    },
    error = function(e) {
      log_message("Enhanced CNN prediction failed, falling back to legacy method", level = "WARNING")

      # Legacy fallback
      predict_cnn_legacy(fasta_content, temp_dir)
    }
  )
}

# Legacy prediction functions (for fallback)
predict_rf_legacy <- function(fasta_content, temp_dir = tempdir()) {
  tmp_fasta <- file.path(temp_dir, paste0("input_", as.integer(Sys.time()), ".fasta"))
  writeLines(fasta_content, tmp_fasta, useBytes = TRUE)

  script_path <- SCRIPT_CONFIG$rf_script_path
  if (!file.exists(script_path)) stop(paste("RF script not found at", script_path))

  tmp_out <- file.path(temp_dir, paste0("output_", as.integer(Sys.time()), ".out"))
  args <- c(script_path, tmp_fasta, tmp_out)

  old_wd <- getwd()
  setwd(dirname(script_path))

  status <- system2("Rscript", args = c(basename(script_path), tmp_fasta, tmp_out), stdout = TRUE, stderr = TRUE)
  setwd(old_wd)

  if (!file.exists(tmp_out)) {
    log_message("RF prediction failed", level = "ERROR", fields = list(
      stderr = paste(status, collapse = "\n"),
      script_path = script_path
    ))
    stop(paste("RF prediction failed. Status:", paste(status, collapse = "\n")))
  }

  df <- tryCatch(
    {
      read.table(tmp_out, header = FALSE, stringsAsFactors = FALSE)
    },
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 3) stop("Invalid RF output format")
  colnames(df) <- c("sequence_name", "prediction", "probability")
  df$probability <- as.numeric(df$probability)
  df$prediction <- as.integer(df$prediction)
  df
}

predict_cnn_legacy <- function(fasta_content, temp_dir = tempdir()) {
  tmp_fasta <- file.path(temp_dir, paste0("input_", as.integer(Sys.time()), ".fasta"))
  writeLines(fasta_content, tmp_fasta, useBytes = TRUE)

  script_path <- SCRIPT_CONFIG$cnn_script_path
  if (!file.exists(script_path)) stop(paste("CNN script not found at", script_path))
  tmp_out <- file.path(temp_dir, paste0("output_", as.integer(Sys.time()), ".out"))
  args <- c(script_path, tmp_fasta, tmp_out)
  status <- system2("Rscript", args = args, stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp_out)) {
    log_message("CNN prediction failed", level = "ERROR", fields = list(stderr = status))
    stop("CNN prediction failed")
  }
  df <- tryCatch(
    {
      read.table(tmp_out, header = FALSE, stringsAsFactors = FALSE)
    },
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 3) stop("Invalid CNN output format")
  colnames(df) <- c("sequence_name", "prediction", "probability")
  df$probability <- as.numeric(df$probability)
  df$prediction <- as.integer(df$prediction)
  df
}

predict_by_method <- function(fasta_content, method = SERVICE_CONFIG$default_method) {
  m <- tolower(method)

  log_message(paste("Starting prediction with method:", m), level = "INFO")

  tryCatch(
    {
      if (m == "rf") {
        if (!MODEL_CONFIG$enable_rf) stop("RF method disabled")
        result <- predict_rf_from_fasta_string(fasta_content)
        log_message("RF prediction completed successfully", level = "INFO")
        return(result)
      } else if (m == "cnn" || m == "deep") {
        if (!MODEL_CONFIG$enable_cnn) stop("CNN method disabled")
        result <- predict_cnn_from_fasta_string(fasta_content)
        log_message("CNN prediction completed successfully", level = "INFO")
        return(result)
      } else if (m == "auto") {
        # Prefer RF; fallback CNN if RF fails and enabled
        if (MODEL_CONFIG$enable_rf) {
          result <- predict_rf_from_fasta_string(fasta_content)
          log_message("Auto prediction completed using RF", level = "INFO")
          return(result)
        } else if (MODEL_CONFIG$enable_cnn) {
          result <- predict_cnn_from_fasta_string(fasta_content)
          log_message("Auto prediction completed using CNN", level = "INFO")
          return(result)
        } else {
          stop("No methods enabled")
        }
      }

      stop(paste("Unsupported method:", method))
    },
    error = function(e) {
      log_message(paste("Prediction failed with method", m), level = "ERROR", fields = list(error = e$message))
      stop(e$message)
    }
  )
}

# Enhanced health check functions
get_model_health_status <- function() {
  rf_status <- if (MODEL_CONFIG$enable_rf) rf_health_check() else list(status = "disabled")
  cnn_status <- if (MODEL_CONFIG$enable_cnn) cnn_health_check() else list(status = "disabled")

  overall_status <- "healthy"
  if (MODEL_CONFIG$enable_rf && rf_status$status != "healthy") overall_status <- "unhealthy"
  if (MODEL_CONFIG$enable_cnn && cnn_status$status != "healthy") overall_status <- "unhealthy"

  return(list(
    overall_status = overall_status,
    rf_predictor = rf_status,
    cnn_predictor = cnn_status,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ))
}

# Cleanup all resources
cleanup_all_resources <- function() {
  cleanup_rf_resources()
  cleanup_cnn_resources()
  log_message("All predictor resources cleaned up", level = "INFO")
}
