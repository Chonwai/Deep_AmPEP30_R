library(plumber)
library(jsonlite)
library(utils)

# Source configuration and dependencies
source(file.path("..", "config", "config.R"))
source(file.path("..", "src", "utils", "helpers.R"))
source(file.path("..", "src", "validation", "fasta_validator.R"))
source(file.path("..", "src", "models", "ampep_predictor.R"))

# Global error handler for graceful API responses
handle_api_error <- function(error, res, code = 500) {
  log_message(paste("API Error:", error$message), level = "ERROR")
  res$status <- code
  return(list(
    status = "error",
    message = error$message,
    code = code,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ))
}

#* Basic health check
#* @get /health
#* @serializer json
function() {
  list(
    status = "healthy",
    service = SERVICE_CONFIG$name,
    version = SERVICE_CONFIG$version,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

#* Detailed health check with model status
#* @get /health/detailed
#* @serializer json
function(res) {
  tryCatch(
    {
      model_status <- get_model_health_status()
      list(
        status = model_status$overall_status,
        service = SERVICE_CONFIG$name,
        version = SERVICE_CONFIG$version,
        models = list(
          rf_predictor = model_status$rf_predictor,
          cnn_predictor = model_status$cnn_predictor
        ),
        timestamp = model_status$timestamp
      )
    },
    error = function(e) {
      handle_api_error(e, res, 503)
    }
  )
}

#* Predict AMPs from FASTA string (supports both RF and CNN methods)
#* @post /api/predict
#* @serializer json
function(req, res) {
  request_id <- paste0("req_", as.integer(Sys.time()), "_", sample(1000:9999, 1))
  log_message(paste("Processing request:", request_id), level = "INFO")

  tryCatch(
    {
      # Request size validation
      body <- req$postBody
      if (!limit_request_size(body, API_CONFIG$max_request_size_bytes)) {
        return(handle_api_error(list(message = "Request entity too large"), res, 413))
      }

      # JSON parsing
      payload <- tryCatch(jsonlite::fromJSON(body), error = function(e) NULL)
      if (is.null(payload)) {
        return(handle_api_error(list(message = "Invalid JSON payload"), res, 400))
      }

      # Extract parameters
      fasta <- payload$fasta
      method <- if (!is.null(payload$method)) tolower(payload$method) else SERVICE_CONFIG$default_method

      if (is.null(fasta) || nchar(fasta) == 0) {
        return(handle_api_error(list(message = "FASTA content is required"), res, 400))
      }

      # Process input
      fasta <- gsub("\\\\n", "\n", fasta)
      fasta <- sanitize_input(fasta)

      # Validation
      v <- validate_fasta_string(
        fasta,
        min_len = VALIDATION_CONFIG$min_sequence_length,
        max_len = VALIDATION_CONFIG$max_sequence_length,
        allowed_chars = VALIDATION_CONFIG$allowed_amino_acids
      )
      if (!isTRUE(v$valid)) {
        return(handle_api_error(list(message = v$reason), res, 400))
      }

      # Prediction
      start_time <- Sys.time()
      log_message(paste("Starting prediction with method:", method, "for request:", request_id), level = "INFO")

      df <- predict_by_method(fasta, method = method)

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      log_message(paste("Prediction completed for request:", request_id, "in", round(elapsed, 3), "seconds"), level = "INFO")

      return(list(
        status = "success",
        message = "Prediction completed successfully",
        data = list(results = df),
        metadata = list(
          processing_time_seconds = round(elapsed, 3),
          sequences_processed = nrow(df),
          method = method,
          request_id = request_id
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    },
    error = function(e) {
      log_message(paste("Request", request_id, "failed:", e$message), level = "ERROR")
      return(handle_api_error(e, res, 500))
    }
  )
}

#* Predict AMPs using Random Forest method specifically
#* @post /api/predict/rf
#* @serializer json
function(req, res) {
  req_copy <- req
  req_copy$postBody <- gsub('"method"\\s*:\\s*"[^"]*"', '"method":"rf"', req$postBody)
  if (!grepl('"method"', req_copy$postBody)) {
    payload <- jsonlite::fromJSON(req$postBody)
    payload$method <- "rf"
    req_copy$postBody <- jsonlite::toJSON(payload, auto_unbox = TRUE)
  }

  # Call the main predict function with RF method
  predict_api_main(req_copy, res)
}

#* Predict AMPs using Deep Learning/CNN method specifically
#* @post /api/predict/deep
#* @serializer json
function(req, res) {
  req_copy <- req
  req_copy$postBody <- gsub('"method"\\s*:\\s*"[^"]*"', '"method":"cnn"', req$postBody)
  if (!grepl('"method"', req_copy$postBody)) {
    payload <- jsonlite::fromJSON(req$postBody)
    payload$method <- "cnn"
    req_copy$postBody <- jsonlite::toJSON(payload, auto_unbox = TRUE)
  }

  # Call the main predict function with CNN method
  predict_api_main(req_copy, res)
}

# Helper function to avoid code duplication
predict_api_main <- function(req, res) {
  request_id <- paste0("req_", as.integer(Sys.time()), "_", sample(1000:9999, 1))
  log_message(paste("Processing request:", request_id), level = "INFO")

  tryCatch(
    {
      # Request size validation
      body <- req$postBody
      if (!limit_request_size(body, API_CONFIG$max_request_size_bytes)) {
        return(handle_api_error(list(message = "Request entity too large"), res, 413))
      }

      # JSON parsing
      payload <- tryCatch(jsonlite::fromJSON(body), error = function(e) NULL)
      if (is.null(payload)) {
        return(handle_api_error(list(message = "Invalid JSON payload"), res, 400))
      }

      # Extract parameters
      fasta <- payload$fasta
      method <- if (!is.null(payload$method)) tolower(payload$method) else SERVICE_CONFIG$default_method

      if (is.null(fasta) || nchar(fasta) == 0) {
        return(handle_api_error(list(message = "FASTA content is required"), res, 400))
      }

      # Process input
      fasta <- gsub("\\\\n", "\n", fasta)
      fasta <- sanitize_input(fasta)

      # Validation
      v <- validate_fasta_string(
        fasta,
        min_len = VALIDATION_CONFIG$min_sequence_length,
        max_len = VALIDATION_CONFIG$max_sequence_length,
        allowed_chars = VALIDATION_CONFIG$allowed_amino_acids
      )
      if (!isTRUE(v$valid)) {
        return(handle_api_error(list(message = v$reason), res, 400))
      }

      # Prediction
      start_time <- Sys.time()
      log_message(paste("Starting prediction with method:", method, "for request:", request_id), level = "INFO")

      df <- predict_by_method(fasta, method = method)

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      log_message(paste("Prediction completed for request:", request_id, "in", round(elapsed, 3), "seconds"), level = "INFO")

      return(list(
        status = "success",
        message = "Prediction completed successfully",
        data = list(results = df),
        metadata = list(
          processing_time_seconds = round(elapsed, 3),
          sequences_processed = nrow(df),
          method = method,
          request_id = request_id
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    },
    error = function(e) {
      log_message(paste("Request", request_id, "failed:", e$message), level = "ERROR")
      return(handle_api_error(e, res, 500))
    }
  )
}
