# Enhanced RF Predictor with embedded functionality
# Industry-standard error handling and logging

require(randomForest)
require(caret)
require(seqinr)
require(protr)

# Global variables for model caching
.rf_model_cache <- NULL
.rf_model_loaded <- FALSE

# Enhanced logging function
log_rf_operation <- function(message, level = "INFO", error = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [RF-%s] %s", timestamp, level, message)
  if (!is.null(error)) {
    log_entry <- paste(log_entry, "Error:", error$message)
  }
  cat(log_entry, "\n", file = stderr())
}

# Load and cache RF model
load_rf_model_cached <- function() {
  if (.rf_model_loaded && !is.null(.rf_model_cache)) {
    log_rf_operation("Using cached RF model")
    return(.rf_model_cache)
  }

  model_path <- MODEL_CONFIG$rf_model_path
  log_rf_operation(paste("Loading RF model from:", model_path))

  tryCatch(
    {
      if (!file.exists(model_path)) {
        stop(paste("RF model file not found at:", model_path))
      }

      .rf_model_cache <<- readRDS(model_path)
      .rf_model_loaded <<- TRUE
      log_rf_operation("RF model loaded successfully")
      return(.rf_model_cache)
    },
    error = function(e) {
      log_rf_operation("Failed to load RF model", "ERROR", e)
      stop(paste("RF model loading failed:", e$message))
    }
  )
}

# Enhanced sequence processing with better error handling
process_fasta_sequences <- function(fasta_content) {
  tryCatch(
    {
      # Create temporary file
      temp_file <- tempfile(pattern = "rf_input_", fileext = ".fasta")
      on.exit(unlink(temp_file))

      # Write FASTA content
      writeLines(fasta_content, temp_file, useBytes = TRUE)

      # Read sequences using seqinr
      protdata <- read.fasta(temp_file, seqtype = "AA", as.string = TRUE)

      if (length(protdata) == 0) {
        stop("No valid sequences found in FASTA input")
      }

      names_list <- getName(protdata)
      sequences <- unlist(getSequence(protdata, as.string = TRUE), recursive = FALSE)

      # Filter sequences by length (5-30 amino acids)
      seq_lengths <- nchar(sequences)
      valid_indices <- which(seq_lengths >= 5 & seq_lengths <= 30)

      if (length(valid_indices) == 0) {
        warning("No sequences meet length requirements (5-30 amino acids)")
        return(list(
          sequences = character(0),
          names = character(0),
          invalid_sequences = sequences,
          invalid_names = names_list
        ))
      }

      invalid_indices <- setdiff(1:length(sequences), valid_indices)

      return(list(
        sequences = sequences[valid_indices],
        names = names_list[valid_indices],
        invalid_sequences = if (length(invalid_indices) > 0) sequences[invalid_indices] else character(0),
        invalid_names = if (length(invalid_indices) > 0) names_list[invalid_indices] else character(0)
      ))
    },
    error = function(e) {
      log_rf_operation("Error processing FASTA sequences", "ERROR", e)
      stop(paste("FASTA processing failed:", e$message))
    }
  )
}

# Enhanced feature extraction (from RF-AmPEP30.R)
extract_amino_acid_composition <- function(sequences) {
  tryCatch(
    {
      log_rf_operation(paste("Extracting features for", length(sequences), "sequences"))

      feature_matrix <- NULL

      for (i in seq_along(sequences)) {
        seq <- sequences[i]

        # Extract amino acid composition using protr
        aac_features <- extractAAC(seq)

        # Normalize by sequence length (as in original script)
        aac_features <- aac_features * nchar(seq)

        feature_matrix <- rbind(feature_matrix, aac_features)
      }

      rownames(feature_matrix) <- names(sequences)
      log_rf_operation("Feature extraction completed successfully")

      return(data.frame(feature_matrix))
    },
    error = function(e) {
      log_rf_operation("Error in feature extraction", "ERROR", e)
      stop(paste("Feature extraction failed:", e$message))
    }
  )
}

# Make predictions with comprehensive error handling
predict_rf_enhanced <- function(fasta_content) {
  start_time <- Sys.time()
  log_rf_operation("Starting RF prediction")

  tryCatch(
    {
      # Process input sequences
      seq_data <- process_fasta_sequences(fasta_content)

      if (length(seq_data$sequences) == 0) {
        log_rf_operation("No valid sequences to process", "WARNING")
        return(data.frame(
          sequence_name = character(0),
          prediction = integer(0),
          probability = numeric(0),
          sequence = character(0),
          status = character(0)
        ))
      }

      # Load model
      rf_model <- load_rf_model_cached()

      # Extract features
      features <- extract_amino_acid_composition(seq_data$sequences)

      # Make predictions
      log_rf_operation("Making predictions with RF model")
      predictions <- predict(rf_model, newdata = features, type = "prob")

      # Process results
      results <- data.frame(
        sequence_name = seq_data$names,
        prediction = ifelse(predictions[, 2] >= 0.5, 1, 0),
        probability = round(predictions[, 2], 6),
        sequence = seq_data$sequences,
        status = "success"
      )

      # Add invalid sequences with status
      if (length(seq_data$invalid_sequences) > 0) {
        invalid_results <- data.frame(
          sequence_name = seq_data$invalid_names,
          prediction = -1,
          probability = -1,
          sequence = seq_data$invalid_sequences,
          status = "invalid_length"
        )
        results <- rbind(results, invalid_results)
      }

      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      log_rf_operation(paste("RF prediction completed in", round(elapsed_time, 3), "seconds"))

      return(results)
    },
    error = function(e) {
      log_rf_operation("RF prediction failed", "ERROR", e)
      stop(paste("RF prediction error:", e$message))
    }
  )
}

# Cleanup function for memory management
cleanup_rf_resources <- function() {
  tryCatch(
    {
      if (exists(".rf_model_cache") && !is.null(.rf_model_cache)) {
        .rf_model_cache <<- NULL
        .rf_model_loaded <<- FALSE
        gc() # Force garbage collection
        log_rf_operation("RF resources cleaned up")
      }
    },
    error = function(e) {
      log_rf_operation("Error during cleanup", "WARNING", e)
    }
  )
}

# Health check for RF predictor
rf_health_check <- function() {
  tryCatch(
    {
      # Try to load model
      model <- load_rf_model_cached()

      # Test with dummy sequence
      test_fasta <- ">test_sequence\nACDEFGHIKLMNPQRSTVWY"
      result <- predict_rf_enhanced(test_fasta)

      return(list(
        status = "healthy",
        model_loaded = !is.null(model),
        test_prediction = !is.null(result) && nrow(result) > 0
      ))
    },
    error = function(e) {
      return(list(
        status = "unhealthy",
        error = e$message,
        model_loaded = FALSE,
        test_prediction = FALSE
      ))
    }
  )
}
