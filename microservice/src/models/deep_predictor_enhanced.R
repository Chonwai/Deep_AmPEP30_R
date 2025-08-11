# Enhanced Deep Learning Predictor with embedded functionality
# Industry-standard error handling and logging

require(reticulate)

# Ensure reticulate uses the intended Python before loading keras
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
  Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")
}
try(
  {
    reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)
  },
  silent = TRUE
)

require(keras)
require(seqinr)
require(protr)
require(caret)

# Global variables for model caching
.cnn_model_cache <- NULL
.cnn_model_loaded <- FALSE

# Enhanced logging function
log_deep_operation <- function(message, level = "INFO", error = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [DEEP-%s] %s", timestamp, level, message)
  if (!is.null(error)) {
    log_entry <- paste(log_entry, "Error:", error$message)
  }
  cat(log_entry, "\n", file = stderr())
}

# Load and cache CNN model
load_cnn_model_cached <- function() {
  if (.cnn_model_loaded && !is.null(.cnn_model_cache)) {
    log_deep_operation("Using cached CNN model")
    return(.cnn_model_cache)
  }

  model_path <- MODEL_CONFIG$cnn_model_path
  log_deep_operation(paste("Loading CNN model from:", model_path))

  tryCatch(
    {
      if (!file.exists(model_path)) {
        stop(paste("CNN model file not found at:", model_path))
      }

      # Validate Keras/TensorFlow availability (do not install at runtime)
      if (!keras::is_keras_available()) {
        cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
        py_path <- if (!is.null(cfg)) cfg$python else "<unknown>"
        stop(paste(
          "Keras/TensorFlow backend is not available.",
          "Python:", py_path,
          "RETICULATE_PYTHON=", Sys.getenv("RETICULATE_PYTHON", unset = ""),
          "Please ensure tensorflow==2.13.x and h5py are installed in this Python."
        ))
      }

      .cnn_model_cache <<- keras::load_model_hdf5(model_path)
      .cnn_model_loaded <<- TRUE
      log_deep_operation("CNN model loaded successfully")
      return(.cnn_model_cache)
    },
    error = function(e) {
      log_deep_operation("Failed to load CNN model", "ERROR", e)
      stop(paste("CNN model loading failed:", e$message))
    }
  )
}

# Feature extraction for CNN (based on Deep-AmPEP30.R)
# Define amino acid clustering types from original script
get_clustering_types <- function() {
  list(
    type8_cluster17 = c("AT", "C", "DE", "F", "G", "H", "IV", "K", "L", "M", "N", "P", "Q", "R", "S", "V", "W"),
    type3a_cluster19 = c("FA", "P", "G", "S", "T", "D", "E", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "V", "C"),
    type12_cluster18 = c("TVLI", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "D", "E"),
    type7_cluster15 = c("C", "K", "R", "W", "Y", "A", "FILV", "M", "D", "E", "Q", "H", "TP", "GS", "N"),
    type12_cluster17 = c("TVLI", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "DE")
  )
}

# Extract PseKRAAC features (from Deep-AmPEP30.R)
extract_psekraac_features <- function(sequences) {
  tryCatch(
    {
      log_deep_operation(paste("Extracting PseKRAAC features for", length(sequences), "sequences"))

      # Extract basic amino acid composition
      aac_matrix <- NULL
      for (i in seq_along(sequences)) {
        seq <- sequences[i]
        aac_features <- extractAAC(seq)
        aac_features <- aac_features * nchar(seq) # Normalize by length
        aac_matrix <- rbind(aac_matrix, aac_features)
      }
      rownames(aac_matrix) <- names(sequences)
      aac_df <- data.frame(aac_matrix)

      # Apply clustering transformations
      clustering_types <- get_clustering_types()
      descriptors <- unlist(clustering_types, use.names = FALSE)

      feature_matrix <- NULL
      for (descriptor in descriptors) {
        if (nchar(descriptor) > 1) {
          # Multi-character clusters - sum corresponding amino acids
          chars <- strsplit(descriptor, "")[[1]]
          if (all(chars %in% colnames(aac_df))) {
            cluster_sum <- apply(aac_df[, chars, drop = FALSE], 1, sum)
            feature_matrix <- cbind(feature_matrix, cluster_sum)
          }
        } else {
          # Single character
          if (descriptor %in% colnames(aac_df)) {
            feature_matrix <- cbind(feature_matrix, aac_df[, descriptor])
          }
        }
      }

      # Set column names
      colnames(feature_matrix) <- paste0("feature_", 1:ncol(feature_matrix))

      log_deep_operation("PseKRAAC feature extraction completed successfully")
      return(data.frame(feature_matrix))
    },
    error = function(e) {
      log_deep_operation("Error in PseKRAAC feature extraction", "ERROR", e)
      stop(paste("Feature extraction failed:", e$message))
    }
  )
}

# Process sequences for CNN input
process_sequences_for_cnn <- function(fasta_content) {
  tryCatch(
    {
      # Create temporary file
      temp_file <- tempfile(pattern = "cnn_input_", fileext = ".fasta")
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
      log_deep_operation("Error processing sequences for CNN", "ERROR", e)
      stop(paste("Sequence processing failed:", e$message))
    }
  )
}

# Make predictions with CNN model
predict_cnn_enhanced <- function(fasta_content) {
  start_time <- Sys.time()
  log_deep_operation("Starting CNN prediction")

  tryCatch(
    {
      # Process input sequences
      seq_data <- process_sequences_for_cnn(fasta_content)

      if (length(seq_data$sequences) == 0) {
        log_deep_operation("No valid sequences to process", "WARNING")
        return(data.frame(
          sequence_name = character(0),
          prediction = integer(0),
          probability = numeric(0),
          sequence = character(0),
          status = character(0)
        ))
      }

      # Load model
      cnn_model <- load_cnn_model_cached()

      # Extract features
      features <- extract_psekraac_features(seq_data$sequences)

      # Prepare data for CNN (reshape to match expected input)
      feature_matrix <- as.matrix(features)
      n_samples <- nrow(feature_matrix)
      n_features <- ncol(feature_matrix)

      # Reshape for CNN input (samples, features, 1)
      cnn_input <- array(feature_matrix, dim = c(n_samples, n_features, 1))

      # Make predictions
      log_deep_operation("Making predictions with CNN model")
      predictions <- predict(cnn_model, cnn_input)

      # Process results
      results <- data.frame(
        sequence_name = seq_data$names,
        prediction = ifelse(predictions >= 0.5, 1, 0),
        probability = round(as.numeric(predictions), 6),
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
      log_deep_operation(paste("CNN prediction completed in", round(elapsed_time, 3), "seconds"))

      return(results)
    },
    error = function(e) {
      log_deep_operation("CNN prediction failed", "ERROR", e)
      stop(paste("CNN prediction error:", e$message))
    }
  )
}

# Cleanup function for memory management
cleanup_cnn_resources <- function() {
  tryCatch(
    {
      if (exists(".cnn_model_cache") && !is.null(.cnn_model_cache)) {
        .cnn_model_cache <<- NULL
        .cnn_model_loaded <<- FALSE
        gc() # Force garbage collection
        log_deep_operation("CNN resources cleaned up")
      }
    },
    error = function(e) {
      log_deep_operation("Error during cleanup", "WARNING", e)
    }
  )
}

# Health check for CNN predictor
cnn_health_check <- function() {
  tryCatch(
    {
      # Check if Keras is available
      if (!keras::is_keras_available()) {
        return(list(
          status = "unhealthy",
          error = "Keras backend not available",
          model_loaded = FALSE,
          test_prediction = FALSE
        ))
      }

      # Try to load model
      model <- load_cnn_model_cached()

      # Test with dummy sequence
      test_fasta <- ">test_sequence\nACDEFGHIKLMNPQRSTVWY"
      result <- predict_cnn_enhanced(test_fasta)

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
