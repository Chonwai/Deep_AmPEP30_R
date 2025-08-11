#!/usr/bin/env Rscript

# Test RF functionality only
library(plumber)
library(jsonlite)

cat("Loading RF-only test API...\n")

# Source only necessary files
source(file.path("..", "config", "config.R"))
source(file.path("..", "src", "utils", "helpers.R"))
source(file.path("..", "src", "validation", "fasta_validator.R"))

# Source the original RF script for testing
source(file.path("..", "..", "RF-AmPEP30.R"))

pr <- pr() %>%
  pr_get("/health", function() {
    list(
      status = "healthy",
      service = "RF-AmPEP30-Test",
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  }) %>%
  pr_post("/api/predict/rf", function(req, res) {
    tryCatch(
      {
        body <- req$postBody
        payload <- jsonlite::fromJSON(body)
        fasta <- payload$fasta

        if (is.null(fasta) || nchar(fasta) == 0) {
          res$status <- 400
          return(list(status = "error", message = "FASTA content required"))
        }

        # Process input
        fasta <- gsub("\\\\n", "\n", fasta)
        fasta <- sanitize_input(fasta)

        # Simple validation
        v <- validate_fasta_string(
          fasta,
          min_len = 5,
          max_len = 30,
          allowed_chars = "ACDEFGHIKLMNPQRSTVWY"
        )

        if (!isTRUE(v$valid)) {
          res$status <- 400
          return(list(status = "error", message = v$reason))
        }

        # Use original RF script function
        start_time <- Sys.time()

        # Create temp files
        temp_fasta <- tempfile(fileext = ".fasta")
        temp_out <- tempfile(fileext = ".out")

        writeLines(fasta, temp_fasta)

        # Call RF prediction function directly
        result <- tryCatch(
          {
            predict_from_fasta_rf(temp_fasta, temp_out)
            read.table(temp_out, header = FALSE, stringsAsFactors = FALSE)
          },
          error = function(e) {
            list(error = e$message)
          }
        )

        # Cleanup
        unlink(c(temp_fasta, temp_out))

        if (is.list(result) && !is.null(result$error)) {
          res$status <- 500
          return(list(status = "error", message = result$error))
        }

        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

        if (ncol(result) >= 3) {
          colnames(result) <- c("sequence_name", "prediction", "probability")
          result$probability <- as.numeric(result$probability)
          result$prediction <- as.integer(result$prediction)
        }

        return(list(
          status = "success",
          message = "RF prediction completed",
          data = list(results = result),
          metadata = list(
            processing_time_seconds = round(elapsed, 3),
            sequences_processed = nrow(result),
            method = "rf"
          )
        ))
      },
      error = function(e) {
        res$status <- 500
        return(list(status = "error", message = e$message))
      }
    )
  })

cat("Starting RF test API on port 8003...\n")
pr_run(pr, host = "127.0.0.1", port = 8003)
