#!/usr/bin/env Rscript

# Simple test API to verify basic functionality
library(plumber)

cat("Loading simple test API...\n")

# Create simple API
api <- function() {
  list(status = "ok", message = "Simple API is working")
}

pr <- pr() %>%
  pr_get("/health", function() {
    list(status = "healthy", timestamp = Sys.time())
  }) %>%
  pr_get("/test", function() {
    list(status = "ok", message = "Simple test endpoint")
  })

cat("Starting simple API on port 8002...\n")
pr_run(pr, host = "127.0.0.1", port = 8002)
