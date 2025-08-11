library(jsonlite)

log_message <- function(message, level = "INFO", fields = list()) {
  entry <- c(
    list(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      level = level,
      message = as.character(message)
    ),
    fields
  )
  cat(jsonlite::toJSON(entry, auto_unbox = TRUE), "\n")
}

sanitize_input <- function(input_string) {
  if (is.null(input_string)) {
    return("")
  }
  # Don't remove > character as it's needed for FASTA format
  x <- gsub("[<\"'&]", "", input_string)
  Encoding(x) <- "UTF-8"
  x
}

limit_request_size <- function(req_body, max_bytes) {
  nchar(req_body, type = "bytes") <= max_bytes
}

http_error <- function(status, message, code = "ERROR", details = NULL) {
  list(
    status = "error",
    message = message,
    error_code = code,
    details = details,
    http_status = status,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}
