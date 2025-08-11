validate_fasta_string <- function(fasta_content, min_len, max_len, allowed_chars) {
  if (is.null(fasta_content) || !nzchar(fasta_content)) {
    return(list(valid = FALSE, reason = "Empty input"))
  }

  lines <- strsplit(fasta_content, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    return(list(valid = FALSE, reason = "No content"))
  }

  if (!grepl("^>", lines[1])) {
    return(list(valid = FALSE, reason = "First line must start with '>'"))
  }

  # Build sequences
  sequences <- list()
  current_name <- NULL
  current_seq <- ""
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(current_name)) {
        sequences[[length(sequences) + 1]] <- list(name = current_name, sequence = current_seq)
      }
      current_name <- sub("^>", "", line)
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, gsub("\\s+", "", line))
    }
  }
  if (!is.null(current_name)) {
    sequences[[length(sequences) + 1]] <- list(name = current_name, sequence = current_seq)
  }

  if (length(sequences) == 0) {
    return(list(valid = FALSE, reason = "No sequences found"))
  }

  allowed <- strsplit(allowed_chars, "", fixed = TRUE)[[1]]
  for (s in sequences) {
    seq_upper <- toupper(s$sequence)
    chars <- strsplit(seq_upper, "", fixed = TRUE)[[1]]
    if (!all(chars %in% allowed)) {
      return(list(valid = FALSE, reason = paste0("Invalid amino acids in sequence ", s$name)))
    }
    n <- nchar(seq_upper)
    if (n < min_len || n > max_len) {
      return(list(valid = FALSE, reason = paste0("Sequence length out of range for ", s$name)))
    }
  }

  list(valid = TRUE, sequences = sequences)
}
