cat("Checking dependencies...\n")
required_packages <- c("plumber", "jsonlite", "randomForest", "caret", "seqinr", "protr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Missing package:", pkg, "\n")
  } else {
    cat("Found package:", pkg, "\n")
  }
}
