#!/usr/bin/env Rscript

# Format all R/Rmd files in the project using styler (Prettier-like for R)
# Usage:
#   Rscript format_project.R --write            # 寫入格式化（預設）
#   Rscript format_project.R --check            # 僅檢查，如需要格式化則以非零碼結束（適用 CI）
#   Rscript format_project.R --path /path/to/dir

suppressWarnings(suppressMessages({
  # Lazy-install styler if missing
  if (!requireNamespace("styler", quietly = TRUE)) {
    message("Installing 'styler' package ...")
    install.packages("styler", repos = "https://cran.r-project.org")
  }
}))

args <- commandArgs(trailingOnly = TRUE)

mode <- "write" # write | check
target_path <- getwd() # default: current working directory

if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--check")) {
      mode <- "check"
      i <- i + 1
    } else if (identical(arg, "--write")) {
      mode <- "write"
      i <- i + 1
    } else if (identical(arg, "--path")) {
      if (i == length(args)) stop("--path 需要提供一個目錄路徑")
      target_path <- args[[i + 1]]
      i <- i + 2
    } else if (grepl("^-", arg)) {
      stop(sprintf("未知參數: %s", arg))
    } else {
      # positional path
      target_path <- arg
      i <- i + 1
    }
  }
}

if (!dir.exists(target_path)) {
  stop(sprintf("指定的目錄不存在: %s", target_path))
}

cat(sprintf("[format] path = %s\n", normalizePath(target_path, winslash = "/", mustWork = FALSE)))
cat(sprintf("[format] mode = %s\n", mode))

# Map mode to styler dry-run behavior
dry_mode <- switch(mode,
  check = "fail", # error if changes would be made
  write = "off", # write changes
  stop(sprintf("未知模式: %s", mode))
)

# Configure file types we want to format (R and Rmd)
file_types <- c("R", "Rmd")

format_call <- function() {
  # Use styler's recursive directory formatter. It respects .stylerignore if present.
  styler::style_dir(
    path = target_path,
    dry = dry_mode,
    filetype = file_types
    # Default style is tidyverse; customize via styler::tidyverse_style() if needed.
  )
}

result <- try(format_call(), silent = TRUE)

if (inherits(result, "try-error")) {
  msg <- conditionMessage(attr(result, "condition"))
  # If mode == check and changes are needed, styler errors on dry = "fail"
  if (identical(mode, "check")) {
    cat("[format] 需要格式化變更。\n")
    cat(sprintf("[format] styler 訊息: %s\n", msg))
    quit(save = "no", status = 1)
  } else {
    cat("[format] 發生錯誤：\n")
    cat(sprintf("%s\n", msg))
    quit(save = "no", status = 2)
  }
} else {
  if (identical(mode, "check")) {
    cat("[format] 所有檔案皆已符合格式。\n")
  } else {
    cat("[format] 格式化完成。\n")
  }
  quit(save = "no", status = 0)
}
