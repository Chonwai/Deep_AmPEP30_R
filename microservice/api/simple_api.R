#!/usr/bin/env Rscript

# 簡化版本的 AmPEP30 API 服務
cat("啟動簡化版 AmPEP30 API...\n")

# 設定工作目錄
setwd("/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30")

# 載入套件
suppressWarnings(suppressMessages({
  library(plumber)
  library(jsonlite)
  library(seqinr)
  library(caret)
  library(randomForest)
  library(protr)
}))

cat("✓ 套件載入完成\n")

# 預測函數
predict_amp <- function(fasta_content) {
  tryCatch(
    {
      temp_id <- as.integer(as.numeric(Sys.time()) * 1000)
      temp_fasta <- paste0("temp_", temp_id, ".fasta")
      temp_output <- paste0("temp_", temp_id, ".out")

      writeLines(fasta_content, temp_fasta)

      cmd <- paste("Rscript RF-AmPEP30.R", temp_fasta, temp_output)
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

      if (file.exists(temp_output)) {
        results <- read.table(temp_output, header = FALSE, stringsAsFactors = FALSE, sep = " ")
        unlink(c(temp_fasta, temp_output))

        predictions <- list()
        for (i in 1:nrow(results)) {
          seq_name <- results[i, 1]
          prediction <- results[i, 2]
          probability <- as.numeric(results[i, 3])

          predictions[[i]] <- list(
            sequence_name = seq_name,
            prediction = ifelse(prediction == 1, "AMP", "Non-AMP"),
            probability = probability
          )
        }

        return(list(success = TRUE, predictions = predictions))
      } else {
        unlink(c(temp_fasta, temp_output))
        return(list(success = FALSE, error = "預測失敗"))
      }
    },
    error = function(e) {
      return(list(success = FALSE, error = paste("錯誤:", e$message)))
    }
  )
}

#* 健康檢查
#* @get /health
function() {
  list(status = "healthy", service = "AmPEP30 API", timestamp = Sys.time())
}

#* AMP 預測
#* @param sequence 蛋白質序列
#* @post /predict
function(sequence = "") {
  if (sequence == "") {
    return(list(success = FALSE, error = "請提供序列"))
  }

  seq_length <- nchar(sequence)
  if (seq_length < 5 || seq_length > 30) {
    return(list(success = FALSE, error = paste0("序列長度必須在 5-30 之間，當前: ", seq_length)))
  }

  fasta_content <- paste0(">input\n", toupper(sequence))
  result <- predict_amp(fasta_content)

  if (result$success && length(result$predictions) > 0) {
    pred <- result$predictions[[1]]
    return(list(
      success = TRUE,
      sequence = sequence,
      prediction = pred$prediction,
      probability = pred$probability
    ))
  } else {
    return(result)
  }
}

#* FASTA 預測
#* @param fasta FASTA 格式序列
#* @post /predict_fasta
function(fasta = "") {
  if (fasta == "") {
    return(list(success = FALSE, error = "請提供 FASTA 序列"))
  }

  result <- predict_amp(fasta)
  return(result)
}

cat("啟動 API 於 port 8006...\n")
