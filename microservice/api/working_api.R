#!/usr/bin/env Rscript

# 工作版本的 AmPEP30 API 服務
cat("啟動 AmPEP30 API 服務...\n")

# 設定工作目錄
setwd("/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30")

# 載入必要套件
suppressWarnings(suppressMessages({
  library(plumber)
  library(jsonlite)
  library(seqinr)
  library(caret)
  library(randomForest)
  library(protr)
}))

cat("✓ 套件載入完成\n")

# 檢查 RF 模型
rf_model_path <- "AmPEP30-RF-1200tree.mdl"
if (!file.exists(rf_model_path)) {
  stop("RF 模型文件不存在: ", rf_model_path)
}
cat("✓ RF 模型文件確認存在\n")

# 預測函數
predict_amp <- function(fasta_content) {
  tryCatch(
    {
      # 創建臨時文件
      temp_id <- paste0("temp_", as.integer(Sys.time()), "_", sample(1000:9999, 1))
      temp_fasta <- paste0(temp_id, ".fasta")
      temp_output <- paste0(temp_id, ".out")

      # 寫入 fasta 內容
      writeLines(fasta_content, temp_fasta)

      # 執行預測
      cmd <- paste("Rscript RF-AmPEP30.R", temp_fasta, temp_output)
      system_result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

      if (system_result != 0) {
        unlink(c(temp_fasta, temp_output))
        return(list(success = FALSE, error = "預測執行失敗"))
      }

      # 讀取結果
      if (!file.exists(temp_output)) {
        unlink(c(temp_fasta, temp_output))
        return(list(success = FALSE, error = "預測結果文件未生成"))
      }

      results <- read.table(temp_output, header = FALSE, stringsAsFactors = FALSE, sep = " ")

      # 清理臨時文件
      unlink(c(temp_fasta, temp_output))

      # 格式化結果
      predictions <- list()
      for (i in 1:nrow(results)) {
        seq_name <- results[i, 1]
        prediction <- results[i, 2]
        probability <- as.numeric(results[i, 3])

        # 保持數字格式：1=AMP, 0=Non-AMP, -1=Invalid/Error
        predictions[[i]] <- list(
          sequence_name = seq_name,
          prediction = prediction, # 直接使用數字
          amp_probability = probability,
          confidence = if (prediction == 1) probability else if (prediction == 0) (1 - probability) else -1
        )
      }

      return(list(success = TRUE, predictions = predictions))
    },
    error = function(e) {
      return(list(success = FALSE, error = paste("預測錯誤:", e$message)))
    }
  )
}

# 創建 plumber API
#* @apiTitle AmPEP30 預測 API
#* @apiDescription 使用隨機森林模型預測抗菌胜肽 (AMP)

#* 健康檢查
#* @get /health
function() {
  list(
    status = "healthy",
    service = "AmPEP30 API",
    version = "1.0",
    timestamp = Sys.time()
  )
}

#* AMP 預測 - 單一序列
#* @param sequence:str 蛋白質序列 (5-30 氨基酸)
#* @post /predict/single
function(sequence = "") {
  if (sequence == "") {
    return(list(success = FALSE, error = "請提供蛋白質序列"))
  }

  # 檢查序列長度
  seq_length <- nchar(sequence)
  if (seq_length < 5 || seq_length > 30) {
    return(list(
      success = FALSE,
      error = paste0("序列長度必須在 5-30 氨基酸之間，當前長度: ", seq_length)
    ))
  }

  # 檢查序列字符
  valid_chars <- "ACDEFGHIKLMNPQRSTVWY"
  if (!all(strsplit(toupper(sequence), "")[[1]] %in% strsplit(valid_chars, "")[[1]])) {
    return(list(success = FALSE, error = "序列包含無效氨基酸字符"))
  }

  # 創建 FASTA 格式
  fasta_content <- paste0(">input_sequence\n", toupper(sequence))

  # 執行預測
  result <- predict_amp(fasta_content)

  if (result$success && length(result$predictions) > 0) {
    pred <- result$predictions[[1]]
    return(list(
      success = TRUE,
      sequence = sequence,
      length = seq_length,
      prediction = pred$prediction,
      amp_probability = pred$amp_probability,
      confidence = pred$confidence
    ))
  } else {
    return(result)
  }
}

#* AMP 預測 - FASTA 格式
#* @param fasta:str FASTA 格式的序列
#* @post /predict/fasta
function(fasta = "") {
  if (fasta == "") {
    return(list(success = FALSE, error = "請提供 FASTA 格式序列"))
  }

  # 驗證 FASTA 格式
  lines <- strsplit(fasta, "\n")[[1]]
  lines <- lines[lines != ""] # 移除空行

  if (length(lines) < 2 || !grepl("^>", lines[1])) {
    return(list(success = FALSE, error = "無效的 FASTA 格式"))
  }

  # 執行預測
  result <- predict_amp(fasta)
  return(result)
}

#* 取得模型資訊
#* @get /info
function() {
  list(
    model_type = "Random Forest",
    model_file = rf_model_path,
    sequence_length_range = "5-30 amino acids",
    amino_acids = "ACDEFGHIKLMNPQRSTVWY",
    description = "AmPEP30 antimicrobial peptide prediction model"
  )
}

# 啟動 API 服務
cat("✓ API 路由設定完成\n")
cat("啟動服務於 http://localhost:8005\n")
cat("API 文檔: http://localhost:8005/__docs__/\n")
cat("健康檢查: http://localhost:8005/health\n")
cat("按 Ctrl+C 停止服務\n\n")

# 啟動 plumber
pr() %>%
  pr_set_api_spec(function(spec) {
    spec$info$title <- "AmPEP30 預測 API"
    spec$info$description <- "抗菌胜肽 (AMP) 預測服務"
    spec
  }) %>%
  pr_run(host = "0.0.0.0", port = 8005)
