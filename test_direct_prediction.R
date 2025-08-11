#!/usr/bin/env Rscript

# 直接測試預測功能（不使用 API）
cat("=== 直接預測測試 ===\n\n")

setwd("/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30")

# 載入套件
suppressWarnings(suppressMessages({
  library(seqinr)
  library(caret)
  library(randomForest)
  library(protr)
}))

# 直接預測函數
direct_predict <- function(sequence, seq_name = "test") {
  tryCatch(
    {
      # 檢查序列長度
      seq_length <- nchar(sequence)
      if (seq_length < 5 || seq_length > 30) {
        return(list(
          success = FALSE,
          error = paste0("序列長度超出範圍 (5-30)，當前: ", seq_length)
        ))
      }

      # 生成唯一的文件名
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      random_num <- sample(1000:9999, 1)
      temp_fasta <- paste0("temp_", timestamp, "_", random_num, ".fasta")
      temp_output <- paste0("temp_", timestamp, "_", random_num, ".out")

      # 創建 FASTA 內容
      fasta_content <- paste0(">", seq_name, "\n", toupper(sequence))
      writeLines(fasta_content, temp_fasta)

      # 執行預測命令
      cmd <- paste("Rscript RF-AmPEP30.R", temp_fasta, temp_output)

      cat(sprintf("執行命令: %s\n", cmd))
      system_result <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

      if (system_result != 0) {
        unlink(c(temp_fasta, temp_output))
        return(list(success = FALSE, error = paste("命令執行失敗，退出碼:", system_result)))
      }

      # 檢查輸出文件
      if (!file.exists(temp_output)) {
        unlink(c(temp_fasta, temp_output))
        return(list(success = FALSE, error = "輸出文件未生成"))
      }

      # 讀取結果
      results <- read.table(temp_output, header = FALSE, stringsAsFactors = FALSE, sep = " ")

      # 清理文件
      unlink(c(temp_fasta, temp_output))

      if (nrow(results) > 0) {
        seq_name_result <- results[1, 1]
        prediction <- results[1, 2]
        probability <- as.numeric(results[1, 3])

        pred_label <- if (prediction == 1) {
          "AMP"
        } else if (prediction == 0) {
          "Non-AMP"
        } else {
          "Error"
        }

        return(list(
          success = TRUE,
          sequence_name = seq_name_result,
          prediction = pred_label,
          amp_probability = probability,
          confidence = if (prediction == 1) probability else (1 - probability)
        ))
      } else {
        return(list(success = FALSE, error = "結果文件為空"))
      }
    },
    error = function(e) {
      return(list(success = FALSE, error = paste("預測錯誤:", e$message)))
    }
  )
}

# 測試序列
test_cases <- list(
  list(name = "Magainin", seq = "GIGKFLHSAKKFGKAFVGEIMNS", expected = "AMP"),
  list(name = "Melittin", seq = "GIGAVLKVLTTGLPALISWIKRKRQQ", expected = "AMP"),
  list(name = "Short_valid", seq = "KLMNP", expected = "Unknown"),
  list(name = "Medium_test", seq = "ACDEFGHIKL", expected = "Unknown"),
  list(name = "Too_long", seq = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY", expected = "Error"),
  list(name = "Too_short", seq = "ACE", expected = "Error")
)

cat("開始測試...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

for (i in 1:length(test_cases)) {
  test_case <- test_cases[[i]]

  cat(sprintf("測試 %d: %s\n", i, test_case$name))
  cat(sprintf("序列: %s (長度: %d)\n", test_case$seq, nchar(test_case$seq)))
  cat(sprintf("預期: %s\n", test_case$expected))

  result <- direct_predict(test_case$seq, test_case$name)

  if (result$success) {
    cat(sprintf(
      "✓ 結果: %s (AMP 機率: %.3f, 信心度: %.3f)\n",
      result$prediction, result$amp_probability, result$confidence
    ))
  } else {
    cat(sprintf("✗ 錯誤: %s\n", result$error))
  }

  cat(paste(rep("-", 50), collapse = ""), "\n\n")
}

cat("測試完成！\n")
