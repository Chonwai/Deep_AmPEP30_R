#!/usr/bin/env Rscript

# 工作演示版本 - 修正路徑問題
cat("=== AmPEP30 工作演示 ===\n\n")

# 設定工作目錄
setwd("/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30")

cat("1. 檢查必要文件...\n")

# 檢查 RF 模型文件
rf_model_path <- "AmPEP30-RF-1200tree.mdl"
if (file.exists(rf_model_path)) {
  cat("✓ RF 模型文件存在:", rf_model_path, "\n")
} else {
  cat("✗ RF 模型文件不存在:", rf_model_path, "\n")
  quit(save = "no", status = 1)
}

# 檢查 RF 腳本
rf_script_path <- "RF-AmPEP30.R"
if (file.exists(rf_script_path)) {
  cat("✓ RF 腳本存在:", rf_script_path, "\n")
} else {
  cat("✗ RF 腳本不存在:", rf_script_path, "\n")
  quit(save = "no", status = 1)
}

cat("\n2. 載入必要套件...\n")

# 載入套件
tryCatch(
  {
    suppressWarnings(suppressMessages({
      library(seqinr)
      library(caret)
      library(randomForest)
      library(protr)
    }))
    cat("✓ 所有套件載入成功\n")
  },
  error = function(e) {
    cat("✗ 套件載入失敗:", e$message, "\n")
    quit(save = "no", status = 1)
  }
)

cat("\n3. 定義簡化預測函數...\n")

# 簡化的預測函數
simple_predict <- function(fasta_content, output_file) {
  tryCatch(
    {
      # 寫入臨時 fasta 文件
      temp_fasta <- "temp_input.fasta"
      writeLines(fasta_content, temp_fasta)

      # 使用 RF 腳本進行預測
      cmd <- paste("Rscript", rf_script_path, temp_fasta, output_file)
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

      # 清理臨時文件
      unlink(temp_fasta)

      return(file.exists(output_file))
    },
    error = function(e) {
      cat("預測錯誤:", e$message, "\n")
      return(FALSE)
    }
  )
}

cat("\n4. 執行測試預測...\n")

# 測試序列
test_sequences <- c(
  ">test_seq_1\nACDEFGHIKLMNPQRSTVWY",
  ">test_seq_2\nKLMNPQRSTVWY",
  ">test_seq_3\nACDEFGHIKL"
)

for (i in 1:length(test_sequences)) {
  cat(sprintf("測試序列 %d...\n", i))

  output_file <- sprintf("test_output_%d.txt", i)

  success <- simple_predict(test_sequences[i], output_file)

  if (success && file.exists(output_file)) {
    cat("✓ 預測成功\n")

    # 讀取結果
    results <- tryCatch(
      {
        read.table(output_file, header = FALSE, stringsAsFactors = FALSE, sep = " ")
      },
      error = function(e) {
        cat("讀取結果失敗:", e$message, "\n")
        return(NULL)
      }
    )

    if (!is.null(results) && nrow(results) > 0) {
      cat("結果:\n")
      for (j in 1:nrow(results)) {
        seq_name <- results[j, 1]
        prediction <- results[j, 2]
        probability <- results[j, 3]

        pred_text <- ifelse(prediction == 1, "AMP", ifelse(prediction == 0, "Non-AMP", "Error"))
        cat(sprintf("  %s: %s (機率: %s)\n", seq_name, pred_text, probability))
      }
    } else {
      cat("  結果文件為空或格式錯誤\n")
    }

    # 清理輸出文件
    unlink(output_file)
  } else {
    cat("✗ 預測失敗\n")
  }

  cat("\n")
}

cat("演示完成！\n")
