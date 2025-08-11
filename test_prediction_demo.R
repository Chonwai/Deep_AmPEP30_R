#!/usr/bin/env Rscript

# 直接演示預測功能
cat("=== AmPEP30 預測功能演示 ===\n\n")

# 1. RF 預測演示
cat("1. 隨機森林 (RF) 預測演示\n")
cat("載入 RF 腳本...\n")

tryCatch(
  {
    source("RF-AmPEP30.R")
    cat("✓ RF 腳本載入成功\n")

    # 創建測試序列
    test_fasta <- ">test_sequence_1
ACDEFGHIKLMNPQRSTVWY
>test_sequence_2
GHRKLMNPQRSTVWY
>test_sequence_3
ACDEFG"

    cat("測試序列:\n", test_fasta, "\n\n")

    # 創建臨時文件
    temp_input <- "temp_test.fasta"
    temp_output <- "temp_test.out"

    writeLines(test_fasta, temp_input)

    cat("執行 RF 預測...\n")
    predict_from_fasta_rf(temp_input, temp_output)

    if (file.exists(temp_output)) {
      cat("✓ RF 預測完成\n")
      results <- read.table(temp_output, header = FALSE, stringsAsFactors = FALSE)
      cat("預測結果:\n")
      for (i in 1:nrow(results)) {
        cat(sprintf(
          "序列 %s: 預測=%s, 機率=%.3f\n",
          results[i, 1], results[i, 2], as.numeric(results[i, 3])
        ))
      }
    } else {
      cat("✗ 預測輸出文件未生成\n")
    }

    # 清理
    unlink(c(temp_input, temp_output))
  },
  error = function(e) {
    cat("✗ RF 預測錯誤:", e$message, "\n")
  }
)

cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")

# 2. Deep Learning 預測演示
cat("2. 深度學習 (Deep Learning) 預測演示\n")
cat("載入 Deep Learning 腳本...\n")

tryCatch(
  {
    source("Deep-AmPEP30.R")
    cat("✓ Deep Learning 腳本載入成功\n")

    # 創建測試序列
    test_fasta2 <- ">test_deep_1
ACDEFGHIKLMNPQRSTVWY
>test_deep_2
KLMNPQRSTVWY"

    cat("測試序列:\n", test_fasta2, "\n\n")

    # 創建臨時文件
    temp_input2 <- "temp_deep.fasta"
    temp_output2 <- "temp_deep.out"

    writeLines(test_fasta2, temp_input2)

    cat("執行深度學習預測...\n")
    predict_from_fasta_deep(temp_input2, temp_output2)

    if (file.exists(temp_output2)) {
      cat("✓ 深度學習預測完成\n")
      results2 <- read.table(temp_output2, header = FALSE, stringsAsFactors = FALSE)
      cat("預測結果:\n")
      for (i in 1:nrow(results2)) {
        cat(sprintf(
          "序列 %s: 預測=%s, 機率=%.3f\n",
          results2[i, 1], results2[i, 2], as.numeric(results2[i, 3])
        ))
      }
    } else {
      cat("✗ 深度學習預測輸出文件未生成\n")
    }

    # 清理
    unlink(c(temp_input2, temp_output2))
  },
  error = function(e) {
    cat("✗ Deep Learning 預測錯誤:", e$message, "\n")
  }
)

cat("\n演示完成！\n")
