#!/usr/bin/env Rscript

# ========================================
# AmPEP30 抗菌胜肽預測系統演示
# ========================================
#
# 此演示展示如何使用 AmPEP30 系統預測蛋白質序列是否為抗菌胜肽 (AMP)
# 系統基於隨機森林模型，支援 5-30 氨基酸長度的序列
#

cat("
========================================
    AmPEP30 抗菌胜肽預測系統演示
========================================

功能：預測蛋白質序列是否為抗菌胜肽 (AMP)
模型：隨機森林 (Random Forest)
序列長度：5-30 氨基酸
支援氨基酸：ACDEFGHIKLMNPQRSTVWY

")

# 設定工作環境（容器與本機皆可用）
app_root <- Sys.getenv("APP_ROOT", unset = "/app")
if (dir.exists(app_root)) {
  try(setwd(app_root), silent = TRUE)
}

# 載入必要套件
cat("正在載入必要套件...\n")
suppressWarnings(suppressMessages({
  library(seqinr)
  library(caret)
  library(randomForest)
  library(protr)
}))
cat("✓ 套件載入完成\n\n")

# 預測函數
predict_peptide <- function(sequence, sequence_name = "query", method = "rf") {
  # 預測蛋白質序列是否為 AMP
  #
  # 參數:
  #   sequence: 蛋白質序列 (字符串)
  #   sequence_name: 序列名稱 (可選)
  #
  # 返回:
  #   列表包含預測結果和機率

  # 驗證序列
  sequence <- toupper(trimws(sequence))
  seq_length <- nchar(sequence)

  if (seq_length < 5 || seq_length > 30) {
    return(list(
      success = FALSE,
      error = paste0("序列長度必須在 5-30 氨基酸之間，當前長度: ", seq_length)
    ))
  }

  # 檢查氨基酸字符
  valid_aa <- "ACDEFGHIKLMNPQRSTVWY"
  invalid_chars <- setdiff(strsplit(sequence, "")[[1]], strsplit(valid_aa, "")[[1]])
  if (length(invalid_chars) > 0) {
    return(list(
      success = FALSE,
      error = paste0("序列包含無效字符: ", paste(invalid_chars, collapse = ", "))
    ))
  }

  tryCatch(
    {
      # 生成臨時文件名（使用可寫的臨時目錄，並使用腳本的絕對路徑）
      app_root <- Sys.getenv("APP_ROOT", unset = "/app")

      tmp_dir <- tempdir()
      temp_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1))
      temp_fasta <- file.path(tmp_dir, paste0("temp_", temp_id, ".fasta"))
      temp_output <- file.path(tmp_dir, paste0("temp_", temp_id, ".out"))

      # 創建 FASTA 文件
      fasta_content <- paste0(">", sequence_name, "\n", sequence)
      writeLines(fasta_content, temp_fasta)

      # Rscript 可執行檔
      rscript_bin <- unname(Sys.which("Rscript"))
      if (is.na(rscript_bin) || rscript_bin == "") rscript_bin <- "Rscript"

      method_norm <- tolower(trimws(method))
      if (method_norm %in% c("rf", "random_forest", "random-forest")) {
        # RF 流程
        rf_candidates <- c(
          file.path(getwd(), "RF-AmPEP30.R"),
          file.path(app_root, "RF-AmPEP30.R"),
          "RF-AmPEP30.R"
        )
        rf_existing <- rf_candidates[file.exists(rf_candidates)]
        if (length(rf_existing) == 0) {
          unlink(c(temp_fasta, temp_output), force = TRUE)
          return(list(success = FALSE, error = "找不到 RF 預測腳本 RF-AmPEP30.R"))
        }
        script_path <- rf_existing[[1]]
        model_used <- "rf"
      } else if (method_norm %in% c("cnn", "deep", "deep-ampep30", "deep_ampep30")) {
        # CNN 流程
        cnn_candidates <- c(
          file.path(getwd(), "Deep-AmPEP30.R"),
          file.path(app_root, "Deep-AmPEP30.R"),
          "Deep-AmPEP30.R"
        )
        cnn_existing <- cnn_candidates[file.exists(cnn_candidates)]
        if (length(cnn_existing) == 0) {
          unlink(c(temp_fasta, temp_output), force = TRUE)
          return(list(success = FALSE, error = "找不到 CNN 預測腳本 Deep-AmPEP30.R"))
        }
        script_path <- cnn_existing[[1]]
        model_used <- "cnn"
      } else {
        unlink(c(temp_fasta, temp_output), force = TRUE)
        return(list(success = FALSE, error = paste0("未知的 method: ", method, "，可用值: rf, cnn")))
      }

      cmd <- paste(shQuote(rscript_bin), shQuote(script_path), shQuote(temp_fasta), shQuote(temp_output))
      system_result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

      # 檢查執行結果
      if (system_result != 0 || !file.exists(temp_output)) {
        unlink(c(temp_fasta, temp_output), force = TRUE)
        return(list(success = FALSE, error = "預測執行失敗"))
      }

      # 讀取結果
      results <- read.table(temp_output, header = FALSE, stringsAsFactors = FALSE, sep = " ")
      unlink(c(temp_fasta, temp_output), force = TRUE)

      if (nrow(results) > 0) {
        prediction <- results[1, 2]
        probability <- as.numeric(results[1, 3])

        return(list(
          success = TRUE,
          sequence = sequence,
          length = seq_length,
          prediction = prediction, # 保持數字格式：1=AMP, 0=Non-AMP, -1=Invalid
          amp_probability = probability,
          confidence = ifelse(prediction == 1, probability, ifelse(prediction == 0, 1 - probability, -1)),
          model_used = model_used,
          interpretation = ifelse(
            prediction == 1,
            paste0("此序列很可能是抗菌胜肽 (機率: ", sprintf("%.1f%%", probability * 100), ")"),
            ifelse(
              prediction == 0,
              paste0("此序列不太可能是抗菌胜肽 (非 AMP 機率: ", sprintf("%.1f%%", (1 - probability) * 100), ")"),
              "序列長度不符合要求 (5-30 個氨基酸)"
            )
          )
        ))
      } else {
        return(list(success = FALSE, error = "無法解析預測結果"))
      }
    },
    error = function(e) {
      return(list(success = FALSE, error = paste("預測錯誤:", e$message)))
    }
  )
}

# 演示案例
cat("=== 演示案例 ===\n\n")

demo_sequences <- list(
  list(
    name = "Magainin-2",
    sequence = "GIGKFLHSAKKFGKAFVGEIMNS",
    description = "已知的青蛙抗菌胜肽"
  ),
  list(
    name = "Melittin",
    sequence = "GIGAVLKVLTTGLPALISWIKRKRQQ",
    description = "蜜蜂毒液中的抗菌胜肽"
  ),
  list(
    name = "Random_peptide",
    sequence = "ACDEFGHIKL",
    description = "隨機氨基酸序列"
  ),
  list(
    name = "Short_peptide",
    sequence = "KLMNP",
    description = "短肽序列"
  )
)

for (i in 1:length(demo_sequences)) {
  demo <- demo_sequences[[i]]

  cat(sprintf("%d. %s\n", i, demo$name))
  cat(sprintf("   描述: %s\n", demo$description))
  cat(sprintf("   序列: %s\n", demo$sequence))
  cat(sprintf("   長度: %d 氨基酸\n", nchar(demo$sequence)))

  result <- predict_peptide(demo$sequence, demo$name)

  if (result$success) {
    cat(sprintf("   預測: %s\n", result$prediction))
    cat(sprintf("   AMP 機率: %.3f (%.1f%%)\n", result$amp_probability, result$amp_probability * 100))
    cat(sprintf("   信心度: %.3f\n", result$confidence))
    cat(sprintf("   解釋: %s\n", result$interpretation))
  } else {
    cat(sprintf("   錯誤: %s\n", result$error))
  }

  cat("\n")
}

cat("=== 使用指南 ===\n\n")
cat("1. 序列要求：\n")
cat("   - 長度：5-30 氨基酸\n")
cat("   - 字符：標準 20 種氨基酸 (ACDEFGHIKLMNPQRSTVWY)\n")
cat("   - 格式：大小寫不敏感，會自動轉換為大寫\n\n")

cat("2. 結果解釋：\n")
cat("   - AMP 機率 > 0.5：預測為抗菌胜肽\n")
cat("   - AMP 機率 < 0.5：預測為非抗菌胜肽\n")
cat("   - 信心度：預測的可信程度\n\n")

cat("3. 模型性能：\n")
cat("   - 基於隨機森林演算法\n")
cat("   - 訓練數據：3298 個序列\n")
cat("   - 特徵：PSEKRAAc (偽氨基酸組成)\n\n")

cat("演示完成！\n")
cat("如需預測自己的序列，請調用 predict_peptide(your_sequence) 函數。\n")
