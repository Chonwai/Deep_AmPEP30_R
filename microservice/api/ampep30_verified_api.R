# AmPEP30 驗證API - 基於已驗證的本地預測功能
# 直接整合 RF-AmPEP30.R 和 AmPEP30_Demo.R 的核心功能

library(plumber)
library(jsonlite)

# 載入必要套件
suppressMessages({
  library(randomForest)
  library(seqinr)
  library(protr)
})

# 全域變量
RF_MODEL <- NULL
MODEL_LOADED <- FALSE

# 載入RF模型
load_rf_model <- function() {
  if (!MODEL_LOADED) {
    # 嘗試多個可能的模型路徑
    possible_paths <- c(
      "AmPEP30-RF-1200tree.mdl",
      "/app/AmPEP30-RF-1200tree.mdl",
      file.path(getwd(), "AmPEP30-RF-1200tree.mdl"),
      "../../AmPEP30-RF-1200tree.mdl",
      "/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30/AmPEP30-RF-1200tree.mdl"
    )

    model_path <- NULL
    for (path in possible_paths) {
      if (file.exists(path)) {
        model_path <- path
        break
      }
    }

    if (is.null(model_path)) {
      stop(paste("模型文件不存在於任何路徑:", paste(possible_paths, collapse = ", ")))
    }

    cat("正在載入RF模型...\n")
    RF_MODEL <<- readRDS(model_path)
    MODEL_LOADED <<- TRUE
    cat("✓ RF模型載入成功\n")
  }
  return(RF_MODEL)
}

# 核心預測函數 - 直接來自我們驗證的代碼
predict_peptide_core <- function(sequence, sequence_name = "query") {
  # 驗證序列
  sequence <- toupper(trimws(sequence))
  seq_length <- nchar(sequence)

  # 檢查長度
  if (seq_length < 5 || seq_length > 30) {
    stop(paste("序列長度必須在5-30之間，當前長度:", seq_length))
  }

  # 檢查氨基酸
  valid_aa <- "ACDEFGHIKLMNPQRSTVWY"
  invalid_chars <- setdiff(strsplit(sequence, "")[[1]], strsplit(valid_aa, "")[[1]])
  if (length(invalid_chars) > 0) {
    stop(paste("包含無效氨基酸:", paste(invalid_chars, collapse = ", ")))
  }

  # 載入模型
  model <- load_rf_model()

  # 提取特徵 (PSEKRAAc)
  tryCatch(
    {
      feature_vector <- psekraac(sequence, 3, .4, 6)
      prediction <- predict(model, feature_vector, type = "prob")

      # 提取機率
      amp_prob <- as.numeric(prediction[1, "AMP"])
      non_amp_prob <- as.numeric(prediction[1, "non-AMP"])

      # 預測類別
      predicted_class <- ifelse(amp_prob > 0.5, "AMP", "Non-AMP")
      confidence <- max(amp_prob, non_amp_prob)

      # 解釋
      if (predicted_class == "AMP") {
        interpretation <- sprintf("此序列很可能是抗菌胜肽 (機率: %.1f%%)", amp_prob * 100)
      } else {
        interpretation <- sprintf("此序列不太可能是抗菌胜肽 (非 AMP 機率: %.1f%%)", non_amp_prob * 100)
      }

      return(list(
        sequence_name = sequence_name,
        sequence = sequence,
        length = seq_length,
        prediction = predicted_class,
        amp_probability = round(amp_prob, 3),
        non_amp_probability = round(non_amp_prob, 3),
        confidence = round(confidence, 3),
        interpretation = interpretation,
        status = "success"
      ))
    },
    error = function(e) {
      return(list(
        sequence_name = sequence_name,
        sequence = sequence,
        length = seq_length,
        error = e$message,
        status = "error"
      ))
    }
  )
}

# 處理FASTA格式的函數
parse_fasta_content <- function(fasta_content) {
  lines <- strsplit(fasta_content, "\n")[[1]]
  lines <- lines[lines != ""]

  sequences <- list()
  current_name <- NULL
  current_seq <- ""

  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, ">")) {
      # 保存前一個序列
      if (!is.null(current_name) && current_seq != "") {
        sequences[[length(sequences) + 1]] <- list(
          name = current_name,
          sequence = current_seq
        )
      }
      # 開始新序列
      current_name <- substr(line, 2, nchar(line))
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }

  # 保存最後一個序列
  if (!is.null(current_name) && current_seq != "") {
    sequences[[length(sequences) + 1]] <- list(
      name = current_name,
      sequence = current_seq
    )
  }

  return(sequences)
}

#* Health check endpoint
#* @get /health
function() {
  list(
    status = "healthy",
    service = "AmPEP30-Verified-API",
    version = "1.0.0",
    model_loaded = MODEL_LOADED,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

#* 預測單一序列
#* @post /predict/single
#* @param sequence:str 蛋白質序列
#* @param name:str 序列名稱 (可選)
function(sequence, name = "query") {
  tryCatch(
    {
      result <- predict_peptide_core(sequence, name)
      return(result)
    },
    error = function(e) {
      return(list(
        status = "error",
        message = e$message,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    }
  )
}

#* 預測FASTA格式的多個序列
#* @post /predict/fasta
#* @param fasta_content:str FASTA格式的序列內容
function(fasta_content) {
  tryCatch(
    {
      sequences <- parse_fasta_content(fasta_content)

      if (length(sequences) == 0) {
        return(list(
          status = "error",
          message = "未找到有效的FASTA序列",
          timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        ))
      }

      results <- list()
      for (i in seq_along(sequences)) {
        seq_info <- sequences[[i]]
        result <- predict_peptide_core(seq_info$sequence, seq_info$name)
        results[[i]] <- result
      }

      return(list(
        status = "success",
        total_sequences = length(sequences),
        results = results,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    },
    error = function(e) {
      return(list(
        status = "error",
        message = e$message,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    }
  )
}

#* 獲取模型資訊
#* @get /model/info
function() {
  list(
    model_type = "Random Forest",
    model_file = "AmPEP30-RF-1200tree.mdl",
    model_loaded = MODEL_LOADED,
    features = "PSEKRAAc (偽氨基酸組成)",
    sequence_length_range = "5-30 amino acids",
    supported_amino_acids = "ACDEFGHIKLMNPQRSTVWY",
    training_data_size = "3298 sequences",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

#* 測試端點 - 使用已知的測試序列
#* @get /test/demo
function() {
  # 使用我們已經驗證的測試序列
  test_sequences <- list(
    list(name = "Magainin-2", sequence = "GIGKFLHSAKKFGKAFVGEIMNS"),
    list(name = "Melittin", sequence = "GIGAVLKVLTTGLPALISWIKRKRQQ"),
    list(name = "Random_peptide", sequence = "ACDEFGHIKL")
  )

  tryCatch(
    {
      results <- list()
      for (i in seq_along(test_sequences)) {
        seq_info <- test_sequences[[i]]
        result <- predict_peptide_core(seq_info$sequence, seq_info$name)
        results[[i]] <- result
      }

      return(list(
        status = "success",
        message = "Demo test completed successfully",
        total_sequences = length(test_sequences),
        results = results,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    },
    error = function(e) {
      return(list(
        status = "error",
        message = paste("Demo test failed:", e$message),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ))
    }
  )
}
