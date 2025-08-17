# AmPEP30 最終API - 基於100%驗證可工作的代碼

library(plumber)
library(jsonlite)

# 嚴格正規化與驗證 method
normalize_method <- function(method_raw) {
  if (is.null(method_raw) || length(method_raw) == 0) {
    return("rf")
  }
  m <- tolower(trimws(as.character(method_raw)))
  if (m %in% c("rf", "random_forest", "random-forest")) {
    return("rf")
  }
  if (m %in% c("cnn", "deep", "deep-ampep30", "deep_ampep30")) {
    return("cnn")
  }
  return(NA_character_)
}

# 處理 precision 的安全函數
normalize_precision <- function(p) {
  if (is.null(p) || is.na(suppressWarnings(as.numeric(p)))) {
    return(3L)
  }
  v <- as.integer(round(as.numeric(p)))
  if (v < 0) v <- 0L
  if (v > 6) v <- 6L
  v
}

# 直接使用已驗證的預測函數
# 嘗試多個可能的路徑
possible_demo_paths <- c(
  "AmPEP30_Demo.R",
  "../../AmPEP30_Demo.R",
  "/app/AmPEP30_Demo.R",
  "/Users/chonwai/Desktop/Self/AxPEP/Deep-AmPEP30/AmPEP30_Demo.R"
)

demo_loaded <- FALSE
for (path in possible_demo_paths) {
  if (file.exists(path)) {
    cat("載入演示函數從:", path, "\n")
    source(path)
    demo_loaded <- TRUE
    break
  }
}

if (!demo_loaded) {
  stop("無法找到 AmPEP30_Demo.R 文件")
}

#* Health check endpoint
#* @get /health
function() {
  list(
    status = "healthy",
    service = "AmPEP30-Final-API",
    version = "1.0.0",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

#* 預測單一序列
#* @post /predict/single
#* @param sequence:str 蛋白質序列
#* @param name:str 序列名稱 (可選)
#* @param method:str 預測方法，可選 rf 或 cnn（預設 rf）
#* @param precision:int 小數點位數（預設 3）
function(sequence, name = "query", method = "rf", precision = 3) {
  tryCatch(
    {
      # 驗證與正規化參數
      method_norm <- normalize_method(method)
      if (is.na(method_norm)) {
        return(list(
          status = "error",
          message = paste0("未知的 method: ", method, "，可用值: rf, cnn"),
          timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        ))
      }
      precision <- normalize_precision(precision)

      # 使用已驗證的預測函數
      result <- predict_peptide(sequence, name, method_norm)

      # 轉換為統一的API格式
      if (result$success) {
        pr <- as.numeric(result$amp_probability)
        api_result <- list(
          sequence_name = name,
          sequence = result$sequence,
          length = result$length,
          prediction = result$prediction,
          amp_probability = round(pr, precision),
          non_amp_probability = round(1 - pr, precision),
          confidence = round(as.numeric(result$confidence), precision),
          model_used = result$model_used,
          interpretation = result$interpretation,
          status = "success"
        )
      } else {
        # 統一錯誤格式 - 使用與成功響應相同的結構
        api_result <- list(
          sequence_name = name,
          sequence = sequence,
          length = if (nchar(sequence) > 0) nchar(sequence) else NA,
          prediction = NA,
          amp_probability = NA,
          non_amp_probability = NA,
          confidence = NA,
          model_used = method_norm,
          error = result$error,
          status = "error"
        )
      }

      return(api_result)
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
#* @param method:str 預測方法，可選 rf 或 cnn（預設 rf）
#* @param precision:int 小數點位數（預設 3）
function(fasta_content, method = "rf", precision = 3) {
  tryCatch(
    {
      # 驗證與正規化參數
      method_norm <- normalize_method(method)
      if (is.na(method_norm)) {
        return(list(
          status = "error",
          message = paste0("未知的 method: ", method, "，可用值: rf, cnn"),
          timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        ))
      }
      precision <- normalize_precision(precision)

      # 解析FASTA內容
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

      if (length(sequences) == 0) {
        return(list(
          status = "error",
          message = "未找到有效的FASTA序列",
          timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        ))
      }

      # 對每個序列進行預測
      results <- list()
      for (i in seq_along(sequences)) {
        seq_info <- sequences[[i]]
        result <- predict_peptide(seq_info$sequence, seq_info$name, method_norm)

        if (result$success) {
          pr <- as.numeric(result$amp_probability)
          api_result <- list(
            sequence_name = seq_info$name,
            sequence = result$sequence,
            length = result$length,
            prediction = result$prediction,
            amp_probability = round(pr, precision),
            non_amp_probability = round(1 - pr, precision),
            confidence = round(as.numeric(result$confidence), precision),
            model_used = result$model_used,
            interpretation = result$interpretation,
            status = "success"
          )
        } else {
          api_result <- list(
            sequence_name = seq_info$name,
            sequence = seq_info$sequence,
            error = result$error,
            status = "error"
          )
        }

        results[[i]] <- api_result
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
    model_type = "Random Forest / CNN (Keras)",
    model_file = "RF: AmPEP30-RF-1200tree.mdl; CNN: AmPEP30-CNN.mdl",
    features = "RF: PSEKRAAc (偽氨基酸組成); CNN: PSEKRAAc + CNN",
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
        result <- predict_peptide(seq_info$sequence, seq_info$name)

        if (result$success) {
          api_result <- list(
            sequence_name = seq_info$name,
            sequence = result$sequence,
            length = result$length,
            prediction = result$prediction,
            amp_probability = round(result$amp_probability, 3),
            non_amp_probability = round(1 - result$amp_probability, 3),
            confidence = round(result$confidence, 3),
            interpretation = result$interpretation,
            status = "success"
          )
        } else {
          api_result <- list(
            sequence_name = seq_info$name,
            sequence = seq_info$sequence,
            error = result$error,
            status = "error"
          )
        }

        results[[i]] <- api_result
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
