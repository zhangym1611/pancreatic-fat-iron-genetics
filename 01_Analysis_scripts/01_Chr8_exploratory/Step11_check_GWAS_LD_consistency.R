(function() {

  # ============================================================
  # 11修正版：GWAS-LD一致性诊断
  # 使用上一次已保存的频率QC后输入
  # ============================================================

  pkgs <- c("data.table", "susieR")
  ok <- vapply(
    pkgs, requireNamespace, logical(1), quietly = TRUE
  )

  if (!all(ok)) {
    stop("无法加载：", paste(pkgs[!ok], collapse = ", "))
  }

  base_dir <- path.expand("~/胰腺")

  input_file <- file.path(
    base_dir,
    "11_GWAS_LD_diagnostics",
    "chr8_diagnostic_c35bc64d238b",
    "LD_GWAS_after_frequency_QC.rds"
  )

  if (!file.exists(input_file)) {
    stop("未找到上次保存的输入：", input_file)
  }

  obj <- readRDS(input_file)

  traits <- c(
    "pancreas_fat", "pancreas_iron", "pancreas_volume"
  )

  if (!all(traits %in% names(obj$aligned))) {
    stop("输入中缺少三种表型。")
  }

  parent_dir <- file.path(
    base_dir, "11_GWAS_LD_diagnostics"
  )

  out_dir <- tempfile(
    "chr8_diagnostic_fixed_",
    tmpdir = parent_dir
  )

  if (!dir.create(out_dir)) {
    stop("无法创建输出目录。")
  }

  sink(file.path(out_dir, "run_log.txt"), split = TRUE)
  on.exit(sink(), add = TRUE)

  cat("结果目录：", out_dir, "\n")
  cat(
    "susieR版本：",
    as.character(utils::packageVersion("susieR")), "\n"
  )

  LD <- obj$LD
  aligned <- obj$aligned

  if (!is.matrix(LD) ||
      nrow(LD) != ncol(LD) ||
      is.null(colnames(LD)) ||
      !identical(rownames(LD), colnames(LD)) ||
      anyDuplicated(colnames(LD))) {
    stop("LD矩阵结构或变异名称异常。")
  }

  if ("rs200708125" %in% colnames(LD)) {
    stop("频率异常变异尚未排除，请核查输入。")
  }

  if (any(!is.finite(LD)) ||
      max(abs(LD - t(LD))) > 1e-8 ||
      max(abs(diag(LD) - 1)) > 1e-8) {
    stop("LD矩阵基本数值检查失败。")
  }

  for (trait in traits) {

    d <- aligned[[trait]]

    if (!identical(d$snp, colnames(LD))) {
      stop("GWAS与LD顺序不一致：", trait)
    }

    if (any(!is.finite(d$beta)) ||
        any(!is.finite(d$varbeta)) ||
        any(d$varbeta <= 0) ||
        length(d$N) != 1L ||
        !is.finite(d$N) || d$N <= 2) {
      stop("GWAS输入异常：", trait)
    }
  }

  key_snps <- c(
    "rs4733612", "rs1561927", "rs7824269",
    "rs7819154", "rs12676584"
  )

  stopifnot(all(key_snps %in% colnames(LD)))

  cat("保留变异数：", ncol(LD), "\n")

  # 缓存特征分解，后续如中途停止可复用本文件
  cat("\n正在计算LD特征分解……\n")
  flush.console()

  eig <- eigen(LD, symmetric = TRUE)

  cat("最小特征值：", min(eig$values), "\n")
  cat(
    "大于1e-8的特征值数量：",
    sum(eig$values > 1e-8), "\n"
  )

  if (min(eig$values) < -1e-8) {
    stop("LD存在超出容差的负特征值。")
  }

  saveRDS(
    eig,
    file.path(out_dir, "LD_eigendecomposition.rds")
  )

  attr(LD, "eigen") <- eig

  warning_log <- character()
  summary_rows <- list()

  save_plot <- function(plot_object, filename) {
    grDevices::pdf(
      filename, width = 7, height = 6,
      useDingbats = FALSE
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot_object)
  }

  # 保存警告的函数，保证中途停止前已有记录可保留
  save_warnings <- function() {
    writeLines(
      if (length(warning_log)) {
        warning_log
      } else {
        "No R warnings captured."
      },
      file.path(out_dir, "warnings.txt")
    )
  }

  for (trait in traits) {

    cat("\n===== 正在诊断：", trait, " =====\n")
    flush.console()

    d <- aligned[[trait]]
    z <- d$beta / sqrt(d$varbeta)

    result <- withCallingHandlers({

      s_est <- susieR::estimate_s_rss(
        z = z,
        R = LD,
        n = d$N,
        method = "null-mle"
      )

      if (length(s_est) != 1L ||
          !is.finite(s_est) ||
          s_est < 0 || s_est > 1) {
        stop("不一致性参数s异常：", trait)
      }

      cat("估计的不一致性参数s：", s_est, "\n")
      flush.console()

      if (s_est == 0) {
        list(
          s = s_est,
          kriging = NULL,
          status = "s_zero_boundary_local_diagnostic_not_run"
        )
      } else {
        list(
          s = s_est,
          kriging = susieR::kriging_rss(
            z = z,
            R = LD,
            n = d$N,
            s = s_est
          ),
          status = "completed"
        )
      }

    }, warning = function(w) {

      warning_log <<- c(
        warning_log,
        paste(trait, conditionMessage(w), sep = ": ")
      )

      invokeRestart("muffleWarning")
    })

    # 先保存完整对象，再解析返回结构
    saveRDS(
      result,
      file.path(out_dir, paste0(trait, "_diagnostic_raw.rds"))
    )

    save_warnings()

    n_flags <- NA_integer_

    if (!is.null(result$kriging)) {

      kr <- result$kriging

      cat(
        "kriging返回字段：",
        paste(names(kr), collapse = ", "), "\n"
      )

      # 首选官方源码字段；兼容可能存在的其他版本
      if (!is.null(kr$conditional_dist)) {
        tab <- as.data.frame(kr$conditional_dist)
      } else if (!is.null(kr$table)) {
        tab <- as.data.frame(kr$table)
      } else {
        stop(
          "未找到诊断表字段。完整结果已保存，",
          "请返回kriging返回字段。"
        )
      }

      if (nrow(tab) != length(z)) {
        stop(
          "诊断表有", nrow(tab),
          "行，但输入有", length(z),
          "个变异。完整对象已保存。"
        )
      }

      needed <- c(
        "z", "condmean", "condvar", "z_std_diff", "logLR"
      )

      if (!all(needed %in% names(tab))) {
        stop(
          "诊断表列名与预期不同：",
          paste(names(tab), collapse = ", "),
          "。完整对象已保存。"
        )
      }

      tab <- cbind(
        data.frame(
          SNP = d$snp,
          position = d$position,
          input_z = z,
          stringsAsFactors = FALSE
        ),
        tab
      )

      # 使用函数文档中的诊断标记规则
      # 仅标记待复核，不自动判定错误或删除
      tab$review_flag <- (
        is.finite(tab$logLR) &
        is.finite(tab$z) &
        tab$logLR > 2 &
        abs(tab$z) > 2
      )

      n_flags <- sum(tab$review_flag)

      data.table::fwrite(
        tab,
        file.path(out_dir, paste0(trait, "_kriging.tsv")),
        sep = "\t"
      )

      flagged <- tab[tab$review_flag, , drop = FALSE]

      if (nrow(flagged)) {
        flagged <- flagged[
          order(flagged$logLR, decreasing = TRUE),
          , drop = FALSE
        ]
      }

      data.table::fwrite(
        flagged,
        file.path(out_dir, paste0(trait, "_flagged_variants.tsv")),
        sep = "\t"
      )

      key_tab <- tab[
        match(key_snps, tab$SNP),
        , drop = FALSE
      ]

      data.table::fwrite(
        key_tab,
        file.path(
          out_dir,
          paste0(trait, "_five_candidate_diagnostics.tsv")
        ),
        sep = "\t"
      )

      if (!is.null(kr$plot)) {
        save_plot(
          kr$plot,
          file.path(out_dir, paste0(trait, "_kriging.pdf"))
        )
      }

      cat("\n五个关键候选诊断：\n")
      print(key_tab, row.names = FALSE, digits = 6)

      cat("\n需要复核的标记变异数：", n_flags, "\n")

      if (n_flags > 0L) {
        cat("按logLR排序的前10条标记记录：\n")
        print(
          head(flagged, 10L),
          row.names = FALSE,
          digits = 6
        )
      }
    }

    summary_rows[[trait]] <- data.frame(
      trait = trait,
      GWAS_N = d$N,
      reference_N = obj$reference_n,
      n_variants = length(z),
      estimated_s = result$s,
      n_review_flags = n_flags,
      status = result$status,
      stringsAsFactors = FALSE
    )

    # 每完成一个表型就保存一次汇总
    summary_table <- do.call(rbind, summary_rows)
    rownames(summary_table) <- NULL

    data.table::fwrite(
      summary_table,
      file.path(out_dir, "GWAS_LD_diagnostic_summary.tsv"),
      sep = "\t"
    )
  }

  save_warnings()

  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )

  cat("\n===== GWAS-LD一致性诊断汇总 =====\n")
  print(summary_table, row.names = FALSE, digits = 6)

  cat("\n诊断标记不等于已确认的等位基因错误。\n")
  cat("未自动翻转或删除其他变异，尚未运行SuSiE。\n")
  cat("结果目录：", out_dir, "\n")

})()